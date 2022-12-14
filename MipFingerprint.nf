#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Utils modules
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'
include ExportParams as Workflow_ExportParams from './NextflowModules/Utils/workflow.nf'

// Mapping modules
include MEM as BWA_MEM from './NextflowModules/BWA/0.7.17/MEM.nf' params(genome:"$params.genome", optional: '-c 100 -M')
include ViewSort as Sambamba_ViewSort from './NextflowModules/Sambamba/0.7.0/ViewSort.nf'

// Fingerprint modules
include UnifiedGenotyper as GATK_UnifiedGenotyper from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/UnifiedGenotyper.nf' params(gatk_path: "$params.gatk_path", genome:"$params.genome", optional: "--intervals $params.dxtracks_path/$params.fingerprint_target --output_mode EMIT_ALL_SITES")

// QC Modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include Flagstat as Sambamba_Flagstat from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include MultiQC from './NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:"--config $baseDir/assets/multiqc_config.yaml")

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def samples = fastq_files.map({it.flatten()}).groupTuple(by:[0], sort:true)
def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Trim and Dedup MIP fastq files
    MipsTrimDedup(samples)

    //Mapping
    BWA_MEM(
        MipsTrimDedup.out.map{sample_id, rg_id, r1_fastq, r2_fastq -> [sample_id, rg_id, [r1_fastq, r2_fastq]]}
    )
    Sambamba_ViewSort(BWA_MEM.out)

    // Fingerprint
    GATK_UnifiedGenotyper(Sambamba_ViewSort.out.map{sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())
    CheckFingerprintVCF(GATK_UnifiedGenotyper.out.map{sample_id, vcf_file -> [vcf_file]}.collect())

    // QC
    FastQC(fastq_files)
    Sambamba_Flagstat(Sambamba_ViewSort.out.map{sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple())
    MultiQC(analysis_id, Channel.empty().mix(FastQC.out.collect()))

    // Create log files: Repository versions and Workflow params
    VersionLog()
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "MIP Fingerprint Workflow Successful: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html")
    } else {
        def subject = "MIP Fingerprint Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}

// Custom processes
process MipsTrimDedup {
    // Custom process to run MIPS TrimDedup
    tag {"MIPS TrimDedup ${sample_id} - ${rg_id}"}
    label 'MIPS_1_0_1'
    label 'MIPS_1_0_1_TrimDedup'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(sample_id, rg_id, path(r1_fastqs), path(r2_fastqs))

    output:
        tuple(sample_id, rg_id, path('*_LMergedTrimmedDedup_R1_*.fastq.gz'), path('*_LMergedTrimmedDedup_R2_*.fastq.gz'), emit: fastq_files)

    script:
        def r1_args = r1_fastqs.collect{ "$it" }.join(" ")
        def r2_args = r2_fastqs.collect{ "$it" }.join(" ")

        rg_id = "${sample_id}_MergedTrimmedDedup"

        """
        python2 ${params.mips_trim_dedup_path}/mips_trim_dedup.py -d ${params.dxtracks_path}/${params.mips_design_file}  -l ${params.mips_uuid_length} -ur ${params.mips_uuid_read} -r1 ${r1_args} -r2 ${r2_args}
        """
}

process CheckFingerprintVCF {
    // Custom process to check fingerprint vcf files
    tag {"CheckFingerprintVCF"}
    label 'CheckFingerprintVCF'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(vcf_files)

    output:
        tuple(path('disapprovedVCFs'), path('approvedVCFs/*.vcf'), emit: vcf_files)
        path('logbook.txt', emit: logbook)


    script:
        """
        python2 ${baseDir}/assets/check_fingerprint_vcf.py ${vcf_files} > logbook.txt
        """
}

process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log', emit: log_file)

    script:
        """
        echo 'DxNextflowMIP' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'Dx_tracks' >> repository_version.log
        git --git-dir=${params.dxtracks_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'MipsTrimDedup' >> repository_version.log
        git --git-dir=${params.mips_trim_dedup_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}
