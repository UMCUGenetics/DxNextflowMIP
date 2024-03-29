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

// CustomModules
include CheckFingerprintVCF from './CustomModules/CheckFingerprintVCF/CheckFingerprintVCF.nf'
include CheckQC from './CustomModules/CheckQC/CheckQC.nf'
include MipsTrimDedup from './CustomModules/MipsTrimDedup/MipsTrimDedup.nf'
include VersionLog from './CustomModules/Utils/VersionLog.nf'

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

    // QC - Collect and check
    CheckQC(
        analysis_id, 
        CheckFingerprintVCF.out.logbook
    )

    // Create log files: Repository versions and Workflow params
    VersionLog(
        Channel.of(
            "${workflow.projectDir}/",
            "${params.dxtracks_path}/",
            "${params.mips_trim_dedup_path}/",
        ).collect()
    )
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
