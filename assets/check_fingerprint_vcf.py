import argparse
import os
import shutil


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check fingerprint vcf files.')
    parser.add_argument('fingerprint_vcf_files', type=argparse.FileType('r'), nargs='*', help='Fingerprint VCF')
    arguments = parser.parse_args()

    badfolder = 'disapprovedVCFs'
    os.mkdir(badfolder)

    gender = ['M', 'F', 'O']

    logbook = []

    for vcf_file in arguments.fingerprint_vcf_files:
        refcov = 0  # Total reference reads for all homozygous (1/1) calls
        totalcov = 0  # Total coverage for all homozygous calls
        homaltcount = 0  # Number of homozygous calls
        ycount = 0  # Sum of coverage for two Y SNPs
        lowcovcount = 0  # Number of SNPs with <15X coverage
        disbalancecount = 0  # Number of heterozygous (0/1) calls with allelefrequency <0.2 or >0.8

        for line in vcf_file:
            if not line.startswith('#'):
                line = line.split()

                # Parse Genotype format
                gt_format = line[8].split(':')
                gt_index = gt_format.index('GT')

                # Parse sample genotype
                gt_values = line[9].split(':')
                gt_value = gt_values[gt_index]

                if line[0] == 'Y':
                    if gt_value != './.':
                        ycount += int(gt_values[gt_format.index('DP')])
                elif gt_value == '1/1':
                    homaltcount += 1
                    if int(gt_values[gt_format.index('DP')]) < 15:
                        lowcovcount += 1
                    refcov += int(gt_values[gt_format.index('AD')].split(',')[0])
                    totalcov += int(gt_values[gt_format.index('DP')])
                elif gt_value != '1/1':
                    if gt_value == './.':
                        lowcovcount += 1
                    elif gt_value == '0/0':
                        if int(gt_values[gt_format.index('DP')]) < 15:
                            lowcovcount += 1
                    else:
                        if int(gt_values[gt_format.index('DP')]) < 15:
                            lowcovcount += 1
                if gt_value == '0/1':
                    af_value = int(gt_values[gt_format.index('AD')].split(',')[0]) / float(int(gt_values[gt_format.index('DP')]))
                    if af_value > 0.8 or af_value < 0.2:
                        disbalancecount += 1

        contamination = refcov / float(totalcov)

        result = vcf_file.name, str(lowcovcount), str(homaltcount), str(round(contamination, 6)), vcf_file.name[8], str(ycount), str(disbalancecount)
        logbook.append(result)

        if result[4] not in gender:
            print('### {}: report filename issue to lab'.format(vcf_file.name))
        if int(result[1]) > 15:
            print('### {}: >15 SNPs with <15X coverage ({}) --> disapproved'.format(vcf_file.name, result[1]))
            shutil.move(vcf_file.name, badfolder)
        elif int(result[6]) > 8:
            print('### {}: >8 heterozygous SNPs with <20% MAF ({}) --> disapproved'.format(vcf_file.name, result[6]))
            shutil.move(vcf_file.name, badfolder)
        elif int(result[2]) < 8:
            print('### {}: <8 homozygous ALT SNPs called ({}) --> disapproved'.format(vcf_file.name, result[2]))
            shutil.move(vcf_file.name, badfolder)
        elif result[4] == 'F' and int(result[5]) > 100 or result[4] == 'M' and int(result[5]) < 100:
            print('### {}: gender {} with {} reads on chromosome Y, discuss with lab and disapprove if needed'.format(vcf_file.name, result[4], result[5]))

    for line in logbook:
        print '\t'.join(line)
