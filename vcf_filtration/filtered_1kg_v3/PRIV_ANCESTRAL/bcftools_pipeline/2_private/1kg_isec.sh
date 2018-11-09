#!/bin/sh
# submits a job to filter out private variants from a vcf filtered from 1kg
# USAGE: 1kg_isec.sh POP CHROM
# NOTE: this script uses relative filepaths. It must be run in the subdirectory `2_private` to function properly

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "# submits a job to filter out private variants from a vcf filtered from 1kg
# USAGE: 1kg_isec.sh POP CHROM      
# NOTE: this script uses relative filepaths. It must be run in the subdirectory `2_private` to function properly"
    exit 1	
fi

HAND="extract_private_chr_${2}_1kg_${1}"

echo "Filtering ${1} variants from chromosome ${2}."

if [ ${1} == "AFR" ]; then
	cmd="module load bcftools; bcftools isec -n=1 ../1_filtered/chr_${2}_AFR_all_snps.vcf.gz ../1_filtered/chr_${2}_EAS_all_snps.vcf.gz ../1_filtered/chr_${2}_SAS_all_snps.vcf.gz ../1_filtered/chr_${2}_EUR_all_snps.vcf.gz -w1 -Oz -o chr_${2}_AFR_private.vcf.gz; tabix -p vcf chr_${2}_AFR_private.vcf.gz"
elif [ ${1} == "EUR" ]; then
	cmd="module load bcftools; bcftools isec -n=1 ../1_filtered/chr_${2}_AFR_all_snps.vcf.gz ../1_filtered/chr_${2}_EAS_all_snps.vcf.gz ../1_filtered/chr_${2}_SAS_all_snps.vcf.gz ../1_filtered/chr_${2}_EUR_all_snps.vcf.gz -w4 -Oz -o chr_${2}_EUR_private.vcf.gz; tabix -p vcf chr_${2}_EUR_private.vcf.gz"
elif [ ${1} == "EAS" ]; then
	cmd="module load bcftools; bcftools isec -n=1 ../1_filtered/chr_${2}_AFR_all_snps.vcf.gz ../1_filtered/chr_${2}_EAS_all_snps.vcf.gz ../1_filtered/chr_${2}_SAS_all_snps.vcf.gz ../1_filtered/chr_${2}_EUR_all_snps.vcf.gz -w2 -Oz -o chr_${2}_EAS_private.vcf.gz; tabix -p vcf chr_${2}_EAS_private.vcf.gz"
elif [ ${1} == "SAS" ]; then
	cmd="module load bcftools; bcftools isec -n=1 ../1_filtered/chr_${2}_AFR_all_snps.vcf.gz ../1_filtered/chr_${2}_EAS_all_snps.vcf.gz ../1_filtered/chr_${2}_SAS_all_snps.vcf.gz ../1_filtered/chr_${2}_EUR_all_snps.vcf.gz -w3 -Oz -o chr_${2}_SAS_private.vcf.gz; tabix -p vcf chr_${2}_SAS_private.vcf.gz"
else
	echo "Error: not a valid population code."
	exit 1
fi

echo ${cmd}
bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
