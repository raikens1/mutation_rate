#!/bin/sh
# submits a job to filter out private variants from a vcf filtered from 1kg
# USAGE: SGDP_intersect.sh POP CHROM

HAND="extract_private_chr_${2}_SGDP_${1}"

echo "Filtering ${1} variants from chromosome ${2}."

if [ ${1} == "AFR" ]; then
	cmd="module load bcftools; bcftools isec -n=1 AFR/chr_${2}_SGDP_AFR_all_snps.vcf.gz EAS/chr_${2}_SGDP_EAS_all_snps.vcf.gz SAS/chr_${2}_SGDP_SAS_all_snps.vcf.gz EUR/chr_${2}_SGDP_EUR_all_snps.vcf.gz -w1 -Oz -o chr_${2}_SGDP_AFR_private.vcf.gz"
elif [ ${1} == "EUR" ]; then
	cmd="module load bcftools; bcftools isec -n=1 AFR/chr_${2}_SGDP_AFR_all_snps.vcf.gz EAS/chr_${2}_SGDP_EAS_all_snps.vcf.gz SAS/chr_${2}_SGDP_SAS_all_snps.vcf.gz EUR/chr_${2}_SGDP_EUR_all_snps.vcf.gz -w4 -Oz -o chr_${2}_SGDP_EUR_private.vcf.gz"
elif [ ${1} == "EAS" ]; then
	cmd="module load bcftools; bcftools isec -n=1 AFR/chr_${2}_SGDP_AFR_all_snps.vcf.gz EAS/chr_${2}_SGDP_EAS_all_snps.vcf.gz SAS/chr_${2}_SGDP_SAS_all_snps.vcf.gz EUR/chr_${2}_SGDP_EUR_all_snps.vcf.gz -w2 -Oz -o chr_${2}_SGDP_EAS_private.vcf.gz"
elif [ ${1} == "SAS" ]; then
	cmd="module load bcftools; bcftools isec -n=1 AFR/chr_${2}_SGDP_AFR_all_snps.vcf.gz EAS/chr_${2}_SGDP_EAS_all_snps.vcf.gz SAS/chr_${2}_SGDP_SAS_all_snps.vcf.gz EUR/chr_${2}_SGDP_EUR_all_snps.vcf.gz -w3 -Oz -o chr_${2}_SGDP_SAS_private.vcf.gz"
else
	echo "Error: not a valid population code."
	exit 1
fi

echo ${cmd}
bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
