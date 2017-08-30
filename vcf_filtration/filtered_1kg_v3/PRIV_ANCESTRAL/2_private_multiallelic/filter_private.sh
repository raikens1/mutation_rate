#!/bin/sh
# submits a job to filter out private variants from a vcf filtered from 1kg
# USAGE: filter_private.sh POP CHROM

HAND="chr${2}_${1}_priv"

echo "Filtering ${1} variants from chromosome ${2}." 

if [ ${1} == "AFR" ]; then
	cmd="vcf-isec -c ../1_filtered_vcfs/chr${2}_AFR.vcf.gz ../1_filtered_vcfs/chr${2}_EAS.vcf.gz ../1_filtered_vcfs/chr${2}_SAS.vcf.gz ../1_filtered_vcfs/chr${2}_EUR.vcf.gz -f | bgzip -c > chr${2}_AFR_private_multi.vcf.gz"
elif [ ${1} == "EUR" ]; then
	cmd="vcf-isec -c ../1_filtered_vcfs/chr${2}_EUR.vcf.gz ../1_filtered_vcfs/chr${2}_EAS.vcf.gz ../1_filtered_vcfs/chr${2}_SAS.vcf.gz ../1_filtered_vcfs/chr${2}_AFR.vcf.gz -f | bgzip -c > chr${2}_EUR_private_multi.vcf.gz"
elif [ ${1} == "EAS" ]; then
	cmd="vcf-isec -c ../1_filtered_vcfs/chr${2}_EAS.vcf.gz ../1_filtered_vcfs/chr${2}_AFR.vcf.gz ../1_filtered_vcfs/chr${2}_SAS.vcf.gz ../1_filtered_vcfs/chr${2}_EUR.vcf.gz -f | bgzip -c > chr${2}_EAS_private_multi.vcf.gz"
elif [ ${1} == "SAS" ]; then
	cmd="vcf-isec -c ../1_filtered_vcfs/chr${2}_SAS.vcf.gz ../1_filtered_vcfs/chr${2}_EAS.vcf.gz ../1_filtered_vcfs/chr${2}_AFR.vcf.gz ../1_filtered_vcfs/chr${2}_EUR.vcf.gz -f | bgzip -c > chr${2}_SAS_private_multi.vcf.gz"
else
	echo "Error: not a valid population code."
	exit 1
fi

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"

