#!/bin/sh
# submits a job to filter out cosmopolitan variants from a vcf filtered from 1kg
# USAGE: cosmo_multiallelic_v2.sh CHROM

HAND="chr${1}_cosmo"

echo "Filtering cosmopolitan variants from chromosome ${1}." 

cmd="vcf-isec -n +2 -f ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_AFR.vcf.gz ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_EAS.vcf.gz ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_SAS.vcf.gz ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_EUR.vcf.gz | bgzip -c > chr${1}_COSMO_private_multi.vcf.gz"

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
