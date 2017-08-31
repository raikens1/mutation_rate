#!/bin/sh
# submits a job to filter out variants shared by all ancestral continental groups from vcfs filtered from 1kg
# USAGE: ALL_multiallelic.sh CHROM

HAND="chr${1}_ALL"

echo "Filtering variants shared by all populations on chromosome ${1}." 

cmd="vcf-isec -n +4 -f ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_AFR.vcf.gz ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_EAS.vcf.gz ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_SAS.vcf.gz ../../PRIV_ANCESTRAL/1_filtered_vcfs/chr${1}_EUR.vcf.gz | bgzip -c > chr${1}_ALL_multi.vcf.gz"

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
