#!/bin/sh
# submits a job to filter out AMR 'private' variants - those that are seen in AMR and not Cosmopolitan
# USAGE: filter_private_v3.sh POP CHROM

HAND="chr${2}_${1}"

echo "Filtering variants private to ${1} on chromosome ${1}."

cmd="vcf-isec -c ../2_biallelic_variants/chr${2}_${1}.vcf.gz ../../../COSMO/2_cosmo_multiallelic/chr${2}_COSMO_private_multi.vcf.gz -f | bgzip -c > chr${2}_${1}_private.gz"

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
