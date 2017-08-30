#!/bin/sh
# submits a job to the cluster to filter out multiallelic variants.
# USAGE: awk_multiallelic.sh CHROM

HAND="chr${1}_ALL_awk"

echo "Filtering ALL-type variants from chromosome ${1}." 

cmd="bgzip -cd ../2_multiallelic/chr${1}_ALL_multi.vcf.gz | awk '{if (length(\$5) == 1 || \$1 ~ /^#/) print }' | bgzip -c > chr${1}_ALL.gz"

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
