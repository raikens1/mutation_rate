#!/bin/sh
# submits a job to the cluster to filter out multiallelic variants.
# USAGE: awk_multiallelic.sh POP CHROM

HAND="chr${2}_${1}_awk"

echo "Filtering ${1} variants from chromosome ${2}." 

cmd="bgzip -cd ../2_multiallelic/chr${2}_${1}_multi.vcf.gz | awk '{if (length(\$5) == 1 || \$1 ~ /^#/) print }' | bgzip -c > chr${2}_${1}.gz"

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"
