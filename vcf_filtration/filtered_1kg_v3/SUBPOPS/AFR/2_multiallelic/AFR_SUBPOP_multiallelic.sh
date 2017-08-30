#!/bin/sh
# submits a job to filter variants from a subpopulation out of a Priv multiallelic vcf file.
# USAGE: SUBPOP_multiallelic.sh POP CHROM
# I had to make a special version of this script for AFR because I missed an underscore when I named these files originally and I need to keep things consistent.

HAND="chr${2}_${1}"
SUPERPOP=${1: -3}

echo "Filtering ${1} variants from chr${2}"
cmd="vcftools --gzvcf ../../../PRIV_ANCESTRAL/2_private_multiallelic/chr${2}_${SUPERPOP}_private_multi.vcf.gz --keep ../../../../populations/${1} --mac 1 --recode --recode-INFO AA --stdout |  bgzip -c > chr${2}_${1}multi.vcf.gz"

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"

