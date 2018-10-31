#!/bin/sh
# Merges and filters SGDP files:
#  - extract biallelic sites
#  - remove indels
#  - remove sites with >20% missing data
#  - remove singletons
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: SGDP_filter.sh FILENAME

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "# Merges and filters an SGDP file:
#  - merge vcfs	specified in VCF_LIST
#  - separate my chromosome
#  - extract biallelic sites
#  - remove indels
#  - remove filtered sites
#  - remove regions in nc_bedfiles
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: SGDP_filter.sh FILENAME"
    exit 1
fi

FILENAME=$(basename -- "${1}")
HANDLE="${FILENAME%%.*}"

JOB_NAME="filter_${HANDLE}"
echo "${JOB_NAME}"
cmd="module load bcftools; bcftools view ${1} -v snps -m 2 -M 2 -c 2 | vcftools --vcf - --max-missing 0.2 --recode --stdout | bgzip -c > ${HANDLE}_all_snps.vcf.gz; tabix -p vcf ${HANDLE}_all_snps.vcf.gz"
echo "${cmd}"

bsub -q voight_long -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"

