#!/bin/shc
# Filters an SGDP file:
#  - extract biallelic sites
#  - remove indels
#  - remove filtered sites
#  - remove regions in nc_bedfiles
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: SGDP_filter.sh SGDP_FILE

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "# Filters an SGDP file:
#  - extract biallelic sites
#  - remove indels
#  - remove filtered sites
#  - remove regions in nc_bedfiles
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: SGDP_filter.sh SGDP_FILE"
    exit 1
fi

FILENAME=$(basename -- "${1}")
HANDLE="${FILENAME%%.*}"

for i in {1..1}; do

    JOB_NAME="filter_chr${i}_${HANDLE}"
    echo "${JOB_NAME}"
    cmd="vcftools --gzvcf ${1} --bed ../../nc_bedfiles/nc_regions_for_vcftools.bed --chr ${i} --remove-indels --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout | bgzip -c > chr_${i}_${HANDLE}.vcf.gz"
    echo "${cmd}"

    bsub -q voight_normal -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"

done
