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

    JOB_NAME="bcftools_filter_chr${i}_${HANDLE}"
    echo "${JOB_NAME}"
    cmd="module load bcftools; bcftools view -r ${i} -V indels -m 2 -M 2 -Oz -o  chr_${i}_${HANDLE}_bcftools.vcf.gz ${1}"
    echo "${cmd}"

    bsub -q voight_normal -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"

done
