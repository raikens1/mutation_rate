#!/bin/sh
# Merges and filters SGDP files:
#  - merge vcfs specified in VCF_LIST
#  - separate my chromosome
#  - extract regions in nc_bedfiles
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: SGDP_filter.sh VCF_LIST POP

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "# Merges and filters an SGDP file:
#  - merge vcfs specified in VCF_LIST
#  - separate my chromosome
#  - extract regions in nc_bedfiles
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: SGDP_filter.sh VCF_LIST POP"
    exit 1
fi

POP=${2}
FILENAME=$(basename -- "${1}")
HANDLE="${FILENAME%%.*}"

for i in {1..22}; do

    JOB_NAME="merge_chr${i}_${POP}"
    echo "${JOB_NAME}"
    cmd="module load bcftools; bcftools merge -l ${1} -m snps -R ../../nc_bedfiles/nc_chr${i}_regions_for_vcftools.bed -Oz -o chr_${i}_SGDP_${POP}.vcf.gz; tabix -p vcf chr_${i}_SGDP_${POP}.vcf.gz"
    echo "${cmd}"

    bsub -q voight_long -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"

done

JOB_NAME="merge_chrX_${POP}"
echo "${JOB_NAME}"
cmd="module load bcftools; bcftools merge -l ${1} -m snps -R ../../nc_bedfiles/nc_chrX_regions_for_vcftools.bed -Oz -o chr_X_SGDP_${POP}.vcf.gz; tabix -p vcf chr_X_SGDP_${POP}.vcf.gz"
echo "${cmd}"

bsub -q voight_normal -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"
