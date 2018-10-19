#!/bin/shc
# Merge a list of SGD files, separated by chromosome:
# Takes a list of file endings - this should be the names of the files to be merged with any annotations after the first dot removed.
# E.g. SGDP_AFR_LP60044XXXX.vcf.gz (NOT "chr_1_SGDP_AFR...")
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: merge_by_chrom.sh VCFLIST POP

if [ $# -ne 2 ]
  then
    echo "Incorrect number of arguments supplied."
    echo "# Merge a list of SGD files, separated by chromosome:
# Takes	a list of file endings - this should be the names of the files to be merged with the "chr_N_" prefix removed. 
# E.g. SGDP_AFR_LP60044XXXX.vcf.gz (NOT	"chr_1_SGDP_AFR...")
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: merge_by_chrom.sh VCFLIST POP"
    exit 1
fi

FILENAME=$(basename -- "${1}")
HANDLE="${FILENAME%%.*}"

for i in {1..22}; do
    # make file of all things to merge
    sed -e "s/^/chr_${i}_/" ${FILENAME} > chr_${i}_${FILENAME}

    JOB_NAME="merge_chr${i}_${HANDLE}"
    echo "${JOB_NAME}"

    cmd="module load bcftools; bcftools merge -l chr_${i}_${FILENAME} -Oz -o chr_${i}_SGDP_${2}.vcf.gz; tabix -p vcf chr_${i}_SGDP_${2}.vcf.gz"
    echo "${cmd}"

    bsub -q voight_normal -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"

done

# make file	of all things to merge
sed -e 's/^/chr_X_/' ${FILENAME} > chr_X_${FILENAME}

JOB_NAME="merge_chrX_${HANDLE}" 
echo "${JOB_NAME}"

cmd="module load bcftools; bcftools merge -l chr_X_${FILENAME} -Oz -o chr_X_SGDP_${2}.vcf.gz; tabix -p vcf chr_X_SGDP_${2}.vcf.gz"                                                         
echo "${cmd}" 

bsub -q voight_normal -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"
