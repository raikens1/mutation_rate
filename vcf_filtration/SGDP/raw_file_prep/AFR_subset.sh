#!/bin/sh
# Filters a 1kg file from the reference:
#  - extract biallelic sites
#  - remove indels
#  - remove filtered sites
#  - extract regions in nc_bedfiles
#  - extract individuals in POP
#  - remove singletons
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: 1kg_filter.sh POP

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "# Filters a 1kg file from the reference:
#  - extract biallelic sites
#  - remove indels
#  - remove filtered sites
#  - extract regions in nc_bedfiles
#  - extract individuals in POP
#  - remove singletons
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: 1kg_filter.sh POP"
    exit 1
fi

POP=${1}

for i in {1..22}; do
	FILENAME="chr_${i}_SGDP_AFR_private.vcf.gz"
	HANDLE="chr_${i}_AFR_subset"
	
	JOB_NAME="${HANDLE}"
	echo "${JOB_NAME}"
	cmd="module load bcftools; bcftools view ${FILENAME} -S ../../AFR_1kg_overlap.txt -Oz -o chr_${i}_AFR_subset.vcf.gz; tabix -p vcf chr_${i}_AFR_subset.vcf.gz"
	echo "${cmd}"

	bsub -q voight_long -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"
done

FILENAME="chr_X_SGDP_AFR_private.vcf.gz"
HANDLE="chr_X_AFR_subset"  
JOB_NAME="${HANDLE}"
echo "${JOB_NAME}"
cmd="module load bcftools; bcftools view ${FILENAME} -S ../../AFR_1kg_overlap.txt -Oz -o chr_X_AFR_subset.vcf.gz; tabix -p vcf chr_X_AFR_subset.vcf.gz"
echo "${cmd}"

bsub -q voight_long -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"
