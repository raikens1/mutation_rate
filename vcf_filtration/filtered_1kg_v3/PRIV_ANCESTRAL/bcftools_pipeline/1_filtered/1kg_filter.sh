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
	FILENAME="/project/voight_datasets/1kg/phaseIII_2013/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
	HANDLE="chr_${i}_${POP}"
	
	JOB_NAME="filter_${HANDLE}"
	echo "${JOB_NAME}"
	cmd="module load bcftools; bcftools view ${FILENAME} -v snps -m 2 -M 2 -c 2 -f .,PASS -R ../../../../nc_bedfiles/nc_chr${i}_regions_for_vcftools.bed -S ../../../../populations/${POP}_pops --force-samples -Oz -o chr_${i}_${POP}_all_snps.vcf.gz; tabix -p vcf chr_${i}_${POP}_all_snps.vcf.gz"
	echo "${cmd}"

	bsub -q voight_long -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"
done

FILENAME="/project/voight_datasets/1kg/phaseIII_2013/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
HANDLE="chr_X_${POP}"  
JOB_NAME="filter_${HANDLE}"
echo "${JOB_NAME}"
cmd="module load bcftools; bcftools view ${FILENAME} -v snps -m 2 -M 2 -c 2 -f .,PASS -R ../../../../nc_bedfiles/nc_chrX_regions_for_vcftools.bed -S ../../../../populations/${POP}_pops --force-samples -Oz -o chr_X_${POP}_all_snps.vcf.gz; tabix -p vcf chr_X_${POP}_all_snps.vcf.gz"
echo "${cmd}"

bsub -q voight_long -o "${JOB_NAME}".out -e "${JOB_NAME}".err "${cmd}"
