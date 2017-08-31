#!/bin/sh
# submits a job to filter a 1kg vcf file for a given population and chromosome
# USAGE: filter_from_1kg.sh POP CHROM

HAND="chr${2}_${1}"

if [ ${2} == "X" ]; then
	echo "Filtering ${1} variants from X chromosome." 
	cmd="vcftools --gzvcf /project/voight_datasets/1kg/phaseIII_2013/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz --keep ../../../populations/${1}_pops --bed ../../../nc_bedfiles/nc_regions_for_vcftools.bed  --mac 2 --remove-indels --remove-filtered-all --recode --recode-INFO AA --recode-INFO ${1}_AF --stdout | bgzip -c > chr${2}_${1}.vcf.gz"
else
	echo "Filtering ${1} variants from chr${2}"
	cmd="vcftools --gzvcf /project/voight_datasets/1kg/phaseIII_2013/ALL.chr${2}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep ../../../populations/${1}_pops --bed ../../../nc_bedfiles/nc_regions_for_vcftools.bed  --mac 2 --remove-indels --remove-filtered-all --recode --recode-INFO AA --recode-INFO ${1}_AF --stdout | bgzip -c > chr${2}_${1}.vcf.gz"
fi

bsub -q voight_normal -o "${HAND}".out -e "${HAND}".err "${cmd}"