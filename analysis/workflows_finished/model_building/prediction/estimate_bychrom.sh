#!/bin/sh

# this script runs predict_polymorphisms on each autosome in AFR

# make array of True counts for each autosome
tcounts=( 519312 292615 332024 364931 424121 318476 375672 379454 323844 335170 325137 292178 303063 203793 143092 185353 144778 211668  98701 161045  91050  69774 )

for i in {1..22}; do
	handle="predict_COSMO_AFR_chr${i}"
	bsub -q voight_normal -o ${handle}.out -e ${handle}.err "predict_polymorphism.py ../../../../vcf_filtration/nc_bedfiles/nc_chr${i}_regions.bed ../results/AFR_3mer_model_params.txt 1 ${tcounts[`expr ${i} - 1`]}"
done
