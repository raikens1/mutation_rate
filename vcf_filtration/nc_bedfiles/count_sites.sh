#!/bin/sh
# runs site_counter for all chromosomes
# USAGE: count_sites.sh FLANK

for i in {1..22}; do
	BED="nc_chr${i}_regions.bed"
	HAND="chr${i}_gw_context_counts"
	bsub -q voight_normal -o "$HAND".out -e "$HAND".err site_counter.py "$BED" $1;

done

BED="nc_chrX_regions.bed"
HAND="chrX_gw_context_counts"
bsub -q voight_normal -o "$HAND".out -e "$HAND".err site_counter.py "$BED" $1;
