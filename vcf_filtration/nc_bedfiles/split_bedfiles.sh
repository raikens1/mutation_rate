#!/bin/sh

# this script is for splitting Varun's nc_regions.bed file by chromosome.

for i in {1..22}; do

	fgrep -w chr${i} nc_regions_for_vcftools.bed > nc_chr${i}_regions_for_vcftools.bed;

done
