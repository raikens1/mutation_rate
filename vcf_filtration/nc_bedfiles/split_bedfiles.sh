#!/bin/sh

# this script is for splitting Varun's nc_regions.bed file by chromosome.

for i in {1..22}; do

	awk -v chrom="${i}" '{if ($1 == chrom) print $0;}' nc_regions_for_vcftools.bed > nc_chr${i}_regions_for_vcftools.bed;

done

awk '{if ($1 == "X") print $0;}' nc_regions_for_vcftools.bed > nc_chrX_regions_for_vcftools.bed;
