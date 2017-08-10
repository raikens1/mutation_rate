#!/bin/sh
# quick script for cutting and pasting together the outputs of site_counter.py

paste nc_chr1_regions_site_counts.txt nc_chr2_regions_site_counts.txt nc_chr3_regions_site_counts.txt nc_chr4_regions_site_counts.txt nc_chr5_regions_site_counts.txt nc_chr6_regions_site_counts.txt nc_chr7_regions_site_counts.txt nc_chr8_regions_site_counts.txt nc_chr9_regions_site_counts.txt  nc_chr10_regions_site_counts.txt  nc_chr11_regions_site_counts.txt  nc_chr12_regions_site_counts.txt  nc_chr13_regions_site_counts.txt  nc_chr14_regions_site_counts.txt  nc_chr15_regions_site_counts.txt  nc_chr16_regions_site_counts.txt  nc_chr17_regions_site_counts.txt  nc_chr18_regions_site_counts.txt  nc_chr19_regions_site_counts.txt  nc_chr20_regions_site_counts.txt  nc_chr21_regions_site_counts.txt  nc_chr22_regions_site_counts.txt nc_chrX_regions_site_counts.txt > temp
cut -f1-2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54 temp > gw_"${1}"mer_chrom_counts.txt
rm temp