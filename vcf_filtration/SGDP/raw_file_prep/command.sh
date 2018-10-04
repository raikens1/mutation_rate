#!/bin/sh
# Filters an SGDP file:
#  - extract biallelic sites
#  - remove singletons
#  - remove indels
#  - remove filtered sites
#  - remove regions in nc_bedfiles
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: command.sh SGDP_FILE

# TODO: fix naming of output file to not have double extension

vcftools --gzvcf ${1} --bed ../../nc_bedfiles/nc_regions_for_vcftools.bed --mac 2 --remove-indels --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout | bgzip -c > ${1}_filtered.vcf.gz
