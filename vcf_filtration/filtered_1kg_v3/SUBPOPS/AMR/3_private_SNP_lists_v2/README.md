private_SNP_lists_v2
#########################

This directory will contain all American variants from monoallelic_variants which are 'private' to AMR based on the following definition:

A variant is 'private' to AMR if it is not one of the cosmopolitan variants (COSMO).
Cosmopolitan variants are those which are found in two or more of: AFR, EUR, EAS, or SAS

The .gz files in this directory were made by running filter_private_v2.sh in this working directory.
The SNPs from these files were counted using shell scripts count_all_3mers.sh, count_all_5mers.sh, count_all_7mers.sh. 
The outputs from these counting exploits can be found in their respective subdirectories.


