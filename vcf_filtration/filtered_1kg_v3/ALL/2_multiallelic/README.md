# multiallelic

The vcfs in this directory contain variants on each chromosome which are found in all four nonadmixed ancestral groups.  They were filtered from the vcfs in PRIV_ANCESTRAL/1_filtered_SNP_lists.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - Found in 'ALL' nonadmixed continental groups. A variant is in this set if it appears on the filtered SNP lists for all four nonadmixed continental groups (AFR, EAS, SAS, EUR)

They are marked as 'multiallelic' because they still contain multiallelic variants, which still need to be filtered out.

## File generation

These files were generated using the vcf tools vcf-isec command (v0.1.12) on the files in PRIV_ANCESTRAL/1_filtered_vcfs.  The commands are saved in ALL_multiallelic.sh.

## Notes:

*purpose* these files will go on to be filtered for multiallelic variants and used as the final SNP sets for counting in ../3_SNP_lists
