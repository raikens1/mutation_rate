# private SNP lists

The files in this directory contain variants for each chromosome which are found in all four nonadmixed continental groups.  They were filtered from the vcfs in ../2_multiallelic on 8/29/2017.  They are in vcf format.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - found (based on all these filtration criteria) in all four nonadmixed continental groups
 - no multiallelic variants

## File generation

These files were generated using the awk and cut on the files in ../2_multiallelic to remove multiallelic variants. The commands are saved in awk_multiallelic.sh

## Notes:

*purpose* these are used as the final SNP sets for counting.
