# private SNP lists

The files in this directory contain population private variants for each nonadmixed continental subpopulation in South Asia.  They were filtered from the vcfs in ../multiallelic around 7/31/2016.  They are not in vcf format.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - population 'private.' Any variant 'private' to an ancestral subpopulation must be present in that population and private to the continental group that it belongs to.
 - no multiallelic variants

## File generation

These files were generated using the awk and cut on the files in ../multiallelic to remove individual sample columns, comments, and multiallelic variants. The commands are saved in awk_SAS.sh

## Notes:

*purpose* these are used as the final SNP sets for counting

*future versions* in the future, I'd like to maintain the vcf formatting of these files instead of removing comments and sample columns.
