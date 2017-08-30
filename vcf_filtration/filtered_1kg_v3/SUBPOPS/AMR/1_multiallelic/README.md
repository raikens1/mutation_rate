# biallelic variants

The files in this directory contain filtered variants found in each admixed American population.  They were filtered from the vcfs in from the 2013 1kg phase III release around 8/4/2016.  They are in vcf format.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed

## File generation

These files were generated using vcftools on the original 1kg vcfs to filter variants based on the criteria above.  The script AMR_multiallelic.sh was originally used to do this.  The commands are the same, but the paths have changed, so this script will no longer run and needs to be rewritten.

## Notes:

*purpose* These SNPs must be filtered to exclude multiallelic variants and Cosmopolitan variants.
