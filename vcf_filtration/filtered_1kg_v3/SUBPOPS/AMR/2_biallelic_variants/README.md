# private SNP lists

The files in this directory contain filtered variants found in each admixed American population.  They were filtered from the vcfs in ../1_multiallelic around 11/26/2016.  They are in vcf format.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - no multiallelic variants

## File generation

These files were generated using awk on the files in ../1_multiallelic to remove multiallelic variants.  The script awk_AMR.sh was originally used to do this work.  The commands listed are the same, however the file paths have changed, so that this shell script will no longer run.

## Notes:

*purpose* These SNPs will be filtered to exclude cosmopolitan variants so that they are AMR 'private.'
