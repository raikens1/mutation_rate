# private SNP lists

The files in this directory contain population private variants for each admixed American population.  They were filtered from the vcfs in ../2_biallelic_variants around 5/12/2017.  Most are not in vcf format (the sample columns have been cut out).

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - population 'private.' Any variant 'private' to an admixed american subpopulation if it is observed in that subpopulation and not counted as cosmopolitan.
 - no multiallelic variants

## File generation

These files were generated using the vcf-isec tool on the files in ../2_biallelic_variants to remove cosmopolitan variants so that these sets are AMR 'private.'  The script filter_private.sh will submit a job to the voightlab cluster to generate these files for a given population and chromosome. (The original script used to do this work is now deprecated. It is saved in ./archive/)

## Notes:

*purpose* these are used as the final SNP sets for counting.

*future versions* in the future, I'd like to maintain the vcf formatting of these files instead of removing comments and sample columns.
