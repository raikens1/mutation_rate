# multiallelic

The vcfs in this directory contain population private variants for each nonadmixed subpopulation in Africa.  They were filtered from the vcfs in PRIV_ANCESTRAL/2_mutliallelic around 7/29/2016.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - only individuals from the appropriate subpopulation (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - population 'private.' Any variant 'private' to an ancestral subpopulation must be present in that population and private to the continental group that it belongs to.

They are marked as 'multiallelic' because they still contain multiallelic variants, which still need to be filtered out.

## File generation

These files were generated using the vcf tools (v0.1.12) on the files in PRIV_ANCESTRAL/filtered_vcfs, extracting private polymorphisms present in each subpopulation.  The script AFR_SUBPOP_multiallelic.sh will submit a job to the cluster to generate this file for a given population and chromosome.

## Notes:

*purpose* these files will go on to be filtered for multiallelic variants and used as the final SNP sets for counting in ../SNP_lists

*future versions* in the future, I'd like to fold multiallelic filtering into the first filtration step from 1kg source files.
