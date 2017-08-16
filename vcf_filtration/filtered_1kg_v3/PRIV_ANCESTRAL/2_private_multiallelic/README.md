# private_multiallelic

The vcfs in this directory contain population private variants for each nonadmixed continental group on each chromosome.  They were filtered from the vcfs in ../filtered_vcfs around 7/29/2016.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - population 'private.' Any variant 'private' to an ancestral continental group must be present in that group and in none of the other three populations.

They are marked as 'multiallelic' because they still contain multiallelic variants, which still need to be filtered out.

## File generation

These files were generated using the vcf-isec command from vcf tools (v0.1.12) on the files in ../filtered_vcfs

## Notes:

*purpose* these files will go on to be filtered for multiallelic variants and used as the final SNP sets for counting in ../private_SNP_lists

*future versions* in the future, I'd like to fold multiallelic filtering into the first filtration step from 1kg source files.