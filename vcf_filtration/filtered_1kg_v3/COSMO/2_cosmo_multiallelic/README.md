# multiallelic

The vcfs in this directory contain cosmopolitan variants for each chromosome.  They were filtered from the vcfs in PRIV_ANCESTRAL/mutliallelic.

These files have been filtered based on the following criteria:
 - minor allele count 2 or higher
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed
 - 'cosmopolitan.' A variant is cosmopolitan if it appears on the filtered SNP lists for two or more nonadmixed continental groups (AFR, EAS, SAS, EUR)

They are marked as 'multiallelic' because they still contain multiallelic variants, which still need to be filtered out.

## File generation

These files were generated using the vcf tools vcf-isec command (v0.1.12) on the files in PRIV_ANCESTRAL/filtered_vcfs.  The commands are saved in cosmo_multiallelic.sh.

## Notes:

*purpose* these files will go on to be filtered for multiallelic variants and used as the final SNP sets for counting in ../SNP_lists

*future versions* in the future, I'd like to fold multiallelic filtering into the first filtration step from 1kg source files.
