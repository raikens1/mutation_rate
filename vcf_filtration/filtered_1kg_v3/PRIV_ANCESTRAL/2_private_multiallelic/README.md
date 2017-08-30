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

These files were generated using the vcf-isec command from vcf tools (v0.1.12) on the files in ../filtered_vcfs.  For example, to filter out European private variants on chromosome 22, one could run:

vcf-isec -c ../1_filtered_vcfs/chr22_EUR.vcf.gz ../1_filtered_vcfs/chr22_EAS.vcf.gz ../1_filtered_vcfs/chr22_SAS.vcf.gz ../1_filtered_vcfs/chr22_AFR.vcf.gz -f | bgzip -c > chr22_EUR_private_test.vcf.gz

The shell script 'filter_private.sh' is meant to submit such a job to the cluster for a fiven population and chromosome. 

All files should be tabix indexed (tabix -p vcf <myfile.vcf.gz>).

## Notes:

*purpose* these files will go on to be filtered for multiallelic variants and used as the final SNP sets for counting in ../private_SNP_lists

*future versions* in the future, I'd like to fold multiallelic filtering into the first filtration step from 1kg source files.

## Bug Report

On 8/16/2017, I noticed that several of these files were truncated in the filtration process.  These files were regenerated on 8/22/2017.  The original (truncated) copies are saved in the "truncated" subdirectory.
