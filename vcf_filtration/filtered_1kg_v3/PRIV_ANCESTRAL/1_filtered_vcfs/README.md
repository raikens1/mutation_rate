# filtered vcfs

These vcfs were filtered from the original 1,000 genomes files from voight_datasets in phaseIII_2013 on 7/25/2016. 

Filtration was done using vcftools(v0.1.12), and the following criteria were used:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed

In addition, ancestral allele (AA) and Allele frequency in the relevant population (POP_AF) was recoded into the INFO column of the vcf.

## File generation

These files were generated with vcftools using the --recode option, with the output piped to stdout using the --stdout option so that they could be compressed with bgzip.  For example, to filter European variants from the chromosome 22 vcf, one would run:

vcftools --gzvcf /project/voight_datasets/1kg/phaseIII_2013/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep ../../../populations/EUR_pops --bed ../../../nc_bedfiles/archive/nc_regions_header.bed  --mac 2 --remove-indels --remove-filtered-all --recode --recode-INFO AA --recode-INFO EUR_AF --stdout | bgzip -c > chr22_EUR_test.vcf.gz

(Here, "/project/.../ALL.chr22.phase3_shapeit2..." is a path to the chr22 vcf from 1kg phase III)

This runs vcftools with the following options:
        --gzvcf /project/voight_datasets/1kg/phaseIII_2013/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
        --keep ../../../populations/EUR_pops
        --mac 2
        --recode
        --remove-filtered-all
        --remove-indels
        --stdout
        --recode-INFO AA
        --recode-INFO EUR_AF
        --bed ../../../nc_bedfiles/archive/nc_regions_header.bed

The script filter_from_1kg.sh in this directory submits a job like this to the cluster for a given pop and chromosome.

All vcfs should be tabix indexed (tabix -p vcf <myfile.vcf.gz>)

## Notes

*Purpose:* These files will go on to be filtered for private variants and have multiallelic sites removed.

*Future Versions:* In the future, I would like to refilter these files with updated versions of the inclusion region bed files, and filtering out multiallelic variants at this step, rather than downstream.  Also, I'd like to remove the recoding of allele frequency in AFR only to avoid confusion. Since 1kg uses a different definition of AFR than we do, this number is not useful to us and will most-likely cause confusion.

## Bug history

Seven files in this directory were truncated during their initial filtration step. The original versions of these files are saved in the "truncated" subdirectory.  The current versions of these files (saved in this directory) were remade on 8/17/2017.  
