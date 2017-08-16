# filtered vcfs

These vcfs were filtered from the original 1,000 genomes files from voight_datasets in phaseIII_2013 on 7/25/2016. 

Filtration was done using vcftools(v0.1.12), and the following criteria were used:
 - minor allele count 2 or higher
 - only individuals from the appropriate continental group (based on identity files in :/vcf_filtration/populations)
 - within regions specified by nc_regions (in :/vcf_filtration/nc_regions/archive; not on github)
 - indels removed

In addition, ancestral allele (AA) and Allele frequency in the relevant population (POP_AF) was recoded into the INFO column of the vcf.

## File generation

Below is an example of the command to perform these filtration steps for chromosome 22 in Europeans:

bsub -q voight_normal -o ~/EUR_22_AA.out -i ~/EUR_22_AA.err "vcftools --gzvcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /project/voight_MR/raikens/SNP_AFs/EUR_pops --mac 2 --bed /project/voight_MR/raikens/SNP_AFs/nc_regions  --remove-indels --recode-INFO AA  --recode-INFO EUR_AF --out /project/voight_MR/raikens/SNP_AFs/chr22_EUR"

## Notes

*Purpose:* These files will go on to be filtered for private variants and have multiallelic sites removed.

*Future Versions:* In the future, I would like to refilter these files with updated versions of the inclusion region bed files, and filtering out multiallelic variants at this step, rather than downstream.  Also, I'd like to remove the recoding of allele frequency in AFR only to avoid confusion. Since 1kg uses a different definition of AFR than we do, this number is not useful to us and will most-likely cause confusion.
