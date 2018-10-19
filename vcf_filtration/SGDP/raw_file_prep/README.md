# Raw File Prep

The purpose of this directory is to hold SGDP files as they are processed into usable-sized vcfs. Not all intermediate files from this process can be saved in this directory due to memory constraints. 

This pipeline is a work in progress.  The steps are:

0. Download a group of files (from the same population, 'POP'):

   A. Using the SGDP data portal, select the desired vcfs and their tabix files.  Get the download links and copy them to a file.

   B. run the command:
        `sh wget_from_file.sh FILE_WITH_LINKS POP`

        Note: vcf downloads on sciget take about 10 minutes; tbi's are almost instantaneous.

1. Filter each file in parallel by chromosome 

   A. run the command:
        `sh SGDP_filter.sh FILENAME`

        Note: this script will queue one job on scisub for each chromosome.  Most jobs seem to take no more than 10 minutes.  Needs to be run in a loop over many files.  For example: for file in SGDP_AFR_*: do sh SGDP_filter.sh ${file}; done

2. Combine files of individuals into a single POP reference file

    A. run the command:
	`sh merge_by_chrom VCFLIST`.  VCFLIST should be a text file of the names of the vcfs to be merged, with the 'chr_N_' prefix removed.  A good way to do this might be: `for file in chr_1_SGDP_*vcf.gz; do echo ${file#chr_1_} >> VCFLIST; done`
	
3. (Still to do)  Remove sites with 20% or more missing data

    A. (need to write this script) probably a python script which eliminates undesired rows while keeping headers.  May subsequently need to re-tabix

3b? FILTER SINGLETONS?  Not sure whether this makes sense given the size of the dataset

4. (script mostly written) Use vcf-isec to retrieve population-specific variants


## Notes

Not certain yet whether it makes sense to remove singletons from this analysis.
