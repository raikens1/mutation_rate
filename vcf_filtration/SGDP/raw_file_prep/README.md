# Raw File Prep

The purpose of this directory is to hold SGDP files as they are processed into usable-sized vcfs. Not all intermediate files from this process can be saved in this directory due to memory constraints. For details on the precise files that were downloaded for this pipeline, see DATA_NOTES.md

This pipeline is a work in progress.  The steps are:

0. Download a group of files (from the same population, 'POP'):

   A. Using the SGDP data portal, select the desired vcfs and their tabix files.  Get the download links and copy them to a file.

   B. run the command:
        `sh wget_from_file.sh FILE_WITH_LINKS POP`

        Note: vcf downloads on sciget take about 10 minutes; tbi's are almost instantaneous.

1. Marge files of each individual into a population reference file (in parallel by chromosome) 

   A. run the command:
        `sh SGDP_merge.sh FILE_LIST POP`

        Notes: 
		- This script will queue one job on scisub for each chromosome.
		- A list of file names to merge is required.  To make this, one might run something like `for file in SGDP*.vcf.gz; do echo ${file} >> FILE_LIST; done
		 - Extracts only the noncoding regions of each chromosome.

2. Filter the chromosome reference files to include only biallelic SNPs with <= 20% missing data and no singletons

    A. run the command:
	`sh SGDP_filter FILENAME`

	Notes: 
		-This should be run in a loop over each chromosome.  To do this, run `for file in chr*_.vcf.gz; do sh SGDP_filer ${file}; done
	
3. Use bcftools isec to retrieve population-specific variants

    A. Go to the subdirectory `2_filtered`

    B. Run the command `sh ../SGDP_bcftools_isec.sh POP CHROM` for all 23 chromosomes and 4 pops.  Best to use a `for` loop over the numeric chromosomes.


