# Raw File Prep

The purpose of this directory is to hold SGDP files as they are processed into usable-sized vcfs. Not all intermediate files from this process can be saved in this directory due to memory constraints. 

This pipeline is a work in progress.  The steps are:

1. Download a group of files (from the same population, 'POP'):

   A. Using the SGDP data portal, select the desired vcfs and their tabix files.  Get the download links and copy them to a file.

   B. run the command:
        `sh wget_from_file.sh FILE_WITH_LINKS POP`

        Note: vcf downloads on sciget take about 10 minutes; tbi's are almost instantaneous.

2. Filter each file in parallel by chromosome 

   A. run the command:
        `sh filter_by_chrom.sh POP`
        Note: this script is not yet complete; filtration for larger chromosomes may take hours

3. Combine files of individuals into a single POP reference file

    B. This script is not yet complete


## Notes

Not certain yet whether it makes sense to remove singletons from this analysis.
