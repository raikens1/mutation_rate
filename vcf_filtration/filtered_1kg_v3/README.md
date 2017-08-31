#filtered_1kg_v3

This directory contains filtered vcfs from the 1kg source files available from the 1,000 genomes website.  The files in these directories are explained more thoroughly in their own READMEs

This version is distinct from the previous two in that it excludes singletons and uses vcf tools to do most of the filtration work.

## Bug

There is an issue with some of the files in this runthrough, discovered 8/16/2017.  The PRIV_ANCESTRAL vcfs for certain chromosomes in certain ancestral populations appear to be truncated somehow due to cluster computing errors.  Since the lists of private SNPs generated from those files is important for defining other subsets of SNPs, many of the other files (e.g. SUBPOP vcfs) are likely to be incorrect as well.  

At the 1_filtered_vcfs step, the following files are suspect:
AFR: chr 4
EUR: chr 3, 4
EAS: chr 3, 4
SAS: chr 3, 4

At the 2_private_multiallelic step, the following files are suspect:
AFR: chr 2, 3, 4, 6
EUR: chr 1, 2, 3, 4, 5, 6
EAS: chr 1, 2, 3, 4, 5, 6
SAS: chr 2, 3, 4, 5, 6

These files are suspect because running the command:
bgzip -cd <file> | wc 
on these files generates the error:
Could not read XXXX bytes.  Error 4

We first noticed this bug because we noted that the numbers of SNPs on chromosomes 2, 3, 4, and 6 were far lower than expected in our probablistic mutational model.


