# data

Here, I've kept all of the main data files that I use for my analyses in R.

## Main File Types

There are a handful of main file types which are important here.  I save all my data as tab delimited text files so that they're easy to read from a simple text viewer.  The first two file types are meant to be read into R as dataframes, while rate matricies should be interpreted as labeled matricies.  Contexts or polymorphism types are always rows, with different information types listed across columns.

- *count dataframes* Usually saved as 'POP_kmer_counts.txt', for a given population (POP), and kmer. These are the capital data type used by my R scripts for most of my analyses. For each possible kmer polymorphism type. This file shows all relevant subcontexts, a genome-wide count of the times these contexts appears, inferred mutation rate, and counts of these polymorphisms on each chromosome (and genome wide total).  They are made by running cut_and_paste_counts.sh and then process_chrom_counts (in R) on the output of count_contexts.py for all 22 autosomes and X within POP.

- *gw_counts* Saved as gw_kmer_counts.txt for a given kmer. These are needed for certain analyses I do in R, and they are passed as an arguement to process_chrom_counts to generate my count dataframes. Each file contains the counts of the number of times each kmer sequence context appears on each chromosome within the included regions of hg19. This data is reformatted and used with permission from counts collected by Varun Aggarwala, former Voight lab member.

- *rate matricies* Saved as 'rates_kmer.txt' for a given kmer.  These are piped into heatmaps2 to make heatmaps of polymorphism types as in Figure 1A. For each kmer polymorphism type, these files show the inferred private mutation rate (per generation per site) across each of the 20 nonadmixed ancestral populations from 1kg (i.e. all populations from EUR, SAS, and EAS, and all from AFR except for ACB and ASW).

## Contents

### 3mer-5mer-7mer

These three directories contain the count dataframes private to each nonadmixed continental group (AFR, EUR, EAS, SAS), and cosmopolitan SNPs (COSMO).

### by_AF

These files contain counts of each 3mer in each continental group across bins of allele frequency constructed to have minimum size 4000.  As of now, they are not featured in any completed analyses for the main paper or supplement.  

### gw_counts

This is where I've kept all my gw_counts files for 3mer, 5mer, and 7mer.

### rate_profiles

Contains all my rate matricies for 3mer, 5mer, and 7mer.

### ref_files

Reference files of the subcontexts associated with each 3mer, 5mer, and 7mer.  These files are required as arguements to process_chrom_counts to construct count dataframes.

### subpops

Count datarames for each kmer for each 1kg population, including Americans.  A SNP is included in the set for a subpop if it is observed in that subpop and present in the 'private' set for its associated continental group.  In opposition to the 1kg designation, ACB and ASW are considered AMR subpopulations.
