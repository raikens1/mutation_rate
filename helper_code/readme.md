# helper_code

This directory contains the code I use to parse vcf files for information about allele frequency and polymorphism counts.

- [ ] add to this readme

## Code

There are a couple main scripts here:

 - **bin_by_freq.py** Given a vcf file which is sorted by allele frequency, write out vcf files binned by allele frequency with a prespecified minimum number of SNPs
 - **count_contexts.py** This is my workhorse script that takes in a vcf file and sequence context window size and prints a text file of the counts of each possible polymorphism type observed in the vcf.  
 - **cut_paste_counts.sh** Since count_contexts.py is made to parallel process by chromosome, this is just a quick shell script to merge the output files from count_contexts.py across all autosomes + X into a single output file.
 - **filter_by_context** given a text file of polymorphism types and an input vcf, this script will generate a .bed file of all the SNPs of those polymorphism types.
 - **class_counter** This defines the counter class relied upon by count_contexts.py
 - **get_context.py** Defines a simple function that returns the sequence context of a variant based on hg19.  I didn't archive this because I thought it might come in handy someday.

## ref_files

This directory contains reference files that I use for some data processing steps.  Each ref file contains an ordered list of polymorphism types under each sequence context paradigm, along with all possible subcontexts.
