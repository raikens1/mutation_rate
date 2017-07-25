# helper_code

This directory contains the code I use to parse vcf files for information about allele frequency and polymorphism counts.

## To Do

- [ ] Update count_contexts.py to a new version which uses biopython rather than Varun's code
- [ ] Move bin_by_freq and Varun's code to archive.

## Code

There are a couple main scripts here:

 - **bin_by_freq.py** Given a vcf file which is sorted by allele frequency, write out vcf files binned by allele frequency with a prespecified minimum number of SNPs
 - **count_contexts.py** This is my workhorse script that takes in a vcf file and sequence context window size and prints a text file of the counts of each possible polymorphism type observed in the vcf.  
 - **cut_paste_counts.sh** Since count_contexts.py is made to parallel process by chromosome, this is just a quick shell script to merge the output files from count_contexts.py across all autosomes + X into a single output file.
 - **filter_by_context** given a text file of polymorphism types and an input vcf, this script will generate a .bed file of all the SNPs of those polymorphism types.

 The remaining scripts in this directory were written by Dr. Varun Aggarwala (a former Voight lab member).  They support the helper function get_seq_context_variant, which is called in some of the scripts above.

 - find_context.py
 - function_wrapper.py
 - modules.py

## ref_files

This directory contains reference files that I use for some data processing steps.  Each ref file contains an ordered list of polymorphism types under each sequence context paradigm, along with all possible subcontexts.
