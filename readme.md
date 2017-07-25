# mutation_rate

This repo contains all of the code and essential raw data
for replicating the central results of my report on differences
in polymorphism patterns across populations

## To Do

- [ ] Make final forms of filtered vcfs available on Ben's website
- [ ] Revise Rmarkdown code so that it runs neatly in this repo
- [ ] Organize, organize, organize (and don't break anything)
- [ ] Generate shell scripts (or examples of scripts) showing how I filter vcfs

## Pipeline

Here is a basic walkthrough of my workflow for filtering vcfs and arriving at
final results

### 1. Downloading raw data from 1kg

On 3/3/2016, I downloaded the most recently available release of the phase III 1kg variant call files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, checking md5 sums.

### 2. Filtering vcf files

Using vcftools v0.1.12b, I filtered these vcfs with the following criteria:
- quality == PASS
- variant type == SNP
- minor allele count >= 2
- within inclusion regions from nc_bedfiles/nc_regions
- (if relevant) observed in population of interest, based on individual codes in
	'vcf_filtration/populations/'

Using the vcftools vcf-isec command, I filtered out the following SNP sets:
- private continental (for most main-paper results)
- subpopulations (for heatmaps and analyses of signal 4)
- cosmopolitan (for model building)

The definitions of these SNP sets are defined in the Methods section of my manuscript.

After this step, I used awk to filter out any multiallelic variants.

Eventually, I will post to 'vcf_filtration' shell scripts (or examples of shell scripts) used to do this part of the data processing.

### 3. Counting vcf files

Using the python script helper_code/count_contexts.py, I compiled counts of each context type in each vcf under the 3mer, 5mer, and 7mer paradigms.

This script was written to parallel process by chromosome.  Once the script is run on each chromosomal file separately, these files could be compressed into a single output using the shell script cut_paste_counts.sh.  I often call the result a chrom_counts file.

### 4. Preprocessing in R

Once the counting has been done, I switch to R for the remainder of my analysis. The R function called 'preprocess_chrom_counts' (saved in 'analysis/R_code') takes in a chrom_counts file and generates a count dataframe, which is the master file type used for almost all of my downstream analyses. For each polymorphism/mutation type, it contains the number of contexts appearing in the human reference genome (from analysis/data/gw_context_counts), inferred genome-wide mutation rate, all subcontexts, and number of observed polymorphic sites in the included regions.  These files can be found in the 'R_data' directory.

The only other input type I've used is the rate matrix format, which is needed to construct heatmaps. These files, along with all my count dataframes, will be organized into a subdirectory of 'analysis/data' soon.

### 5. Analysis in R

Almost all of the analyses for this paper done in R are documented in several Rmarkdown analyses.  These can be found in the 'analysis/results' directory.

## Contents

### analysis

This folder contains markdown files, formatted data files, compiled pdfs and plots replicating all the analyses for my main paper and supplement.

### helper_code

The majority of my workhorse python and shell scripts can be found here. The subdirectory 'ref_files' contains text files for 3mer, 5mer, and 7mer, detailing the subcontexts that correspond to each possible polymorphism/mutation type.

### vcf_filtration

This directory contains reference files and shell scripts used to filter source vcf files from 1kg.
