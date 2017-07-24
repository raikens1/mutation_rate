# mutation_rate

This repo contains all of the code and essential raw data
for replicating the central results of my report on differences
in polymorphism patterns across populations

## Pipeline

Here is a basic walkthrough of my workflow for filtering vcfs and arriving at final results

### 1. Downloading raw data from 1kg

On 3/3/2016, I downloaded the most recently available release of the phase III 1kg variant
call files from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, checking md5 sums

### 2. Filtering vcf files

Using vcftools v0.1.12b, I filtered the raw vcfs from 1kg based on the following criteria:
- quality == PASS
- variant type == SNP
- minor allele count >= 2
- within inclusion regions from nc_bedfiles/nc_regions
- (if relevant) observed in population of interest, based on individual codes in
	populations/

Using the vcftools vcf-isec command, I filtered out the following SNP sets:
- private continental (for most main-paper results)
- subpopulations (for heatmaps and analyses of signal 4)
- cosmopolitan (for model building)

The definitions of these SNP sets are defined in the Methods section of my manuscript.

After this step, I used awk to filter out any multiallelic variants.

### 3. Counting vcf files

Using the python script count_contexts.py, I compiled counts of each context type in each
vcf under the 3mer, 5mer, and 7mer paradigms.

This script was written to parallel process by chromosome.  Once the script is run on each
chromosomal file separately, these files could be compressed into a single output using
the shell script cut_paste_counts.sh.  I often call the result a chrom_counts file

### 4. Preprocessing in R

Once the counting has been done, switch to R for the remainder of my analysis.  The R function called 
'preprocess_chrom_counts' takes in a chrom_counts file and generates a count dataframe, which is the master
file type used for almost all of my downstream analyses. For each polymorphism/mutation type, it contains
the number of contexts appearing in the human reference genome (from gw_context counts), inferred 
genome-wide mutation rate, all subcontexts, and number of observed polymorphic sites in the included regions.

The only other input type I've used is the rate matrix format, which is needed to construct heatmaps.
These files, along with all my count dataframes, will be uploaded soon.

### 5. Analysis in R 

Almost all of the analyses for this paper done in R are documented in several Rmarkdown analyses.  
These will be added to the repo soon.

## Contents

### helper_code

The majority of my workhorse python and shell scripts can be found here.
The subdirectory 'ref_files' contains text files for 3mer, 5mer, and 7mer,
detailing the subcontexts that correspond to each possible polymorphism/mutation
type.

### nc_bedfiles

bedfiles containing noncoding regions from Aggarwalla and Voight, 2016.
These were used as the inclusion regions for the 1,000 genomes vcf files we analyzed.
nc_X_regions.bed contains just the noncoding regions on the X chromosome that we used.

### gw_context_counts

These files contain reformatted data from Varun Aggarwala, former Voight lab member.
Each file contains the counts of the number of times each sequence context appears on each 
chromosome within the included regions of hg19.

### populations

These are reference files of which 1kg individual codes come from which continental group and
subpopulation.

### singletons_excluded (not available in repo)

This directory is not tracked by github and not available on the github repo. It contains all 
of the vcf files filtered down from the raw 1,000 genomes variant call files which were used
this analysis.  Since these files are larger, we host the final versions of them on the voightlab
website.
