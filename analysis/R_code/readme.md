# R Code

This directory contains all of the miscellaneous R scripts that have gone into my analysis over time.  There's some work to do here for organization and annotation.

## To Do

I need to check that all the code here still runs in this newly organized repo, annotate what is important in this readme, and move what is not important to the archive.

 - [ ] Everything in the main directory.

 Also, in data wrangling:

- [ ] Chisq_tables
- [X] import_and_process
- [ ] import_and_process_subpops
- [X] process_chrom_counts
- [ ] process_EAS_not_JPT (what is this anyway?  do I need it?)
- [ ] upload_AF_counts
- [ ] upload_and_process_AF_counts
- [ ] upload_continental
- [ ] upload_subpops
- [x] upload_subpops_3mer
- [x] upload_subpops_5mer
- [x] upload_subpops_7mer

## Contents

Need to populate this section with scripts as I update them.

### Data_wrangling subdirectory

- **import_and_process** uploads all the files necessary to run process_chrom_counts, sources that code, and runs it on a single chrom_counts file. Mostly quick code which I edit a lot for my own purposes.

- **process_chrom_counts** workhorse function for converting chrom_counts files to coutn dataframes.  Requires a gw_counts and a mutations_ref file to run.  If X == F, X will be excluded from the total count and rate calculations.

- **upload_subpops_kmer** uploads all the count dataframes for the 1kg subpops and puts tham in a list.  Meant to be run from the '/data' directory in ':/analaysis'.
