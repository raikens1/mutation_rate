# mutation_rate

This repo contains all of the code and essential raw data
for replicating the central results of my report on differences
in polymorphism patterns across populations

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
