# vcf_filtration

This directory contains all the reference files needed to replicate the steps I take to filter the source vcfs from 1kg.

## To Do

- [ ] rename singletons_excluded to filtered_vcfs_v3
- [ ] add a folder of shell scripts/examples

## Contents

### nc_bedfiles

bedfiles containing noncoding regions from Aggarwalla and Voight, 2016.
These were used as the inclusion regions for the 1,000 genomes vcf files we
analyzed.
nc_X_regions.bed contains just the noncoding regions on the X chromosome that
 we used.

### populations

These are reference files of which 1kg individual codes come from which
continental group and subpopulation.

### filtered_vcfs_v3 (not available in repo)

This directory is not tracked by github and not available on the github repo.
It contains all of the vcf files filtered down from the raw 1,000 genomes
variant call files which were used this analysis.  Since these files are larger,
we'll host the final versions of them on the voightlab website.
