This directory contains files of the 1000 genomes phase III genomes 
from each of the four superpopulations: AFR, EAS, EUR, and SAS

These files are used with vcftools "--keep" flage  in the first step
 which fiters out singletons, noncoding SNPs, and indels from the 
original 1000 genomes .vcf files 

Note: Only the YRI, LWK, GWD, MSL and ESN populations from Africa
were considered in this analysis.  The ASW and ACB groups were
excluded because these populations are likely to be admixed.

Example command to extract European populations file:
awk '{if ($2=="IBS" || $2 == "GBR" || $2 == "FIN" || $2 == "TSI" || $2 == "CEU") print $1;}' pop_codes.txt > EUR_pops

The pop_codes file contains the first two columns of the 1000 genomes
file mapping individuals to populations.

