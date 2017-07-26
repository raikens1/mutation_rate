library(readr)
source('C:/Users/VoightLab/Dropbox/SNP_rates/Code/data_wrangling/process_chrom_counts.R')

X3mer_mutations_ref <- read_delim("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/3mer_mutations_ref.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
gw_3mer_counts <- read_delim("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/gw_counts/gw_3mer_counts.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

setwd("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/singletons_excluded/3mer/")

COSMO_chrom_counts <- read_delim("chrom_counts/COSMO_chrom_counts.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
COSMO_3mer_counts <- process_chrom_counts(COSMO_chrom_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = COSMO_3mer_counts, "COSMO_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)