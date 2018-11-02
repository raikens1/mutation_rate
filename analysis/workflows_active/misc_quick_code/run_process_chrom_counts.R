setwd("~/Documents/Research/mutation_rate/analysis/data")
source("../R_code/data_wrangling/process_chrom_counts.R")

ref <- read_delim("ref_files/3mer_mutations_ref.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

gw_counts <- read_delim("gw_counts/gw_3mer_counts.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

cc <- read.delim("~/Documents/Research/mutation_rate/vcf_filtration/SGDP/raw_file_prep/3_private/SGDP_SAS_chrom_counts.txt")

processed_cc <- process_chrom_counts(chrom_counts = cc,
                     gw_counts = gw_counts, 
                     subcontext_ref = ref)

write.table(processed_cc, "3mer/SGDP_SAS_counts.txt", quote = F, sep = "\t", row.names = F)