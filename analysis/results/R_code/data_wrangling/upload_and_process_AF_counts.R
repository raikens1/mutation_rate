#This script is for uploading the output files from counter > postprocessing
#and reformatting them into a version that's easy to read and use in R

library(readr)
setwd("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/singletons_excluded/by_AF/v1_m_4000")

X3mer_mutations_ref <- read_delim("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/3mer_mutations_ref.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

process_AF_counts<- function(raw.dat, bins.log, mut.ref){
  colnames(raw.dat) <- bins.log$`#BIN`
  raw.dat <- cbind(mut.ref, raw.dat)
  #raw.dat <- rbind(raw.dat, c(NA, NA, bins.log$AVG_AF), c(NA, NA, bins.log$j)) #messes with the formatting.  probably best to save this info separately
  return(raw.dat)
}

#Africa

AFR_AF_counts <- read_delim("AFR_AF_counts_m_40000.txt", " ", escape_double = FALSE, trim_ws = TRUE)
AFR_AF_bins <- read_delim("AFR_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
AFR_AF_counts <- process_AF_counts(AFR_AF_counts, AFR_AF_bins, X3mer_mutations_ref)
write.table(AFR_AF_counts, "AFR_AF_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Europe

EUR_AF_counts <- read_delim("EUR_AF_counts_m_40000.txt", " ", escape_double = FALSE, trim_ws = TRUE)
EUR_AF_bins <- read_delim("EUR_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
EUR_AF_counts <- process_AF_counts(EUR_AF_counts, EUR_AF_bins, X3mer_mutations_ref)
write.table(EUR_AF_counts, "EUR_AF_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#East Asia

EAS_AF_counts <- read_delim("EAS_AF_counts_m_40000.txt", " ", escape_double = FALSE, trim_ws = TRUE)
EAS_AF_bins <- read_delim("EAS_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
EAS_AF_counts <- process_AF_counts(EAS_AF_counts, EAS_AF_bins, X3mer_mutations_ref)
write.table(EAS_AF_counts, "EAS_AF_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#South Asia

SAS_AF_counts <- read_delim("SAS_AF_counts_m_40000.txt", " ", escape_double = FALSE, trim_ws = TRUE)
SAS_AF_bins <- read_delim("SAS_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
SAS_AF_counts <- process_AF_counts(SAS_AF_counts, SAS_AF_bins, X3mer_mutations_ref)
write.table(SAS_AF_counts, "SAS_AF_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
