#This script is for uploading the output files from counter > postprocessing
#and reformatting them into a version that's easy to read and use in R

library(readr)
setwd("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/singletons_excluded/by_AF/v1_m_4000")

#Africa

AFR_AF_bins <- read_delim("AFR_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
AFR_AF_counts <- read_delim("AFR_AF_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#Europe

EUR_AF_bins <- read_delim("EUR_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
EUR_AF_counts <- read_delim("EUR_AF_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#East Asia

EAS_AF_bins <- read_delim("EAS_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
EAS_AF_counts <- read_delim("EAS_AF_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#South Asia

SAS_AF_bins <- read_delim("SAS_sortedSNPs_bins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
SAS_AF_counts <- read_delim("SAS_AF_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
