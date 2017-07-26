library(readr)
setwd("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/singletons_excluded/subpops/3mer/")
#setwd("C:/Users/Rocky/Dropbox/SNP_rates/Data/singletons_excluded/subpops/3mer")

#just some quick code to help me upload a bunch of subpop files

#EUROPE
TSI_EUR_3mer_counts <- read_delim("TSI_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

IBS_EUR_3mer_counts <- read_delim("IBS_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CEU_EUR_3mer_counts <- read_delim("CEU_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GBR_EUR_3mer_counts <- read_delim("GBR_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

FIN_EUR_3mer_counts <- read_delim("FIN_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#SOUTH ASIA
STU_SAS_3mer_counts <- read_delim("STU_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ITU_SAS_3mer_counts <- read_delim("ITU_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

BEB_SAS_3mer_counts <- read_delim("BEB_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GIH_SAS_3mer_counts <- read_delim("GIH_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PJL_SAS_3mer_counts <- read_delim("PJL_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#AFRICA
ESN_AFR_3mer_counts <- read_delim("ESN_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MSL_AFR_3mer_counts <- read_delim("MSL_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GWD_AFR_3mer_counts <- read_delim("GWD_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

LWK_AFR_3mer_counts <- read_delim("LWK_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

YRI_AFR_3mer_counts <- read_delim("YRI_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#EAST ASIA
CHB_EAS_3mer_counts <- read_delim("CHB_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

JPT_EAS_3mer_counts <- read_delim("JPT_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CHS_EAS_3mer_counts <- read_delim("CHS_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CDX_EAS_3mer_counts <- read_delim("CDX_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

KHV_EAS_3mer_counts   <- read_delim("KHV_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#AMERICA

ACB_AMR_3mer_counts   <- read_delim("ACB_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ASW_AMR_3mer_counts   <- read_delim("ASW_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CLM_AMR_3mer_counts   <- read_delim("CLM_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MXL_AMR_3mer_counts   <- read_delim("MXL_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PEL_AMR_3mer_counts   <- read_delim("PEL_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PUR_AMR_3mer_counts   <- read_delim("PUR_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

subpops.3mer.ancestral <- list(LWK_AFR_3mer_counts, ESN_AFR_3mer_counts, YRI_AFR_3mer_counts, MSL_AFR_3mer_counts, GWD_AFR_3mer_counts, 
                               TSI_EUR_3mer_counts, IBS_EUR_3mer_counts, GBR_EUR_3mer_counts, CEU_EUR_3mer_counts, FIN_EUR_3mer_counts,
                               PJL_SAS_3mer_counts, GIH_SAS_3mer_counts, ITU_SAS_3mer_counts, STU_SAS_3mer_counts, BEB_SAS_3mer_counts, 
                               CDX_EAS_3mer_counts, KHV_EAS_3mer_counts, CHS_EAS_3mer_counts, CHB_EAS_3mer_counts, JPT_EAS_3mer_counts) 


pops.ancestral <- c("LWK", "ESN", "YRI", "MSL", "GWD", 
                    "TSI", "IBS", "GBR", "CEU", "FIN", 
                    "PJL", "GIH", "ITU", "STU", "BEB",
                    "CDX", "KHV", "CHS", "CHB", "JPT")

subpops.3mer.all <- list(LWK_AFR_3mer_counts, ESN_AFR_3mer_counts, YRI_AFR_3mer_counts, MSL_AFR_3mer_counts, GWD_AFR_3mer_counts, 
                         ACB_AMR_3mer_counts, ASW_AMR_3mer_counts, CLM_AMR_3mer_counts, MXL_AMR_3mer_counts, PUR_AMR_3mer_counts, PEL_AMR_3mer_counts,
                         TSI_EUR_3mer_counts, IBS_EUR_3mer_counts, GBR_EUR_3mer_counts, CEU_EUR_3mer_counts, FIN_EUR_3mer_counts,
                         PJL_SAS_3mer_counts, GIH_SAS_3mer_counts, ITU_SAS_3mer_counts, STU_SAS_3mer_counts, BEB_SAS_3mer_counts, 
                         CDX_EAS_3mer_counts, KHV_EAS_3mer_counts, CHS_EAS_3mer_counts, CHB_EAS_3mer_counts, JPT_EAS_3mer_counts) 
                         

pops.all <- c("LWK", "ESN", "YRI", "MSL", "GWD", 
              "ACB", "ASW", "CLM", "MXL", "PUR", "PEL",
              "TSI", "IBS", "GBR", "CEU", "FIN", 
              "PJL", "GIH", "ITU", "STU", "BEB",
              "CDX", "KHV", "CHS", "CHB", "JPT")

# X3mer_subpop_rates <- data.frame(matrix(0, nrow = 1536, ncol = 20))
#  
# for (i in 1:20){
#   X3mer_subpop_rates[,i] <- subpops.3mer[[i]]$Rate
#  }
#  
# write.table(X3mer_subpop_rates, "X3mer_subpop_rates", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
