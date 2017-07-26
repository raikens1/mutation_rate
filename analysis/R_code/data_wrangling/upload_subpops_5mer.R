library(readr)

# just some quick code to help me upload a bunch of subpop files
# meant to be run from the 'data' directory

# EUROPE
TSI_EUR_5mer_counts <- read_delim("subpops/5mer/TSI_EUR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

IBS_EUR_5mer_counts <- read_delim("subpops/5mer/IBS_EUR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CEU_EUR_5mer_counts <- read_delim("subpops/5mer/CEU_EUR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GBR_EUR_5mer_counts <- read_delim("subpops/5mer/GBR_EUR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

FIN_EUR_5mer_counts <- read_delim("subpops/5mer/FIN_EUR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# SOUTH ASIA
STU_SAS_5mer_counts <- read_delim("subpops/5mer/STU_SAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ITU_SAS_5mer_counts <- read_delim("subpops/5mer/ITU_SAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

BEB_SAS_5mer_counts <- read_delim("subpops/5mer/BEB_SAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GIH_SAS_5mer_counts <- read_delim("subpops/5mer/GIH_SAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PJL_SAS_5mer_counts <- read_delim("subpops/5mer/PJL_SAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# AFRICA
ESN_AFR_5mer_counts <- read_delim("subpops/5mer/ESN_AFR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MSL_AFR_5mer_counts <- read_delim("subpops/5mer/MSL_AFR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GWD_AFR_5mer_counts <- read_delim("subpops/5mer/GWD_AFR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

LWK_AFR_5mer_counts <- read_delim("subpops/5mer/LWK_AFR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

YRI_AFR_5mer_counts <- read_delim("subpops/5mer/YRI_AFR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# EAST ASIA
CHB_EAS_5mer_counts <- read_delim("subpops/5mer/CHB_EAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

JPT_EAS_5mer_counts <- read_delim("subpops/5mer/JPT_EAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CHS_EAS_5mer_counts <- read_delim("subpops/5mer/CHS_EAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CDX_EAS_5mer_counts <- read_delim("subpops/5mer/CDX_EAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

KHV_EAS_5mer_counts <- read_delim("subpops/5mer/KHV_EAS_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# AMERICA

ACB_AMR_5mer_counts <- read_delim("subpops/5mer/ACB_AMR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ASW_AMR_5mer_counts <- read_delim("subpops/5mer/ASW_AMR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CLM_AMR_5mer_counts <- read_delim("subpops/5mer/CLM_AMR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MXL_AMR_5mer_counts <- read_delim("subpops/5mer/MXL_AMR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PEL_AMR_5mer_counts <- read_delim("subpops/5mer/PEL_AMR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PUR_AMR_5mer_counts <- read_delim("subpops/5mer/PUR_AMR_5mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

subpops.5mer.ancestral <- list(LWK_AFR_5mer_counts, ESN_AFR_5mer_counts, YRI_AFR_5mer_counts, MSL_AFR_5mer_counts, GWD_AFR_5mer_counts, 
                               TSI_EUR_5mer_counts, IBS_EUR_5mer_counts, GBR_EUR_5mer_counts, CEU_EUR_5mer_counts, FIN_EUR_5mer_counts,
                               PJL_SAS_5mer_counts, GIH_SAS_5mer_counts, ITU_SAS_5mer_counts, STU_SAS_5mer_counts, BEB_SAS_5mer_counts, 
                               CDX_EAS_5mer_counts, KHV_EAS_5mer_counts, CHS_EAS_5mer_counts, CHB_EAS_5mer_counts, JPT_EAS_5mer_counts) 


pops.ancestral <- c("LWK", "ESN", "YRI", "MSL", "GWD", 
                    "TSI", "IBS", "GBR", "CEU", "FIN", 
                    "PJL", "GIH", "ITU", "STU", "BEB",
                    "CDX", "KHV", "CHS", "CHB", "JPT")

subpops.5mer.all <- list(LWK_AFR_5mer_counts, ESN_AFR_5mer_counts, YRI_AFR_5mer_counts, MSL_AFR_5mer_counts, GWD_AFR_5mer_counts, 
                         ACB_AMR_5mer_counts, ASW_AMR_5mer_counts, CLM_AMR_5mer_counts, MXL_AMR_5mer_counts, PUR_AMR_5mer_counts, PEL_AMR_5mer_counts,
                         TSI_EUR_5mer_counts, IBS_EUR_5mer_counts, GBR_EUR_5mer_counts, CEU_EUR_5mer_counts, FIN_EUR_5mer_counts,
                         PJL_SAS_5mer_counts, GIH_SAS_5mer_counts, ITU_SAS_5mer_counts, STU_SAS_5mer_counts, BEB_SAS_5mer_counts, 
                         CDX_EAS_5mer_counts, KHV_EAS_5mer_counts, CHS_EAS_5mer_counts, CHB_EAS_5mer_counts, JPT_EAS_5mer_counts) 
                         

pops.all <- c("LWK", "ESN", "YRI", "MSL", "GWD", 
              "ACB", "ASW", "CLM", "MXL", "PUR", "PEL",
              "TSI", "IBS", "GBR", "CEU", "FIN", 
              "PJL", "GIH", "ITU", "STU", "BEB",
              "CDX", "KHV", "CHS", "CHB", "JPT")

# X5mer_subpop_rates <- data.frame(matrix(0, nrow = 1536, ncol = 20))
#  
# for (i in 1:20){
#   X5mer_subpop_rates[,i] <- subpops.5mer[[i]]$Rate
#  }
#  
# write.table(X5mer_subpop_rates, "X5mer_subpop_rates", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
