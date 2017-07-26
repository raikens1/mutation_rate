library(readr)

# just some quick code to help me upload a bunch of subpop files
# meant to be run from the 'data' directory

# EUROPE
TSI_EUR_7mer_counts <- read_delim("subpops/7mer/TSI_EUR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

IBS_EUR_7mer_counts <- read_delim("subpops/7mer/IBS_EUR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CEU_EUR_7mer_counts <- read_delim("subpops/7mer/CEU_EUR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GBR_EUR_7mer_counts <- read_delim("subpops/7mer/GBR_EUR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

FIN_EUR_7mer_counts <- read_delim("subpops/7mer/FIN_EUR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# SOUTH ASIA
STU_SAS_7mer_counts <- read_delim("subpops/7mer/STU_SAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ITU_SAS_7mer_counts <- read_delim("subpops/7mer/ITU_SAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

BEB_SAS_7mer_counts <- read_delim("subpops/7mer/BEB_SAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GIH_SAS_7mer_counts <- read_delim("subpops/7mer/GIH_SAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PJL_SAS_7mer_counts <- read_delim("subpops/7mer/PJL_SAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# AFRICA
ESN_AFR_7mer_counts <- read_delim("subpops/7mer/ESN_AFR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MSL_AFR_7mer_counts <- read_delim("subpops/7mer/MSL_AFR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GWD_AFR_7mer_counts <- read_delim("subpops/7mer/GWD_AFR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

LWK_AFR_7mer_counts <- read_delim("subpops/7mer/LWK_AFR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

YRI_AFR_7mer_counts <- read_delim("subpops/7mer/YRI_AFR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# EAST ASIA
CHB_EAS_7mer_counts <- read_delim("subpops/7mer/CHB_EAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

JPT_EAS_7mer_counts <- read_delim("subpops/7mer/JPT_EAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CHS_EAS_7mer_counts <- read_delim("subpops/7mer/CHS_EAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CDX_EAS_7mer_counts <- read_delim("subpops/7mer/CDX_EAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

KHV_EAS_7mer_counts <- read_delim("subpops/7mer/KHV_EAS_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# AMERICA

ACB_AMR_7mer_counts <- read_delim("subpops/7mer/ACB_AMR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ASW_AMR_7mer_counts <- read_delim("subpops/7mer/ASW_AMR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CLM_AMR_7mer_counts <- read_delim("subpops/7mer/CLM_AMR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MXL_AMR_7mer_counts <- read_delim("subpops/7mer/MXL_AMR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PEL_AMR_7mer_counts <- read_delim("subpops/7mer/PEL_AMR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PUR_AMR_7mer_counts <- read_delim("subpops/7mer/PUR_AMR_7mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

subpops.7mer.ancestral <- list(LWK_AFR_7mer_counts, ESN_AFR_7mer_counts, YRI_AFR_7mer_counts, MSL_AFR_7mer_counts, GWD_AFR_7mer_counts, 
                               TSI_EUR_7mer_counts, IBS_EUR_7mer_counts, GBR_EUR_7mer_counts, CEU_EUR_7mer_counts, FIN_EUR_7mer_counts,
                               PJL_SAS_7mer_counts, GIH_SAS_7mer_counts, ITU_SAS_7mer_counts, STU_SAS_7mer_counts, BEB_SAS_7mer_counts, 
                               CDX_EAS_7mer_counts, KHV_EAS_7mer_counts, CHS_EAS_7mer_counts, CHB_EAS_7mer_counts, JPT_EAS_7mer_counts) 


pops.ancestral <- c("LWK", "ESN", "YRI", "MSL", "GWD", 
                    "TSI", "IBS", "GBR", "CEU", "FIN", 
                    "PJL", "GIH", "ITU", "STU", "BEB",
                    "CDX", "KHV", "CHS", "CHB", "JPT")

subpops.7mer.all <- list(LWK_AFR_7mer_counts, ESN_AFR_7mer_counts, YRI_AFR_7mer_counts, MSL_AFR_7mer_counts, GWD_AFR_7mer_counts, 
                         ACB_AMR_7mer_counts, ASW_AMR_7mer_counts, CLM_AMR_7mer_counts, MXL_AMR_7mer_counts, PUR_AMR_7mer_counts, PEL_AMR_7mer_counts,
                         TSI_EUR_7mer_counts, IBS_EUR_7mer_counts, GBR_EUR_7mer_counts, CEU_EUR_7mer_counts, FIN_EUR_7mer_counts,
                         PJL_SAS_7mer_counts, GIH_SAS_7mer_counts, ITU_SAS_7mer_counts, STU_SAS_7mer_counts, BEB_SAS_7mer_counts, 
                         CDX_EAS_7mer_counts, KHV_EAS_7mer_counts, CHS_EAS_7mer_counts, CHB_EAS_7mer_counts, JPT_EAS_7mer_counts) 
                         

pops.all <- c("LWK", "ESN", "YRI", "MSL", "GWD", 
              "ACB", "ASW", "CLM", "MXL", "PUR", "PEL",
              "TSI", "IBS", "GBR", "CEU", "FIN", 
              "PJL", "GIH", "ITU", "STU", "BEB",
              "CDX", "KHV", "CHS", "CHB", "JPT")

# X7mer_subpop_rates <- data.frame(matrix(0, nrow = 1536, ncol = 20))
#  
# for (i in 1:20){
#   X7mer_subpop_rates[,i] <- subpops.7mer[[i]]$Rate
#  }
#  
# write.table(X7mer_subpop_rates, "X7mer_subpop_rates", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
