library(readr)

# just some quick code to help me upload a bunch of subpop files
# meant to be run from the 'data' directory

# EUROPE
TSI_EUR_3mer_counts <- read_delim("subpops/3mer/TSI_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

IBS_EUR_3mer_counts <- read_delim("subpops/3mer/IBS_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CEU_EUR_3mer_counts <- read_delim("subpops/3mer/CEU_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GBR_EUR_3mer_counts <- read_delim("subpops/3mer/GBR_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

FIN_EUR_3mer_counts <- read_delim("subpops/3mer/FIN_EUR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# SOUTH ASIA
STU_SAS_3mer_counts <- read_delim("subpops/3mer/STU_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ITU_SAS_3mer_counts <- read_delim("subpops/3mer/ITU_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

BEB_SAS_3mer_counts <- read_delim("subpops/3mer/BEB_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GIH_SAS_3mer_counts <- read_delim("subpops/3mer/GIH_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PJL_SAS_3mer_counts <- read_delim("subpops/3mer/PJL_SAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# AFRICA
ESN_AFR_3mer_counts <- read_delim("subpops/3mer/ESN_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MSL_AFR_3mer_counts <- read_delim("subpops/3mer/MSL_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

GWD_AFR_3mer_counts <- read_delim("subpops/3mer/GWD_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

LWK_AFR_3mer_counts <- read_delim("subpops/3mer/LWK_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

YRI_AFR_3mer_counts <- read_delim("subpops/3mer/YRI_AFR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# EAST ASIA
CHB_EAS_3mer_counts <- read_delim("subpops/3mer/CHB_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

JPT_EAS_3mer_counts <- read_delim("subpops/3mer/JPT_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CHS_EAS_3mer_counts <- read_delim("subpops/3mer/CHS_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CDX_EAS_3mer_counts <- read_delim("subpops/3mer/CDX_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

KHV_EAS_3mer_counts <- read_delim("subpops/3mer/KHV_EAS_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# AMERICA

ACB_AMR_3mer_counts <- read_delim("subpops/3mer/ACB_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ASW_AMR_3mer_counts <- read_delim("subpops/3mer/ASW_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

CLM_AMR_3mer_counts <- read_delim("subpops/3mer/CLM_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

MXL_AMR_3mer_counts <- read_delim("subpops/3mer/MXL_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PEL_AMR_3mer_counts <- read_delim("subpops/3mer/PEL_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

PUR_AMR_3mer_counts <- read_delim("subpops/3mer/PUR_AMR_3mer_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

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

# Code for making rate matricies
# Really should be elsewhere
#
# subpops.names <- c("TSI", "IBS", "CEU", "GBR", "FIN", "STU", "ITU", "BEB", "GIH", "PJL",
#                    "ESN", "GWD", "LWK", "MSL", "YRI", "CDX", "CHB", "CHS", "JPT", "KHV")
# 
# subpops.3mer <- list(TSI_EUR_3mer_counts, IBS_EUR_3mer_counts, CEU_EUR_3mer_counts, GBR_EUR_3mer_counts, FIN_EUR_3mer_counts, 
#                      STU_SAS_3mer_counts, ITU_SAS_3mer_counts, BEB_SAS_3mer_counts, GIH_SAS_3mer_counts, PJL_SAS_3mer_counts,
#                      ESN_AFR_3mer_counts, GWD_AFR_3mer_counts, LWK_AFR_3mer_counts, MSL_AFR_3mer_counts, YRI_AFR_3mer_counts,
#                      CDX_EAS_3mer_counts, CHB_EAS_3mer_counts, CHS_EAS_3mer_counts, JPT_EAS_3mer_counts, KHV_EAS_3mer_counts)
#
# X3mer_subpop_rates <- matrix(0, nrow = 96, ncol = 20)
#   
#   for (i in 1:20){
#     X3mer_subpop_rates[,i] <- subpops.3mer[[i]]$Rate
#   }
#
# colnames(X3mer_subpop_rates) <- subpops.names
# rownames(X3mer_subpop_rates) <- TSI_EUR_3mer_counts$Context
#   
# write.table(X3mer_subpop_rates, "rate_profiles/rates_3mer.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = TRUE)
