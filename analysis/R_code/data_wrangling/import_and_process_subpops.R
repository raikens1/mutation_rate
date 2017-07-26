library(readr)
source('C:/Users/VoightLab/Dropbox/SNP_rates/Code/data_wrangling/process_chrom_counts.R')
X3mer_mutations_ref <- read_delim("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/3mer_mutations_ref.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
gw_3mer_counts <- read_delim("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/gw_counts/gw_3mer_counts.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

setwd("C:/Users/VoightLab/Dropbox/SNP_rates/Raw_Data/singletons_excluded/subpops/3mer/")

#just some quick code to help me upload and process a bunch of subpop files.

#EUROPE
TSI_EUR_3mer_counts <- read_delim("chrom_counts/TSI_EUR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
TSI_EUR_3mer_counts <- process_chrom_counts(TSI_EUR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = TSI_EUR_3mer_counts, "TSI_EUR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

IBS_EUR_3mer_counts <- read_delim("chrom_counts/IBS_EUR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
IBS_EUR_3mer_counts <- process_chrom_counts(IBS_EUR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = IBS_EUR_3mer_counts, "IBS_EUR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

CEU_EUR_3mer_counts <- read_delim("chrom_counts/CEU_EUR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
CEU_EUR_3mer_counts <- process_chrom_counts(CEU_EUR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = CEU_EUR_3mer_counts, "CEU_EUR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

GBR_EUR_3mer_counts <- read_delim("chrom_counts/GBR_EUR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GBR_EUR_3mer_counts <- process_chrom_counts(GBR_EUR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = GBR_EUR_3mer_counts, "GBR_EUR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

FIN_EUR_3mer_counts <- read_delim("chrom_counts/FIN_EUR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
FIN_EUR_3mer_counts <- process_chrom_counts(FIN_EUR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = FIN_EUR_3mer_counts, "FIN_EUR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#SOUTH ASIA
STU_SAS_3mer_counts <- read_delim("chrom_counts/STU_SAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
STU_SAS_3mer_counts <- process_chrom_counts(STU_SAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = STU_SAS_3mer_counts, "STU_SAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

ITU_SAS_3mer_counts <- read_delim("chrom_counts/ITU_SAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ITU_SAS_3mer_counts <- process_chrom_counts(ITU_SAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = ITU_SAS_3mer_counts, "ITU_SAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

BEB_SAS_3mer_counts <- read_delim("chrom_counts/BEB_SAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
BEB_SAS_3mer_counts <- process_chrom_counts(BEB_SAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = BEB_SAS_3mer_counts, "BEB_SAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

GIH_SAS_3mer_counts <- read_delim("chrom_counts/GIH_SAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GIH_SAS_3mer_counts <- process_chrom_counts(GIH_SAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = GIH_SAS_3mer_counts, "GIH_SAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

PJL_SAS_3mer_counts <- read_delim("chrom_counts/PJL_SAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
PJL_SAS_3mer_counts <- process_chrom_counts(PJL_SAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = PJL_SAS_3mer_counts, "PJL_SAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#AFRICA
ESN_AFR_3mer_counts <- read_delim("chrom_counts/ESN_AFR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ESN_AFR_3mer_counts <- process_chrom_counts(ESN_AFR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = ESN_AFR_3mer_counts, "ESN_AFR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

MSL_AFR_3mer_counts <- read_delim("chrom_counts/MSL_AFR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
MSL_AFR_3mer_counts <- process_chrom_counts(MSL_AFR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = MSL_AFR_3mer_counts, "MSL_AFR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

GWD_AFR_3mer_counts <- read_delim("chrom_counts/GWD_AFR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GWD_AFR_3mer_counts <- process_chrom_counts(GWD_AFR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = GWD_AFR_3mer_counts, "GWD_AFR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

LWK_AFR_3mer_counts <- read_delim("chrom_counts/LWK_AFR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
LWK_AFR_3mer_counts <- process_chrom_counts(LWK_AFR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = LWK_AFR_3mer_counts, "LWK_AFR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

YRI_AFR_3mer_counts <- read_delim("chrom_counts/YRI_AFR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
YRI_AFR_3mer_counts <- process_chrom_counts(YRI_AFR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = YRI_AFR_3mer_counts, "YRI_AFR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#EAST ASIA
CHB_EAS_3mer_counts <- read_delim("chrom_counts/CHB_EAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
CHB_EAS_3mer_counts <- process_chrom_counts(CHB_EAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = CHB_EAS_3mer_counts, "CHB_EAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

JPT_EAS_3mer_counts <- read_delim("chrom_counts/JPT_EAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
JPT_EAS_3mer_counts <- process_chrom_counts(JPT_EAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = JPT_EAS_3mer_counts, "JPT_EAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

CHS_EAS_3mer_counts <- read_delim("chrom_counts/CHS_EAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
CHS_EAS_3mer_counts <- process_chrom_counts(CHS_EAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = CHS_EAS_3mer_counts, "CHS_EAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

CDX_EAS_3mer_counts <- read_delim("chrom_counts/CDX_EAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
CDX_EAS_3mer_counts <- process_chrom_counts(CDX_EAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = CDX_EAS_3mer_counts, "CDX_EAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

KHV_EAS_3mer_counts <- read_delim("chrom_counts/KHV_EAS_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
KHV_EAS_3mer_counts <- process_chrom_counts(KHV_EAS_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = KHV_EAS_3mer_counts, "KHV_EAS_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#AMERICA
ACB_AMR_3mer_counts <- read_delim("chrom_counts/ACB_AMR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ACB_AMR_3mer_counts <- process_chrom_counts(ACB_AMR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = ACB_AMR_3mer_counts, "ACB_AMR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

ASW_AMR_3mer_counts <- read_delim("chrom_counts/ASW_AMR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ASW_AMR_3mer_counts <- process_chrom_counts(ASW_AMR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = ASW_AMR_3mer_counts, "ASW_AMR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

CLM_AMR_3mer_counts <- read_delim("chrom_counts/CLM_AMR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
CLM_AMR_3mer_counts <- process_chrom_counts(CLM_AMR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = CLM_AMR_3mer_counts, "CLM_AMR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

MXL_AMR_3mer_counts <- read_delim("chrom_counts/MXL_AMR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
MXL_AMR_3mer_counts <- process_chrom_counts(MXL_AMR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = MXL_AMR_3mer_counts, "MXL_AMR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

PEL_AMR_3mer_counts <- read_delim("chrom_counts/PEL_AMR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
PEL_AMR_3mer_counts <- process_chrom_counts(PEL_AMR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = PEL_AMR_3mer_counts, "PEL_AMR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

PUR_AMR_3mer_counts <- read_delim("chrom_counts/PUR_AMR_chrom_counts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
PUR_AMR_3mer_counts <- process_chrom_counts(PUR_AMR_3mer_counts, gw_3mer_counts, subcontext_ref = X3mer_mutations_ref)
write.table(x = PUR_AMR_3mer_counts, "PUR_AMR_3mer_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

subpops.3mer <- list(TSI_EUR_3mer_counts, IBS_EUR_3mer_counts, CEU_EUR_3mer_counts, GBR_EUR_3mer_counts, FIN_EUR_3mer_counts, STU_SAS_3mer_counts, ITU_SAS_3mer_counts, BEB_SAS_3mer_counts, GIH_SAS_3mer_counts, PJL_SAS_3mer_counts, ESN_AFR_3mer_counts, GWD_AFR_3mer_counts, LWK_AFR_3mer_counts, MSL_AFR_3mer_counts, YRI_AFR_3mer_counts, CDX_EAS_3mer_counts, CHB_EAS_3mer_counts, CHS_EAS_3mer_counts, JPT_EAS_3mer_counts, KHV_EAS_3mer_counts)

