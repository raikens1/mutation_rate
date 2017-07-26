#This file contains all the steps I used to run an FDR analysis
#library(qvalue)

#ps.7mer is a dataframe output from Fourway.Chi
ps.7mer <- fourway.chi(AFR_7mer_counts, EUR_7mer_counts, EAS_7mer_counts, SAS_7mer_counts)

ps.7mer <- ps.7mer[complete.cases(ps.7mer),]

#need to check range of p values to set lambda
range(ps.7mer$p) # 0 - 0.4

l <- seq(0.05, 0.1, length.out = 30)

#This uses Story et al.'s local fdr method
#ps.7mer$fdr <- lfdr(ps.7mer$p, lambda = l)

#This uses Benjamini-Hoochberg, which I'm favoring because I'm concerned by selecting lambda
ps.7mer$fdr <- p.adjust(ps.7mer$p, method = "fdr")
ps.7mer$holm <- p.adjust(ps.7mer$p, method = "holm")

#qs.7mer <- qvalue(ps.7mer$p, lambda = l)

#plot(qs.7mer)

alpha = 0.05/length(ps.7mer$p)
X7mer <- ps.7mer[complete.cases(ps.7mer), ]

sig <- c(sum(X7mer$p < alpha), sum(X7mer$holm< 0.05), sum(X7mer$fdr< 0.1), sum(X7mer$fdr< 0.05), sum(X7mer$fdr< 0.01), sum(X7mer$fdr< 0.001))
 
names(sig) <- c("Bonferroni", "Holm", "FDR<0.1", "FDR<0.05", "FDR<0.01", "FDR<0.001")
