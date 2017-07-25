EASnotJPT <- KHV_EAS_7mer_counts[,-c(1, 25)] + CHB_EAS_7mer_counts[,-c(1, 25)] + 
              CHS_EAS_7mer_counts[,-c(1, 25)] + CDX_EAS_7mer_counts[,-c(1, 25)]
EASnotJPTchrom_counts <- cbind(JPT_EAS_7mer_counts$Context, EASnotJPT, JPT_EAS_7mer_counts$X1mer)
colnames(EASnotJPTchrom_counts) <- colnames(KHV_EAS_7mer_counts)
EASnotJPT_7mer_counts <- process_chrom_counts(EASnotJPTchrom_counts, gw_7mer_counts, X7mer_mutations_ref)
write.table(EASnotJPT_7mer_counts, "EASnotJPT_7mer_counts.txt", quote = F, sep = "\t", row.names = F)
