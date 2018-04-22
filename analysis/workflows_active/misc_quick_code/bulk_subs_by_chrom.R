
# Given a 3mer counts file, plot bulk private substitution probability by chromosome
rate_by_chrom <- function(Counts, pop = "My Population"){
  subsSums <- colSums(Counts[,-c(1:4)])
  subsSums[1] <- subsSums[1]/3.0
  gwSums <- colSums(gw_3mer_counts[,-c(1,2)])
  plot(subsSums[-1]/gwSums, xlab = "Chromosome", ylab = "Private substitution probability", main = pop)
}
par(mfrow=c(2,2))
rate_by_chrom(AFR_3mer_counts, "AFR")
rate_by_chrom(EUR_3mer_counts, "EUR")
rate_by_chrom(EAS_3mer_counts, "EAS")
rate_by_chrom(SAS_3mer_counts, "SAS")
