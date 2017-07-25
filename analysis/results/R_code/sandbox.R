
for (i in 1:32){
  result$One_mer[i] <- substr(result$Context[i],2,2) 
}

chi <- -2*(sum(df$likelihood_0)-sum(df$likelihood_1))
pchisq(chi, df = 95, lower.tail = FALSE, log.p = TRUE)*log10(exp(1))
