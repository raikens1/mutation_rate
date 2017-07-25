#Make supplementary table for Pairwise comparison
Pairwise_table <- function(counts.1, counts.2){
  n.contexts = length(counts.1$Context)
  result <- data.frame(matrix(ncol=8,nrow=n.contexts))
  colnames(result) <- c("Context", "One_mer", "m.p1", "m.p2", "m.compliment.p1", "m.compliment.p2", "chi.stat", "p")
  
  result$Context <- counts.1$Context
  result$One_mer <- counts.1$One_mer
  sum.1 <- sum(counts.1$Count)
  sum.2 <- sum(counts.2$Count)
  
  for (i in 1:n.contexts){
    c.a <- c(counts.1$Count[i], counts.2$Count[i])
    c.b <- c(sum.1, sum.2) - c.a
    data <- cbind(c.a, c.b)
    result$m.p1[i] <- c.a[1]
    result$m.p2[i] <- c.a[2]
    result$m.compliment.p1[i] <- c.b[1]
    result$m.compliment.p2[i] <- c.b[2]
    chi.test <- chisq.test(data) #continuity correction?
    result$chi.stat[i] <- chi.test$statistic
    result$p[i] <- chi.test$p.value
  }
  
  return(result)
}



addtotals <- function(df){
  df$Counts = rep(0, length(df$Context))
  for (i in 1:length(df$Context)){
    df$Counts[i] <- sum(df[i,][-c(1,25)])
  }
  colnames(df)<- c("Context", "chr1", "chr2", 'chr3', 'chr4', "chr5", "chr6", 'chr7', 'chr8', "chr9", "chr10", 'chr11', 'chr12', "chr13", "chr14", 'chr15', 'chr16', "chr17", "chr18", 'chr19', 'chr20', "chr21", "chr22", 'chrX', 'One_mer', 'Count')
  return(df)
}

