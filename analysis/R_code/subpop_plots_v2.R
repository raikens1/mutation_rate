library(ggplot2)
library(dplyr)
require(reshape2)

colors.ancestral <- c(rep("forestgreen", 5), rep("darkblue", 5), rep("magenta", 5), rep("red", 5))
colors.all <- c(rep("forestgreen", 5), rep("orange", 6), rep("darkblue", 5), rep("magenta", 5), rep("red", 5))

count.table.subpop <- function(subpops, kmer, popnames){
  a <- (kmer+1)/2
  keep <- c(1:a, a+3)
  ct <- subpops[[1]][,keep]
  for (i in 1:length(subpops)){
    ct[a+1 + i] = subpops[[i]][,a+1]
  }
  colnames(ct)[-c(1:(a+1))] <- popnames
  return(ct)
}

CI.plot.subpop <- function(subpops, mut, popnames, colors = colors.ancestral){
  n.pops <- length(subpops)
  mut.i <- which(subpops[[1]]$Context == mut)
  
  #have to do a silly workaround or R will sort popnames alphanumerically
  poplabs <- factor(popnames, levels= popnames)
  
  #cycle through subpops and get counts for mut
  counts <- rep(0, n.pops)
  sums <- rep(0, n.pops)
  for (i in 1:n.pops){
    counts[i] <- subpops[[i]]$Count[mut.i]
    sums[i]<- sum(subpops[[i]]$Count)
  }
  
  
  #estimate substitution probability
  N.c <- subpops[[1]]$context_in_genome[mut.i] 
  theta <- counts/N.c
  L <- theta - 1.96*sqrt(theta*(1-theta)/N.c)
  U <- theta + 1.96*sqrt(theta*(1-theta)/N.c)
  
  #normalize to rate estimate; assume genome wide subsitution probability is measured without error
  norm <- 1.2E-8*(sum(as.numeric(subpops[[1]]$context_in_genome))/3)/sums
  
  
  df <- data.frame(cbind(poplabs, theta*norm, L*norm, U*norm))
  plotcol <- reorder(colors, theta*norm)
  
  ggplot(df, aes(reorder(popnames, poplabs), theta*norm)) +
    geom_point(size = 4, color = plotcol) +
    geom_errorbar(aes(ymax = U*norm, ymin = L*norm), color = plotcol)+
    labs(y = paste("Estimated mutation rate of ", mut,"\n")) + #y axis label
    theme(axis.text.x = element_text(size = rel(1.4), angle = 45, hjust = 1), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)),
          legend.position = 'none')
}
