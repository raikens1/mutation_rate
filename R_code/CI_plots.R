library(ggplot2)
library(dplyr)
require(reshape2)

## Superpop_CIs
#---------------------------------------------
# CI.plot
#   plots approximate CIs for superpopulations
#---------------------------------------------

CI.plot <- function(AFR, EUR, EAS, SAS, mut) {
  
  mut.i <- which(AFR$Context == mut)
  popnames <- c("European", "South\nAsian", "African", "East\nAsian")
  colors <- c("darkblue", "magenta", "forestgreen", "red")
  
  #have to do a silly workaround or R will sort popnames alphanumerically
  poplabs <- factor(popnames, levels= popnames)
  
  #cycle through pops and get counts for mut
  counts <- rep(0, 4)
  sums <- rep(0, 4)
  pops <- list(EUR, SAS, AFR, EAS)
  
  for (i in 1:length(pops)){
    counts[i] <- pops[[i]]$Count[mut.i]
    sums[i]<- sum(pops[[i]]$Count)
  }
  
  #estimate substitution probability
  N.c <- AFR$context_in_genome[mut.i] 
  theta <- counts/N.c
  L <- theta - 1.96*sqrt(theta*(1-theta)/N.c)
  U <- theta + 1.96*sqrt(theta*(1-theta)/N.c)
  
  #normalize to rate estimate; assume genome wide subsitution probability is measured without error
  norm <- 1.2E-8*(sum(as.numeric(AFR$context_in_genome))/3)/sums
  
  df <- data.frame(cbind(popnames, theta*norm, L*norm, U*norm))
  plotcol <- reorder(colors, theta*norm)
  
  CI.plot <- ggplot(df, aes(reorder(popnames, c(1,2,3,4)), theta*norm)) +
    geom_point(size = 4, color = plotcol) +
    geom_errorbar(aes(ymax = U*norm, ymin = L*norm), color = plotcol, size =1)+
    labs(y = paste("Estimated mutation rate of ", mut,"\n")) + #y axis label
    theme(axis.text.x = element_text(size = rel(1.4)), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)),
          legend.position = 'none')
  return(CI.plot)

}
