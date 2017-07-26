library(ggplot2)

CI.plot.bygroup <- function(AFR, EUR, EAS, SAS, muts, groupname = "mutation group") {
  #NOTE: muts must all be from different contexts.  Otherwise this will do weird things.
  
  muts.i <- which(is.element(AFR$Context, muts))
  popnames <- c("Africa","Europe", "South\nAsia", "East\nAsia")
  colors <- c("forestgreen", "darkblue", "magenta","red")
  
  #have to do a silly workaround or R will sort popnames alphanumerically
  poplabs <- factor(popnames, levels= popnames)
  
  #cycle through pops and get counts for mut
  counts <- rep(0, 4)
  sums <- rep(0, 4)
  pops <- list(AFR, EUR, SAS, EAS)
  
  for (i in 1:length(pops)){
    counts[i] <- sum(pops[[i]]$Count[muts.i]) #number of obersvations of muts in pop
    sums[i]<- sum(pops[[i]]$Count)#total polymorphisms in pop
  }
  
  #estimate substitution probability
  N.c <- sum(AFR$context_in_genome[muts.i])
  theta <- counts/N.c
  L <- theta - 1.96*sqrt(theta*(1-theta)/N.c)
  U <- theta + 1.96*sqrt(theta*(1-theta)/N.c)
  
  #normalize to rate estimate; assume genome wide subsitution probability is measured without error
  norm <- 1.2E-8*(sum(as.numeric(AFR$context_in_genome))/3)/sums
  
  df <- data.frame(cbind(popnames, theta*norm, L*norm, U*norm))
  plotcol <- reorder(colors, theta*norm)
  
  CI.plot <- ggplot(df, aes(reorder(popnames, c(1,2,3,4)), theta*norm)) +
    geom_point(size = 3, color = plotcol) +
    geom_errorbar(aes(ymax = U*norm, ymin = L*norm), color = plotcol, size =.75)+
    labs(title = paste("Estimated mutation\nrate of", groupname), y = NULL) + #y axis label
    theme(axis.text.x = element_text(size = rel(.9)), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(.9)), axis.text.y = element_text(size = rel(.9), angle = 0, hjust = 0.5), title = element_text(size = rel(.7)),
          legend.position = 'none')
  return(CI.plot)
}


colors.ancestral <- c(rep("forestgreen", 5), rep("darkblue", 5), 
                      rep("magenta", 5), rep("red", 5))
colors.all <- c(rep("forestgreen", 5), rep("orange", 6), 
                rep("darkblue", 5), rep("magenta", 5), rep("red", 5))


CI.plot.subpop.bygroup <- function(subpops, muts, popnames, colors = colors.ancestral, groupname = "mutation group"){
  n.pops <- length(subpops)
  mut.i <- which(is.element(subpops[[1]]$Context, muts))
  
  #have to do a silly workaround or R will sort popnames alphanumerically
  poplabs <- factor(popnames, levels= popnames)
  
  #cycle through subpops and get counts for mut
  counts <- rep(0, n.pops)
  sums <- rep(0, n.pops)
  for (i in 1:n.pops){
    counts[i] <- sum(subpops[[i]]$Count[mut.i])
    sums[i]<- sum(subpops[[i]]$Count)
  }
  
  #estimate substitution probability
  N.c <- sum(subpops[[1]]$context_in_genome[mut.i]) 
  theta <- counts/N.c
  L <- theta - 1.96*sqrt(theta*(1-theta)/N.c)
  U <- theta + 1.96*sqrt(theta*(1-theta)/N.c)
  
  #normalize to rate estimate; assume genome wide subsitution probability is measured without error
  norm <- 1.2E-8*(sum(as.numeric(subpops[[1]]$context_in_genome))/3)/sums
  
  
  df <- data.frame(cbind(poplabs, theta*norm, L*norm, U*norm))
  plotcol <- reorder(colors, theta*norm)
  
  ggplot(df, aes(reorder(popnames, poplabs), theta*norm)) +
    geom_point(size = 3, color = plotcol) +
    geom_errorbar(aes(ymax = U*norm, ymin = L*norm), color = plotcol, size = .75)+
    labs(y = paste("Estimated mutation rate of ", groupname)) + #y axis label
    theme(axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(.8)), axis.text.y = element_text(size = rel(.9)),
          legend.position = 'none')
}