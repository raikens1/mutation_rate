library(DescTools)

## CHISQ_TESTS
##################
# pairwise.chi
#   given: two dfs of count data for all polymorphism types
#   do: return a df of p values by context type
# fourway.chi
#   given: four dfs of count data for all polymorphism types
#   do: return a df of p values by context type
# x.chi
#   given: a counts dataset (POP), and a gw contexts dataset(gw),
#   do: return a df of p values for enrichment on X using a Chi squared test
# get.sigs
#   given:a p value df from a chi squared test function
#   do: return a dataframe of the significant contexts (using bonferroni correction)
##################

#calculates homogeneity test p values for pairwise comparisons of two dfs of counts
pairwise.chi <- function(counts.1, counts.2){
  n.contexts = length(counts.1$Context)
  result <- data.frame(matrix(ncol=2,nrow=n.contexts))
  colnames(result) <- c("Context", "p")
  
  result$Context <- counts.1$Context
  sum.1 <- sum(counts.1$Count)
  sum.2 <- sum(counts.2$Count)

  for (i in 1:n.contexts){
    c.a <- c(counts.1$Count[i], counts.2$Count[i])
    c.b <- c(sum.1, sum.2) - c.a
    data <- cbind(c.a, c.b)
    result$p[i] <- chisq.test(data)$p.value #continuity correction?
  }
  
  return(result)
}


################################################

#calculates homogeneity test p values for Fourway comparisons of counts dfs
fourway.chi <- function(AFR, EUR, EAS, SAS){
  n.contexts = length(AFR$Context)
  
  #make dataframe for results
  result <- data.frame(matrix(ncol=9,nrow=n.contexts))
  colnames(result) <- c("Context", "X5mer","X3mer", "X1mer", "AFR.Count", "EUR.Count", "EAS.Count", "SAS.Count", "p")
  result$Context <- AFR$Context
  result$X5mer <- AFR$X5mer
  result$X3mer <- AFR$X3mer
  result$X1mer <- AFR$X1mer
  result$AFR.Count <- AFR$Count
  result$EUR.Count <- EUR$Count
  result$EAS.Count <- EAS$Count
  result$SAS.Count <- SAS$Count
  
  #start setting up tables
  sum.1 <- sum(AFR$Count)
  sum.2 <- sum(EUR$Count)
  sum.3 <- sum(EAS$Count)
  sum.4 <- sum(SAS$Count)
  
  #set up table and run test for each context
  for (i in 1:n.contexts){
    c.a <- c(AFR$Count[i], EUR$Count[i], EAS$Count[i], SAS$Count[i])
    c.b <- c(sum.1, sum.2, sum.3, sum.4) - c.a
    data <- cbind(c.a, c.b)
    result$p[i] <- chisq.test(data)$p.value #continuity correction?
  }
  
  return(result)
}

#retrieve significant subset of results from chi squared test
get.sigs <- function(ps){
  alpha <- 0.05/length(ps$Context)
  return(subset(ps, ps$p < alpha))
}

#############################################

#Given a counts dataset (POP), and a gw contexts dataset(gw),
#returns a df of p values for enrichment on X using a Chi squared test
x.chi <- function(POP, gw){
  #set up result table
  result <- data.frame(matrix(ncol=6,nrow=length(POP$Context)))
  colnames(result) <- c("Context", "Poly on X", "Sites on X", "Poly on Auto", "Sites on Auto","p")
  result$Context <- POP$Context
  
  #for each context
  for (i in 1:length(gw$GW_total)){
    i.1 <- 1+3*(i-1)
    
    #for each mutation possible from the ith context
    for (j in 1:3){
      result$"Sites on X"[i.1] <- gw$X[i]
      result$"Sites on Auto"[i.1] <- sum(gw[i, 3:24])
      result$"Poly on X"[i.1] <- POP$chrX[i.1]
      result$"Poly on Auto"[i.1] <- sum(POP[i.1, 8:29])
      
      #run chisquared test and fill out table
      data <- cbind(as.numeric(c(result[i.1,2:3])), as.numeric(c(result[i.1, 4:5])))
      result$p[i.1] <- chisq.test(as.matrix(data))$p.value
      i.1 <- i.1+1
    }
  }
  
  return(result)
}
