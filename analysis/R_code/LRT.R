likelihood.test <- function(priv, cosmo, gw, nmer){
  if (nmer == 3){
    dfs <- c(31, 64, 95)
  }
  else if (nmer == 5){
    dfs <- c(511, 1024, 1535)
  }
  else {
    dfs <- c(8191, 16384, 24575)
  }
  n.contexts <- length(gw$Context)
  result <- data.frame(matrix(nrow = n.contexts + 1, ncol = 18))
  colnames(result)<- c('Context', 's1', 's2', 's3', 'p1', 'p2', 'p3', 'same', 
                       'alpha', 'likelihood_0','likelihood_1', 'likelihood_2',
                       'Chi_stat_0v1', 'p_0v1', 'Chi_stat_1v2', 'p_1v2', 'Chi_stat_0v2', 'p_0v2')
  result$Context[-(n.contexts+1)] <- as.character(gw$Context)
  alpha.all <- sum(priv$Count)/sum(cosmo$Count)
  
  for (i in 1:n.contexts){
    #priv and cosmo are indexed differently than the result table, so calculate a second index value
    i.1 <- 1+3*(i-1)
    
    #getting counts for mutation types from input dataframes
    counts.s <- c(cosmo$Count[i.1], cosmo$Count[i.1+1],cosmo$Count[i.1+2])
    counts.p <- c(priv$Count[i.1],priv$Count[i.1+1],priv$Count[i.1+2])
    gw_total <- gw$GW_total[i]
    counts.all <- c(counts.s, counts.p, gw_total - sum(c(counts.s, counts.p)))
    
    #ratio of private to shared SNPs
    alpha <- sum(counts.p)/sum(counts.s)
    
    #parameters for null and alternative hypotheses
    thetas.H2 <- counts.all/gw_total
    thetas.H1 <- c(counts.s, counts.s*alpha, counts.all[7])/gw_total
    thetas.H0 <- c(counts.s, counts.s*alpha.all, counts.all[7])/gw_total
    
    #log likelihood ratio test
    likelihood_2 <- dmultinom(x = counts.all, prob = thetas.H2, log = TRUE)    
    likelihood_1 <- dmultinom(x = counts.all, prob = thetas.H1, log = TRUE)
    likelihood_0 <- dmultinom(x = counts.all, prob = thetas.H0, log = TRUE)
    Chi_0v1 <- NA
    Chi_1v2 <- -2*(likelihood_1-likelihood_2)
    Chi_0v2 <- NA
    p_0v1 <- NA
    p_1v2 <- pchisq(Chi_1v2, df = 2, lower.tail = FALSE, log.p = FALSE)
    p_0v2 <- NA
    
    #populate result table
    result[i,][-1] <- c( counts.all, alpha, likelihood_0, likelihood_1, likelihood_2,
                        Chi_0v1, p_0v1, Chi_1v2, p_1v2, Chi_0v2, p_0v2)
  }
  
  #calculate stats over all contexts
  likelihood_0.all <- sum(as.numeric(result$likelihood_0), na.rm = TRUE)
  likelihood_1.all <- sum(as.numeric(result$likelihood_1), na.rm = TRUE)
  likelihood_2.all <- sum(as.numeric(result$likelihood_2), na.rm = TRUE)
  Chi_0v1.all <- -2*(likelihood_0.all-likelihood_1.all)
  Chi_1v2.all <- -2*(likelihood_1.all-likelihood_2.all)
  Chi_0v2.all <- -2*(likelihood_0.all-likelihood_2.all)
  p_0v1.all <- pchisq(Chi_0v1.all, df = dfs[1], lower.tail = FALSE, log.p = TRUE)*log10(exp(1))
  p_1v2.all <- pchisq(Chi_1v2.all, df = dfs[2], lower.tail = FALSE, log.p = TRUE)*log10(exp(1))
  p_0v2.all <- pchisq(Chi_0v2.all, df = dfs[3], lower.tail = FALSE, log.p = TRUE)*log10(exp(1))
  
  result[n.contexts + 1,] <- c("Summary", rep(NA, 7), alpha.all, likelihood_0.all, likelihood_1.all, likelihood_2.all,
                                   Chi_0v1.all, p_0v1.all, Chi_1v2.all, p_1v2.all, Chi_0v2.all, p_0v2.all)
  
  return(result)
}

#p-value scattter plot for 4-way comparison
#TODO: debug for new input format
LRT_plot <- function(p.values, lab.lim, before){
  n = length(p.values$Context)
  for (i in 1:n){
    p.values$One_mer[i] <- substr(p.values$Context[i],before+1,before+1) 
  }
  
  alpha = 0.05/n
  p.values$p_1v2 <- as.numeric(p.values$p_1v2)
  data <- p.values[order(p.values$p_1v2),][-1,]
  data$Context <- as.character(data$Context)
  data$Context <- factor(data$Context, levels = unique(data$Context))
  sigs <- subset(data, data$p_1v2 < lab.lim)
  my.plot <- qplot(data$Context, -log10(data$p_1v2)) +
    labs(x = "Context", y = "-Log(p-value)\n")+
    scale_color_manual("", values = c("blue", "red", "red", "blue"))+
    geom_point(aes(color=factor(data$One_mer)), size = 3)+
    geom_text(data = sigs, aes(sigs$Context, -log10(sigs$p_1v2), label = sigs$Context, color = sigs$One_mer), 
              vjust = 1, nudge_x = 25, nudge_y = 0, size = 2, check_overlap = TRUE)+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)),
          legend.text = element_text(size = rel(1.25)), legend.position = c(0.8, 0.8))
  
  return(my.plot)
}
