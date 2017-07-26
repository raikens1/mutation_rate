library(ggplot2)
library(dplyr)
require(reshape2)

## PLOTS BY CHROMOSOME
##################
# chrom.scatter
#   given: dataframes for all pops, gw_context dataframe, mutation type
#   do: return a scatterplot of mutation rate by chromosome
#
# chrom.box
#   given: dataframes for all pops, gw_context dataframe, mutation type
#   do: return a boxplot of mutation rate by chromosome
#
# chrom.process.data (helper function)
#   given: dataframes for all pops, gw_context dataframe, mutation type
#   do: return a dataframe of mutation rate by chromosome for plotting
##################


#scatter plot of rate by chrom for a particular sequence context
chrom.scatter <- function(AFR, EUR, EAS, SAS, gw, mut){
  #get data
  chrom.dat <- chrom.process.data(AFR, EUR, EAS, SAS, gw, mut)
  
  #plot
  c_plot <- ggplot(chrom.dat, aes(chrom, rate))+
    labs(x = "\nChromosome", y = paste("Estimated mutation rate of", mut,"\n"))+
    scale_color_manual("", values = c("forest green", "red", "dark blue", 'magenta'))+
    geom_point(aes(color = factor(chrom.dat$pop)), size = 3)+
    theme(axis.text.x = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.9)), #adjust text sizes
          axis.title.y = element_text(size = rel(1.9)), axis.text.y = element_text(size = rel(1.7)), 
          legend.text = element_text(size = rel(1.7)))#legend position
  return(c_plot)
}

#box plot of rate by chrom for a particular sequence context
chrom.box <- function(AFR, EUR, EAS, SAS, gw, mut){
  #get data
  chrom.dat <- chrom.process.data(AFR, EUR, EAS, SAS, gw, mut)
  
  #find outliers; THIS DOESN'T WORK :(
  AFR.rates <- subset(chrom.dat, chrom.dat$pop == "African")
  AFR.outliers <- subset(AFR.rates, is_outlier(AFR.rates$rate)==TRUE)
  EUR.rates <- subset(chrom.dat, chrom.dat$pop == "European")
  EUR.outliers <- subset(EUR.rates, is_outlier(EUR.rates$rate)==TRUE)
  EAS.rates <- subset(chrom.dat, chrom.dat$pop == "East\nAsian")
  EAS.outliers <- subset(EAS.rates, is_outlier(EAS.rates$rate)==TRUE)
  SAS.rates <- subset(chrom.dat, chrom.dat$pop == "South\nAsian")
  SAS.outliers <- subset(SAS.rates, is_outlier(SAS.rates$rate)==TRUE)
  
  #plot
  c_plot <- ggplot(chrom.dat, aes(pop, rate))+
    geom_boxplot(outlier.color = NA, fill = c("palegreen1", "lightcoral", "steelblue1", 'plum1'))+
    labs(x = "\nPopulation", y = paste("Mutation rate of", mut,"by chromosome\n"))+
    scale_color_manual("", values = c("forest green","red", "dark blue", 'magenta'))+
    
    #add outlier labels
    geom_text(data = EUR.outliers, aes(pop, rate, label = chrom), color = "dark blue", nudge_x = 0.1)+
    geom_text(data = AFR.outliers, aes(pop, rate, label = chrom), color = "forest green", nudge_x = 0.1)+
    geom_text(data = EAS.outliers, aes(pop, rate, label = chrom), color = "red", nudge_x = 0.1)+
    geom_text(data = SAS.outliers, aes(pop, rate, label = chrom), color = "magenta", nudge_x = 0.1)+
    
    #add points
    geom_point(aes(color = factor(chrom.dat$pop)), size = 2, position = position_jitter(width = 0.1))+
    theme(axis.text.x = element_text(size = rel(1.4)), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)),
          legend.position = 'none')#legend position
  return(c_plot)
}

#helper function for ploting which formats the data for ggplot
chrom.process.data <- function(AFR, EUR, EAS, SAS, gw, mut){
  #get indicies for mutation and contextAFR, EUR, EAS, SAS, gw, mut
  i <- which(EUR$Context == mut)
  cntxt <- substr(mut, 1, nchar(mut)-3)
  i.cntxt <- which(gw$Context == cntxt)
  col.e <- ncol(AFR)
  col.s <- col.e-22
  
  #trim summary columns from each input dataframe
  AFR <- AFR[c(col.s: col.e)]
  EUR <- EUR[c(col.s: col.e)]
  EAS <- EAS[c(col.s: col.e)]
  SAS <- SAS[c(col.s: col.e)]
  gw <- gw[-c(1,2)]
  gw_totals <- colSums(gw)
  
  #make output dataframe for plot
  dat <- data.frame(matrix(nrow = 23, ncol = 4))
  colnames(dat) <- c("African", "East\nAsian", "European", "South\nAsian")  
  
  #get rates for each pop
  dat$'European' <- t(EUR[i,]/gw[i.cntxt,]*1.2E-8*gw_totals/colSums(EUR))
  dat$'African' <- t(AFR[i,]/gw[i.cntxt,]*1.2E-8*gw_totals/colSums(AFR))
  dat$'East\nAsian' <- t(EAS[i,]/gw[i.cntxt,]*1.2E-8*gw_totals/colSums(EAS))
  dat$'South\nAsian' <- t(SAS[i,]/gw[i.cntxt,]*1.2E-8*gw_totals/colSums(SAS))
  
  #melt data to necessary format
  row.names(dat) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X')
  dat.m <- melt(t(dat))
  colnames(dat.m) <- c("pop", "chrom", "rate")
  
  return(dat.m)
}

#helper function that returns which elements of a vector x are outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}