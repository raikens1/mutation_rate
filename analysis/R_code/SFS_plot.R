library(ggplot2)

#given: AFR, EUR, EAS, SAS, dataframes of count by AF
#       before, the number of flanking nucleotides preceeding the mutation
#       mut, the mutation of interest
#makes a line plot of proportion by AF
SFS_plot <- function(AFR, EUR, EAS, SAS, mut, before){
  i <- which(AFR$Context == mut)
  
  #based on sequence context, figure out which text columns to remove
  remove <- c(seq(before+2)) #TODO: fix this when you fix your off-by-one error
  
  #calculate mutation proportion across AF bin for each pop
  AFR_prop <- AFR[i, -remove]/colSums(AFR[,-remove])
  EUR_prop <- EUR[i, -remove]/colSums(EUR[,-remove])
  EAS_prop <- EAS[i, -remove]/colSums(EAS[,-remove])
  SAS_prop <- SAS[i, -remove]/colSums(SAS[,-remove])
  
  bins <- as.numeric(colnames(AFR[-remove]))
  
  SFSplot <- plot(bins, EUR_prop, type = "l", col = "dark green",
                  ylim = c(0.01, 0.05), xlab = "log Alternate Allele Frequency", ylab = "proportion", lwd = 3)
  lines(bins, SAS_prop, type = "l", col = "magenta", lwd = 3)
  lines(bins, AFR_prop, type = "l", col = "blue", lwd = 3)
  lines(bins, EAS_prop, type = "l", col = "red", lwd = 3)
  legend("topleft", legend = c("EUR", "SAS", "AFR", "EAS"), lty = 1, col=c("dark green", "magenta", "blue", "red"), lwd = 3)
  
  return(SFSplot)
}

#given: AFR, EUR, EAS, SAS, dataframes of count by AF
#       before, the number of flanking nucleotides preceeding the mutation
#       mut, the mutation of interest
#makes a line plot of proportion by AF
SFS_plot_2 <- function(AFR, AFR.bins, EUR, EUR.bins, EAS, EAS.bins, SAS, SAS.bins, mut, before){
  i <- which(AFR$Context == mut)
  
  #based on sequence context, figure out which text columns to remove
  remove <- c(seq(before+1))
  
  #calculate mutation proportion across AF bin for each pop
  AFR_prop <- AFR[i, -remove]/colSums(AFR[,-remove])
  EUR_prop <- EUR[i, -remove]/colSums(EUR[,-remove])
  EAS_prop <- EAS[i, -remove]/colSums(EAS[,-remove])
  SAS_prop <- SAS[i, -remove]/colSums(SAS[,-remove])
  
  props <- as.numeric(c(AFR_prop, EUR_prop, EAS_prop, SAS_prop))
  bins <- c(AFR.bins$AVG_AF, EUR.bins$AVG_AF, EAS.bins$AVG_AF, SAS.bins$AVG_AF)
  
  labels <- c(rep("AFR", length(AFR_prop)), rep("EUR", length(EUR_prop)), rep("EAS", length(EAS_prop)), rep("SAS", length(SAS_prop)))
  dat <- data.frame(cbind(props, bins, labels))
  names(dat)<- c("proportions", "AF", "POP")
  
  YU = max(props)
  YL = min(props)
  XU = max(log10(bins))
  XL = min(log10(bins))
  
  #ggplot method
  SFSplot <- ggplot(dat, aes(x = bins, y = props, group = labels, color = labels),
                    xlim = c(XL, XU), ylim = c(YL, YU))+
              geom_line(size = 1.5)+
              scale_y_continuous(name = "Proportion of private polymorphism")+
              scale_x_log10(name = "Alternate allele frequency")+
              scale_color_manual(values = c("forestgreen", "red", "darkblue", "magenta"))

  #SFSplot <- plot(log10(EUR.bins$AVG_AF), EUR_prop, type = "l", col = "dark blue", 
  #               ylim = c(YL, YU), xlim = c(XL, XU), xlab = "log Alternate Allele Frequency", ylab = "proportion", lwd = 3)
  #lines(log10(SAS.bins$AVG_AF), SAS_prop, type = "l", col = "magenta", lwd = 3)
  #lines(log10(AFR.bins$AVG_AF), AFR_prop, type = "l", col = "forestgreen", lwd = 3)
  #lines(log10(EAS.bins$AVG_AF), EAS_prop, type = "l", col = "red", lwd = 3)
  #legend("bottomleft", legend = c("EUR", "SAS", "AFR", "EAS"), lty = 1, col=c("dark blue", "magenta", "forestgreen", "red"), lwd = 3)
  
  return(SFSplot)
}