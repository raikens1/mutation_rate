#this code is modified from http://www.gettinggeneticsdone.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
library(ggplot2)

#color vector so that plot colorings match Harris
h.colors <- c("dark blue", "magenta", "purple", "forest green", "#0099FF", "red")

qq.labels <- function(pdata, lab.lim, title="Quantile-quantile plot of p-values", NoGGA = F) {
  # Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
  #preprocess data to order and remove nas 
  pdata <- pdata[complete.cases(pdata$p),]
  pdata <- pdata[order(pdata$p),]
  pdata$Context <- as.character(pdata$Context)
  pdata$Context <- factor(pdata$Context, levels = unique(pdata$Context))
  
  #vector to record which SNPs to label
  include <- (pdata$p < lab.lim)
  if (NoGGA){
    include <- include & (pdata$X3mer != "GGA->A")
  }
  sigs <- subset(pdata, include)
  
  #x and y coordinates in qq plot
  o <- -log10(sort(pdata$p,decreasing=F))
  e <- -log10( 1:length(o)/length(o) )
  
  
  #plot
  plot <- qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o)))+
    scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
    scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))+
    scale_color_manual("", values = h.colors)+
    geom_point(aes(color=factor(pdata$X1mer)), size = 1.75)+
    geom_text(data = sigs, 
              aes(subset(e, include), subset(o, include), label = sigs$Context, color = sigs$X1mer), 
              hjust = 1, nudge_x = -0.03, nudge_y = .5, size = 4, check_overlap = T)+
    labs(title = title)+
    theme(axis.text.x = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)), #adjust text sizes
          axis.title.y = element_text(size = rel(1.3)), axis.text.y = element_text(size = rel(1.35)), 
          legend.text = element_text(size = rel(1.2)),
          title = element_text(size = rel(1.2))) +
    geom_abline(intercept=0,slope=1, col="red")
  
  return(plot)
}

qq.generic <- function(pvector, title="Quantile-quantile plot of p-values", spartan=F) {
  # Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
  o <- -log10(sort(pvector,decreasing=F))
  e <- -log10( 1:length(o)/length(o) )
  # you could use base graphics
  
  #You'll need ggplot2 installed to do the rest
  plot <- qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o)))+
    scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
    scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))+
    labs(title = title)+
    theme(axis.text.x = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)), #adjust text sizes
          axis.title.y = element_text(size = rel(1.3)), axis.text.y = element_text(size = rel(1.35)), 
          legend.text = element_text(size = rel(1.2)),
          title = element_text(size = rel(1.2))) +
    geom_abline(intercept=0,slope=1, col="red")
  
  return(plot)
}