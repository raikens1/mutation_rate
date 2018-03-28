library(ggplot2)
library(plyr)
require(reshape2)
library(gridExtra)
library(grid)
library(lattice)

## CHISQ_PLOTS
##################
# hom.test.plot
#   given: homogeneity test p values and a lower limit for text labeling points
#   do: return an ordered plot of homogeneity test p-value
# sigs.plot
#   given: homogeneity test p values and a lower limit for text labeling points
#   do: return an ordered plot of homogeneity test p-value which shows only the bonferroni-corrected significant results
# volcano.plot
#   given: 2 count dataframes, a dataframe of 2-way X^2 tests, and the lower limit for text labeling poitns
#   do: return a volcano plot a la Harris 2015
##################

#color vector so that plot colorings match Harris
h.colors <- c("dark blue", "magenta", "purple", "forest green", "#0099FF", "red")

#p-value scattter plot for 4-way comparison
hom.test.plot <- function(p.values, lab.lim, NoGGA = F){
  n = length(p.values$Context)
  alpha = 0.05/n
  
  #preprocessing data
  data <- p.values[order(p.values$p),]
  data$Context <- as.character(data$Context)
  data$Context <- factor(data$Context, levels = unique(data$Context))
  
  #vector to decide which points to label
  include <- (data$p < lab.lim)
  if (NoGGA){
    include <- include & (data$X3mer != "GGA->A")
  }
  sigs <- subset(data, include)
  
  #plot
  my.plot <- qplot(data$Context, -log10(data$p)) +
    labs(x = "Polymorphism type", y = expression(-log[10](italic(p))))+
    scale_color_manual("", values = h.colors)+
    geom_point(aes(color=factor(data$X1mer)), size = 3.5)+
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = 2, size = 1)+
    geom_text(data = sigs, aes(sigs$Context, -log10(sigs$p), label = sigs$Context, color = sigs$X1mer), 
              hjust = 0, nudge_x = 1.3, nudge_y = 2.5, size = 5, check_overlap = TRUE)+
    geom_text(aes(15, -log10(alpha)+4, label = "significant (p < 5e-4)", vjust = 0), color = "black", size = 5)+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)),
          legend.text = element_text(size = rel(1.25)), legend.position = c(0.8, 0.8))
  
  return(my.plot)
}

#p-value scattter plot for 4-way comparison which just shows significant results
sigs.plot <- function(p.vals, lab.lim, NoGGA = F){
  n = length(p.vals$Context)
  alpha = 0.05/n
  
  #preprocessing data
  p.values <- subset(p.vals, p.vals$p < alpha)
  data <- p.values[order(p.values$p),]
  data$Context <- as.character(data$Context)
  data$Context <- factor(data$Context, levels = unique(data$Context))
  
  
  #vector to decide which points to label
  include <- (data$p < lab.lim)
  if (NoGGA){
    include <- include & (data$X3mer != "GGA->A")
  }
  sigs <- subset(data, include)
  
  my.plot <- qplot(data$Context, -log10(data$p)) +
    labs(x = "Polymorphism type", y = expression(-log[10](italic(p))))+
    scale_color_manual("", values = h.colors)+
    geom_point(aes(color=factor(data$X1mer)), size = 3)+
    geom_text(data = sigs, aes(sigs$Context, -log10(sigs$p), label = sigs$Context, color = sigs$X1mer), 
              hjust = 0, nudge_x = 1.5, nudge_y = .5, size = 4.5, check_overlap = TRUE)+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)),
          legend.text = element_text(size = rel(1.25)), legend.position = c(0.8, 0.8))
  
  return(my.plot)
}

#volcano plot as in Harris 2015
volcano.plot <- function(counts.1, counts.2, 
                         p.vals, lab.lim, 
                         pops = c("Population 1", "Population 2"), 
                         legpos = c(0.2, 0.75)){
  
  f.1 <- counts.1$Count/sum(counts.1$Count)
  f.2 <- counts.2$Count/sum(counts.2$Count)
  f.diff <- (f.1 -f.2)/f.2
  alpha <- 0.05/(6*length(counts.1$Count))
  
  sigs <- subset(counts.1, p.vals$p < lab.lim)
  f.diff.sigs <- subset(f.diff, p.vals$p < lab.lim)
  p.sigs <- subset(p.vals, p.vals$p < lab.lim)
  
  v_plot <- qplot(f.diff, -log10(p.vals$p)) +
    labs(x = paste("Fold excess in", pops[1]) , 
         y = expression(-log[10](italic(p))), 
         title = paste(pops[1], "v", pops[2])) +
    scale_color_manual("", values = h.colors)+ # set colors
    geom_text(data = sigs, aes(f.diff.sigs, -log10(p.sigs$p), label = sigs$Context, color = sigs$X1mer), 
              vjust = 1, nudge_y =-1, hjust = 1, nudge_x = -.03, size = 3, check_overlap = TRUE)+ # label highly significant
    geom_point(aes(color=factor(counts.1$X1mer)), size = 3)+ # Color pts by one_mer
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = 2)+# significance line
    theme(axis.text.x = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.1)), # adjust text sizes
          axis.title.y = element_text(size = rel(1.1)), axis.text.y = element_text(size = rel(1.35)), 
          legend.text = element_text(size = rel(1.1)), legend.position = legpos,
          title = element_text(size = rel(1.2))) + # legend position and size 
    xlim(c(-.4, .6)) + ylim(c(0, 400))
  
  return(v_plot)
}