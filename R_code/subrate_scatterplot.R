library(ggplot2)
# given: count1 and count2: pop-specific polymorphism dataframes; Mut: 3mer subcontext of interest,
# generates a scatter plot of rate for mut in population 1 by population 2
subrate.scplot <- function(count1, count2, mut){
  pop1.mut <- subset(count1, count1$X3mer == mut)
  pop2.mut <- subset(count2, count2$X3mer == mut)
  #dat <- data.frame(cbind(pop1.mut$Context, pop1.mut$Rate, pop2.mut$Rate))
  sample.size <- pop1.mut$Count + pop2.mut$Count
  #lab <- subset(dat, pop1.mut$Rate/pop2.mut$Rate > 4 | pop1.mut$Rate/pop2.mut$Rate < 0.25)
  
  p_plot <- qplot(pop1.mut$Rate, pop2.mut$Rate, xlim = c(0,1.4e-8), ylim = c(0,1.4e-8)) +
    labs(x = "AAC->C rate in KHV", y = "AAC->C rate in JPT") +
    theme(legend.title = element_blank(), 
          legend.justification= c(1,0), 
          legend.position=c(1,0),
          axis.text=element_text(size=10),
          text = element_text(size = 17))+
    geom_point(aes(colour = log(base = 10, x = sample.size)), size = 3) + 
    scale_colour_gradient(low = "orange", high = "blue")+
    geom_abline(aes(intercept=0,slope=1,show.legend = FALSE), size = 1) + 
    geom_text(aes(pop1.mut$Rate, pop2.mut$Rate), label = pop1.mut$Context, check_overlap = T,
              hjust = 0, nudge_x = 5e-10)+
    coord_fixed()
  
  return(p_plot)
}

# given: count files for each pop and a given 3mer type, 
# plots all 5mers with that 3mer subcontext on a line plot
line.subcontext <- function(AFR, EUR, EAS, SAS, mut){
  i <- which(AFR$X3mer == mut)
  
  AFR.dat <- cbind(AFR$Rate[i],AFR$Context[i], rep("AFR", 16))
  colnames(AFR.dat)<- c("Rate", "Context", "Pop")
  EUR.dat <- cbind(EUR$Rate[i],EUR$Context[i], rep("EUR", 16))
  colnames(EUR.dat)<- c("Rate", "Context", "Pop")
  EAS.dat <- cbind(EAS$Rate[i],EAS$Context[i], rep("EAS", 16))
  colnames(EAS.dat)<- c("Rate", "Context", "Pop")
  SAS.dat <- cbind(SAS$Rate[i],SAS$Context[i], rep("SAS", 16))
  colnames(SAS.dat)<- c("Rate", "Context", "Pop")
  
  data <- data.frame(rbind(AFR.dat, EUR.dat, EAS.dat, SAS.dat))
  
  myplot <- ggplot(data, aes(Pop, as.numeric(as.character(Rate)), group = Context, color = Context)) +
    geom_point(size = 4) +
    geom_line(size = 1.2) +
    labs(y = paste("Estimated mutation rate of ", mut,"\n")) + #y axis label
    theme(axis.text.x = element_text(size = rel(1.4)), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)))
  
  return(myplot)
}

line.subcontext.3mer <- function(AFR, EUR, EAS, SAS, mut){
  i <- which(AFR$X1mer == mut)
  
  AFR.dat <- cbind(AFR$Rate[i],AFR$Context[i], rep("AFR", 16))
  colnames(AFR.dat)<- c("Rate", "Context", "Pop")
  EUR.dat <- cbind(EUR$Rate[i],EUR$Context[i], rep("EUR", 16))
  colnames(EUR.dat)<- c("Rate", "Context", "Pop")
  EAS.dat <- cbind(EAS$Rate[i],EAS$Context[i], rep("EAS", 16))
  colnames(EAS.dat)<- c("Rate", "Context", "Pop")
  SAS.dat <- cbind(SAS$Rate[i],SAS$Context[i], rep("SAS", 16))
  colnames(SAS.dat)<- c("Rate", "Context", "Pop")
  
  data <- data.frame(rbind(AFR.dat, EUR.dat, EAS.dat, SAS.dat))
  
  myplot <- ggplot(data, aes(Pop, as.numeric(as.character(Rate)), group = Context, color = Context)) +
    geom_point(size = 4) +
    geom_line(size = 1.2)
  
  return(myplot)
}

line.subcontext.7mer <- function(AFR, EUR, EAS, SAS, mut){
  i <- which(AFR$X5mer == mut)
  
  AFR.dat <- cbind(AFR$Rate[i],AFR$Context[i], rep("AFR", 16))
  colnames(AFR.dat)<- c("Rate", "Context", "Pop")
  EUR.dat <- cbind(EUR$Rate[i],EUR$Context[i], rep("EUR", 16))
  colnames(EUR.dat)<- c("Rate", "Context", "Pop")
  EAS.dat <- cbind(EAS$Rate[i],EAS$Context[i], rep("EAS", 16))
  colnames(EAS.dat)<- c("Rate", "Context", "Pop")
  SAS.dat <- cbind(SAS$Rate[i],SAS$Context[i], rep("SAS", 16))
  colnames(SAS.dat)<- c("Rate", "Context", "Pop")
  
  data <- data.frame(rbind(AFR.dat, EUR.dat, EAS.dat, SAS.dat))
  
  myplot <- ggplot(data, aes(Pop, as.numeric(as.character(Rate)), group = Context, color = Context)) +
    geom_point(size = 4) +
    geom_line(size = 1.2) +
    labs(y = paste("Estimated mutation rate of ", mut,"\n")) + #y axis label
    theme(axis.text.x = element_text(size = rel(1.4)), axis.title.x = element_blank(), #adjust text sizes
          axis.title.y = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.4)))
  
  return(myplot)
}