library(ggplot2)
library(plyr)
require(reshape2)
library(gridExtra)
library(grid)
library(lattice)

#Returns a histogram of all the Xmers for a given threemer subcontext for one population
subcontext.3mer.hist <- function(priv, mut, wdth){
  subs <- subset(priv, priv$X3mer == mut)
  result <- hist(subs$Rate, 
                 col = "gray", border = "white")
  return(result)
}

#Returns a ggplot histogram of all the Xmers for a given threement subcontext for ALL populations.
#wdth <- bin width
#xint and yint <- x and y axis limits, respectively
subcontext.3mer.gghist <- function(AFR, EUR, EAS, SAS, mut, wdth, xint, yint){
  AFR.mut <- subset(AFR, AFR$X3mer == mut)
  EUR.mut <- subset(EUR, EUR$X3mer == mut)
  EAS.mut <- subset(EAS, EAS$X3mer == mut)
  SAS.mut <- subset(SAS, SAS$X3mer == mut)
  
  #AFR
  AFR.hist <- qplot(AFR.mut$Rate, geom="histogram", binwidth = wdth, xlim = xint, ylim = yint,
                    main = mut, xlab = "AFR mutation rate", col = I("dark grey")) +
    annotate("text", x=2/3*xint[2], y=2/3*yint[2],
             label = paste("median =", signif(median(AFR.mut$Rate), 3),
                           "\n mean =", signif(mean(AFR.mut$Rate), 3), 
                           "\n stdev =", signif(sd(AFR.mut$Rate), 3))) + geom_vline(xintercept = median(AFR.mut$Rate), colour = "red") + geom_vline(xintercept = mean(AFR.mut$Rate), colour = "blue") 
  
  #EUR 
  EUR.hist <- qplot(EUR.mut$Rate, geom="histogram", binwidth = wdth, xlim = xint, ylim = yint,
                    xlab = "EUR mutation rate", col = I("dark grey")) +
    annotate("text", x=2/3*xint[2], y=2/3*yint[2], 
             label = paste("median =", signif(median(EUR.mut$Rate), 3),
                           "\n mean =", signif(mean(EUR.mut$Rate), 3), 
                           "\n stdev =", signif(sd(EUR.mut$Rate), 3))) + geom_vline(xintercept = median(EUR.mut$Rate), colour = "red") + geom_vline(xintercept = mean(EUR.mut$Rate), colour = "blue") 
  
  #EAS
  EAS.hist <- qplot(EAS.mut$Rate, geom="histogram", binwidth = wdth, xlim = xint, ylim = yint,
                    xlab = "EAS mutation rate", col = I("dark grey")) +
    annotate("text", x=2/3*xint[2], y=2/3*yint[2], 
             label = paste("median =", signif(median(EAS.mut$Rate), 3),
                           "\n mean =", signif(mean(EAS.mut$Rate), 3), 
                           "\n stdev =", signif(sd(EAS.mut$Rate), 3))) + geom_vline(xintercept = median(EAS.mut$Rate), colour = "red") + geom_vline(xintercept = mean(EAS.mut$Rate), colour = "blue") 
  
  #SAS
  SAS.hist <- qplot(SAS.mut$Rate, geom="histogram", binwidth = wdth, xlim = xint, ylim = yint, 
                    xlab = "SAS mutation rate", col = I("dark grey")) +
    annotate("text", x=2/3*xint[2], y=2/3*yint[2], 
             label = paste("median =", signif(median(SAS.mut$Rate), 3),
                           "\n mean =", signif(mean(SAS.mut$Rate), 3), 
                           "\n stdev =", signif(sd(SAS.mut$Rate), 3))) + geom_vline(xintercept = median(SAS.mut$Rate), colour = "red") + geom_vline(xintercept = mean(SAS.mut$Rate), colour = "blue") 
  return(grid.arrange(AFR.hist, EUR.hist, EAS.hist, SAS.hist, ncol=1))
}