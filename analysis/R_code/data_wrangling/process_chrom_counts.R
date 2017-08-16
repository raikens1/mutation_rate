
# makes a better formatted table of info from a chrom_counts file
# gw_counts should have the counts of # of genomewide appearances for the context of interest
# subcontext_ref should be a dataframe from an nmer_mutations_ref file
# NOTE: even if X is not to be included, the gw_counts df passed in may include X.  (Should work either way... ?)

process_chrom_counts <- function(chrom_counts, gw_counts, subcontext_ref, X = T){
  with_totals <- addtotals(chrom_counts, X)
  with_gw_counts <- add_gw_counts(with_totals, gw_counts, X)
  with_gw_counts$Rate <- head(get_rate(with_gw_counts$Count, with_gw_counts$context_in_genome), -1)
  with_gw_counts$One_mer <- NULL
  
  # reorder to push summary columns to the left
  if (X) {
    with_gw_counts <- with_gw_counts[,c("Context", "Count", "Rate", "context_in_genome","chr1", "chr2", 'chr3', 'chr4', "chr5", "chr6", 'chr7', 'chr8', "chr9", "chr10", 'chr11', 'chr12', "chr13", "chr14", 'chr15', 'chr16', "chr17", "chr18", 'chr19', 'chr20', "chr21", "chr22", 'chrX')]
  }
  else {
    with_gw_counts <- with_gw_counts[,c("Context", "Count", "Rate", "context_in_genome","chr1", "chr2", 'chr3', 'chr4', "chr5", "chr6", 'chr7', 'chr8', "chr9", "chr10", 'chr11', 'chr12', "chr13", "chr14", 'chr15', 'chr16', "chr17", "chr18", 'chr19', 'chr20', "chr21", "chr22")]
  }
 return(cbind(subcontext_ref, with_gw_counts[,-1]))
}

# renames columns and adds a "Count" column with the genome wide count of polymorphism
addtotals <- function(df, X){
  df$Counts = rep(0, length(df$Context))
  if (X){
    for (i in 1:length(df$Context)){
      df$Counts[i] <- sum(df[i,][-c(1,25)])
    }
    colnames(df)<- c("Context", "chr1", "chr2", 'chr3', 'chr4', "chr5", "chr6", 'chr7', 'chr8', "chr9", "chr10", 'chr11', 'chr12', "chr13", "chr14", 'chr15', 'chr16', "chr17", "chr18", 'chr19', 'chr20', "chr21", "chr22", 'chrX', 'One_mer', 'Count')
  }
  else {
    df <- df[,-c(24)]
    for (i in 1:length(df$Context)){
      df$Counts[i] <- sum(df[i,][-c(1,24)])
    }
    colnames(df)<- c("Context", "chr1", "chr2", 'chr3', 'chr4', "chr5", "chr6", 'chr7', 'chr8', "chr9", "chr10", 'chr11', 'chr12', "chr13", "chr14", 'chr15', 'chr16', "chr17", "chr18", 'chr19', 'chr20', "chr21", "chr22", 'One_mer', 'Count')
  }
  return(df)
}

# add column with number of appearances in gw_counts
add_gw_counts <- function(priv, gw, X){
  if (!X) {
    gw$GW_total <- gw$GW_total - gw$X
  }
  priv$context_in_genome <- rep(0, length(priv$Context)) 
  for (i in 1:length(gw$GW_total)){
    i.1 <- 1+3*(i-1)
    priv$context_in_genome[i.1] <- gw$GW_total[i]
    priv$context_in_genome[i.1+1] <- gw$GW_total[i]
    priv$context_in_genome[i.1+2] <- gw$GW_total[i]
  }
  return(priv)
}

# outputs a vector of mutation rate estimates, with the scaling factor, l, appended
# counts <- vector of private polymorphism counts for a pop
# gw <- vector. for each mutation in counts, gives # of appearances of contexts in genome
get_rate <- function(counts, gw){
  # find scaling factor, l
  l <- 1.2E-8/(sum(counts)/sum(gw/3))
  return(c(l*counts/gw, l))
}