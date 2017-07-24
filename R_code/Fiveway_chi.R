Fiveway.Chi <- function(ONE, TWO, THREE, FOUR, FIVE){
  n.contexts = length(ONE$Context)
  
  #make datONEame for results
  result <- data.frame(matrix(ncol=10,nrow=n.contexts))
  colnames(result) <- c("Context", "X5mer","X3mer", "X1mer", "ONE.Count", "TWO.Count", "THREE.Count", "FOUR.Count","FIVE.Count","p")
  result$Context <- ONE$Context
  result$X5mer <- ONE$X5mer
  result$X3mer <- ONE$X3mer
  result$X1mer <- ONE$X1mer
  result$ONE.Count <- ONE$Count
  result$TWO.Count <- TWO$Count
  result$THREE.Count <- THREE$Count
  result$FOUR.Count <- FOUR$Count
  result$FIVE.Count <- FIVE$Count
  
  #start setting up tables
  sum.1 <- sum(ONE$Count)
  sum.2 <- sum(TWO$Count)
  sum.3 <- sum(THREE$Count)
  sum.4 <- sum(FOUR$Count)
  sum.5 <- sum(FIVE$Count)
  
  #set up table and run test for each context
  for (i in 1:n.contexts){
    c.a <- c(ONE$Count[i], TWO$Count[i], THREE$Count[i], FOUR$Count[i], FIVE$Count[i])
    c.b <- c(sum.1, sum.2, sum.3, sum.4, sum.5) - c.a
    data <- cbind(c.a, c.b)
    result$p[i] <- chisq.test(data)$p.value #continuity correction?
  }
  
  return(result)
}
