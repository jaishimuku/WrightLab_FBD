


SFBD <- function(n_extant, n_trees, df,  mu, lambda, psi, row, mers){
  
  num <- max(c(length(mu), length(lambda), length(psi)))
  
  if(length(mu) == 1){
    m1 <- rep(df[row, mu], times = num)} else{
      m1 <- c()
      for(i in mu){m1 <- append(m1, df[row, i])}}
  
  if(length(lambda) == 1){
    l1 <- rep(df[row, lambda], times = num)} else{
      l1 <- c()
      for(i in lambda){l1 <- append(l1, df[row, i])}}
  
  if(length(psi) == 1){
    p1 <- rep(df[row, psi], times = num)} else{
      p1 <- c()
      for(i in psi){p1 <- append(p1, df[row, i])}}
  
  modelledfoss <- sim.fbd.rateshift.taxa(n_extant, n_trees, l1, m1, p1, times = mers, complete = TRUE)
  return(modelledfoss)
}
