SampFBD <- function(n_extant, n_trees, df,  mu, lambda, psi, rows, mers){
  m1 <- df[rows, mu]
  l1 <- df[rows, lambda]
  p1 <- df[rows, psi]
  
  modelledfoss <- sim.fbd.rateshift.taxa(n_extant, n_trees, l1, m1, p1, mers, complete = TRUE)
  return(modelledfoss)
}

