
#' @title sim_single_row: Simulating trees from a single row of the posterior trace
#'
#' @description Runs modified_sfbdrt with parameter specifications to simulate
#' tree from a single, specified row from the phylogenetic posterior prediction.
#'
#' @param n_extant Number of extant taxa
#' @param n_trees Number of trees to simulate
#' @param df data frame
#' @param mu Extinction rate
#' @param lambda Speciation rate
#' @param psi Fossil sampling rate. Input as a vector
#' @param row Row to sample
#' @param mers Vector of time intervals or time bins
#'
#' @return modelledfoss, the simulated list of phylo objects
#' with parameter specifications
#'
#'



sim_single_row <- function(n_extant, n_trees, df,  mu, lambda, psi, row, mers){

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

  modelledfoss <- modified_sfbdrt(n_extant, n_trees, l1, m1, p1, times = mers, complete = TRUE)
  return(modelledfoss)
}


