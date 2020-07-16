
#' @title sim_multiple_rows: Simulate tree from multiple rows and apply summary statistics
#' @description Function runs multiple simulations. Can be specified by row.
#' Also provides summary statistics like total fossil count, fossil
#' count per interval, and origin time.
#'
#' @param data data frame
#' @param extant Number of extant taxa
#' @param mus Extinction rate
#' @param lambdas Speciation rate
#' @param psis Fossil sampling rate. Input as a vector
#' @param intervals Vector of time intervals or time bins
#' @param by Selects for specific rows. For example, by = 100 refers to sampling
#' every 100th row. Defaults to FALSE. If FALSE, every row is sampled to produce
#' a simulated phylo object.
#' @param path Specify path to save
#'
#' @return summlist. An object describing useful summary statistics like
#' total fossil count, fossil count per interval, and origin time
#'
#'


sim_multiple_rows <- function(data, extant, mus, lambdas, psis, intervals, by = FALSE, path = c()){

  num <- length(intervals) + 2 #length of time bins plus 2 (the mu and lambda)

  summlist <- list()
  for(i in c(1:num)){summlist[i + 1] <- c()}

  if(by == FALSE){
    for(i in c(1:nrow(data))){

      tempTree <- sim_single_row(n_extant = extant, n_trees = 1, df = data, mu = mus, lambda = lambdas, psi = psis, row = i, mers = intervals)
      tree_stats <- summarize_tree(tempTree, intervals)

      for(i in c(1:length(tree_stats))){
        summlist[[i]] <- append(summlist[[i]], tree_stats[i])
      }
    }
  }

  else{
    for(i in which(c(1:nrow(data) %% by == 0))){

      tempTree <- sim_single_row(n_extant = extant, n_trees = 1, df = data, mu = mus, lambda = lambdas, psi = psis, row = i, mers = intervals)
      tree_stats <- summarize_tree(tempTree, intervals)

      for(i in c(1:length(tree_stats))){
        summlist[[i]] <- append(summlist[[i]], tree_stats[i]) #add to summlist the tree_stats value
      }
    }
  }

  summlist <- as.data.frame(summlist)

  for(i in c(1:ncol(summlist))){
    names(summlist)[1] <- "OriginTime"
    names(summlist)[2] <- "Total_Foss_Count"

    if(i > 2){
      names(summlist)[i] <- paste("Int", i - 2, "Fossils" ,sep = "_")
    }
  }

  if(length(path) != 0){
    write.csv(summlist, file = path)
  }

  return(summlist)
}

