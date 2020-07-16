
#' @title summarize_tree: Summarizes simulated tree
#' @description summarize_tree calculates useful statistics like
#' the origin time, total fossils, and fossil count per time interval for
#' a particular simulation.
#'
#' @param phylolist simulated list of phylo objects
#' @param interval.ages Vector of time intervals
#'
#' @return A vector of origin time, total fossils count, and fossil count per interval
#'
#'
#'
summarize_tree <- function(phylolist, interval.ages){
  tree_origin <- tree.max(phylolist[[2]][[1]])
  num_fossils <- count.fossils(phylolist[[1]])
  num_fossils_binned <- count.fossils.binned(phylolist[[1]], interval.ages)
  return(c(tree_origin, num_fossils, num_fossils_binned))
}

