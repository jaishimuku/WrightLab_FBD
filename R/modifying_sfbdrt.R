#Modifying the sim.fbd.rateshift.taxa function from FossilSim to output
#Fossil object and SAtree

#' @title n.ages: Accompanying function for sim.fbd.rateshift.taxa from FossilSim package
#'
#' @param tree
#'
#' @return node.ages
#'
#'
#'

n.ages <- function(tree){

  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))

  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)

  return(node.ages)
}


#' @title modified_sfbdrt: Modified the sfbdrt(sim.fbd.rateshift.taxa) function
#'  to return fossil and SAtree objects
#'
#' @description sim.fbd.rateshift.taxa function from the FossilSim Package,
#' but it returns a list of the fossil along with the SAtree object.
#' It also adds an (h) column to the fossil dataframe so that the fossil count
#' functions (summary statistic) will work.
#'
#' @param n Number of extant taxa
#' @param numbsim Number of trees to simulate
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param psi Fossil sampling rate. Input as a vector.
#' @param times Vector of time intervals.
#' @param complete Defaults to FALSE. If TRUE, the tree including the extint
#' lineages and non-sampled lineages is returned.
#'
#' @return list(f,trees). Returns fossil object, f, and simulated SAtrees with n extant
#' sampled tips
#'
#'
#'
modified_sfbdrt <- function (n, numbsim, lambda, mu, psi, times, complete = FALSE)
{
  if (length(lambda) != length(times))
    stop("Length mismatch between rate shift times and birth rates")
  if (length(mu) != length(times))
    stop("Length mismatch between rate shift times and death rates")
  if (length(psi) != length(times))
    stop("Length mismatch between rate shift times and sampling rates")
  trees = TreeSim::sim.rateshift.taxa(n, numbsim, lambda,
                                      mu, rep(1, length(times)), times, complete = TRUE)
  for (i in 1:length(trees)) {
    t = trees[[i]]
    origin = max(n.ages(t)) + t$root.edge
    horizons = c(times, origin)
    f <- sim.fossils.intervals(tree = t, interval.ages = horizons,
                               rates = psi)
    tree = SAtree.from.fossils(t, f)
    node.ages = n.ages(tree)
    if (complete == FALSE) {
      fossil.tips = is.extinct(tree, tol = 1e-06)
      sa.tips = tree$tip.label[tree$edge[, 2][(tree$edge[,
                                                         2] %in% 1:length(tree$tip.label)) & (tree$edge.length ==
                                                                                                0)]]
      unsampled.tips = fossil.tips[!(fossil.tips %in%
                                       sa.tips)]
      tree = ape::drop.tip(tree, unsampled.tips)
      node.ages = n.ages(tree)
    }
    trees[[i]] = tree
    trees[[i]]$root.edge = origin - max(node.ages)
    trees[[i]] = SAtree(trees[[i]], complete)
  }

  f$h <- (f$hmin + f$hmax)/2

  return(list(f, trees))
}

