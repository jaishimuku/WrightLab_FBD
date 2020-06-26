

Summ_Samples <- function(data, extant, mus, lambdas, psis, intervals, by = FALSE, path = c()){

num <- length(intervals) + 2

summlist <- list()
for(i in c(1:num)){summlist[i + 1] <- c()}
  


if(by == FALSE){
  for(i in c(1:nrow(data))){
    
    tempTree <- SFB(n_extant = extant, n_trees = 1, df = data, mu = mus, lambda = lambdas, psi = psis, row = i, mers = intervals)
    tree_stats <- SummPhylo(tempTree, intervals)
    
    for(i in c(1:length(tree_stats))){
      summlist[[i]] <- append(summlist[[i]], tree_stats[i])
    }
  }
}

else{
  for(i in which(c(1:nrow(data) %% by == 0))){
  
    tempTree <- SFB(n_extant = extant, n_trees = 1, df = data, mu = mus, lambda = lambdas, psi = psis, row = i, mers = intervals)
    tree_stats <- SummPhylo(tempTree, intervals)
  
    for(i in c(1:length(tree_stats))){
      summlist[[i]] <- append(summlist[[i]], tree_stats[i])
    }
  }
}



summlist <- as.data.frame(summlist)

for(i in c(1:ncol(summlist))){
  names(summlist)[1] <- "OriginTime"
  names(summlist)[2] <- "Total_Foss_Count"
  
  if(i > 2){
    names(summlist)[i] <- paste("Int", i - 2, "Fossils", sep = "_")
  }
}




if(length(path) != 0){
  write.csv(summlist, file = path)
} 

  
return(summlist)
}






