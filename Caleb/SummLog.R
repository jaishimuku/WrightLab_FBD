
GetFossilStats <- function(data, Intervals, by = FALSE){

#Get length of dataframe need
num <- length(Intervals) + 2
#Construct precursor list for dataframe
summlist <- list()
for(i in c(1:num)){summlist[i + 1] <- c()}
#Get row in original dataframe corresponding to simulated fossils
fossils <- which(substr(names(data), start = 1, stop = 7) == "fossil.")
#Count total fossils
tot_foss_count <- length(fossils)
#Get rows in original dataframe corresponding to simulated fossil dates
Times <- which(substr(names(data), start = 1, stop = 2) == "t.")

#Get data if using all rows
if(by == FALSE){
  for(i in c(1:nrow(data))){
    #Get origin time at proper row
    Origin <- data$origin_time[i]
    #Get interval horizons
    Horizons <- c(Intervals, Origin)
    for(z in c((1:(length(Intervals))))){
        
      counts <- length(which(data[i, Times] > Horizons[z] & data[i, Times] < Horizons[(z + 1)]))
      
      summlist[[z]] <- append(summlist[[z]], counts)
    }
      summlist[[length(Horizons)]] <- append(summlist[[length(Horizons)]], Origin)
  }
  nsims <- nrow(data)
}
else{
  for(i in which(c(1:nrow(data) %% by == 0))){
    #Get origin time at proper row
    Origin <- data$origin_time[i]
    #Get interval horizons
    Horizons <- c(Intervals, Origin)
    for(z in c((1:(length(Intervals))))){
      
      counts <- length(which(data[i, Times] > Horizons[z] & data[i, Times] < Horizons[(z + 1)]))
      
      summlist[[z]] <- append(summlist[[z]], counts)
    }
      summlist[[length(Horizons)]] <- append(summlist[[length(Horizons)]], Origin)
  }
  nsims <- length(which(c(1:nrow(data) %% by == 0)))
}


summlist[[length(Intervals) + 2]] <- rep(tot_foss_count, 108)
summlist <- as.data.frame(summlist)

for(i in c(1:length(Intervals))){
  names(summlist)[i] <- paste("Interval", i, "Fossils", sep = "_")
}

names(summlist)[ncol(summlist)] <- "Total_Fossil_Count"
names(summlist)[ncol(summlist) - 1] <- "Origin_Time"

return(summlist)
}
