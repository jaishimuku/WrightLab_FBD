
mynewfoss <- read.table("C:\\Users\\shrey\\Desktop\\WrightLab\\shreya_foss_3.log", head = TRUE)
library(FossilSim)
library(tidyverse)


#Simulate a 100 trees from one row

num_sims <- 100
sim_output <- rep(NA, times = num_sims) #a vector of num_sims NA values

for(i in 1:num_sims){
  resultFBD <- modifiedSFBD(666, 1, mynewfoss, lambda = "lambda", mu = "mu", psi = c("psi.1.", "psi.2."), 1426, c(0, 61))
  
  origin_time <-  tree.max(resultFBD[[2]][[1]])
  sim_output[i] <- origin_time
}

df_sim_origintime <- data.frame(sim_output) #Put in data frame 


#Incorporating fossil count (to see how the for loop works, see above chunk). Basically removed the sim_output[i]

num_sims <- 3
origin_time <- rep(NA, times = num_sims) #This creates a vector of num_sims NA values
total_foss <- rep(NA, times = num_sims) 

for(i in 1:num_sims){
  resultFBD <- modifiedSFBD(666, 1, mynewfoss, lambda = "lambda", mu = "mu", psi = c("psi.1.", "psi.2."), 1426, c(0, 61))
  
  origin_time[i] <-  tree.max(resultFBD[[2]][[1]])
  total_foss[i] <- count.fossils(resultFBD[[1]])
  
}

my_df <- data.frame(origin_time, total_foss) #Put in data frame (to work for ggplot)


#Histogram

p <- ggplot(data = df_sim_origintime, aes(x = sim_output)) + 
  geom_histogram(binwidth = 4, color = "darkblue", fill = "lightblue") + 
  geom_vline(xintercept = 131, color = "red") +
  scale_x_reverse() + #seq in the -x axis
  labs(title = "Origin Time Histogram", x = "Million of Years Ago", y= "Number of Samples") + #axis labels
  theme(plot.title = element_text(hjust = 0.5)) #centers title





