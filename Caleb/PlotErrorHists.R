library(tidyverse)
library(lubridate)
library(FossilSim)
library(gridExtra)

table <- read.csv("./Wright Lab/sampedtreestats.csv")
mydata <- read.table("./Wright Lab/caleb_foss.log", header = TRUE)

SummedLog <- GetFossilStats(mydata, Intervals = c(0, 60, 80), by = 100)
newFrame <- data.frame(table, SummedLog)

newFrame$OriginError <- abs(((newFrame$OriginTime - newFrame$Origin_Time)/newFrame$Origin_Time)*100)
newFrame$TotalFossilError <- abs(((newFrame$Total_Foss_Count - newFrame$Total_Fossil_Count)/newFrame$Total_Fossil_Count)*100)
newFrame$IntOneError <- abs(((newFrame$Int_1_Fossils - newFrame$Interval_1_Fossils)/newFrame$Interval_1_Fossils)*100)
newFrame$IntTwoError <- abs(((newFrame$Int_2_Fossils - newFrame$Interval_2_Fossils)/newFrame$Interval_2_Fossils)*100)
newFrame$IntThreeError <- abs(((newFrame$Int_3_Fossils - newFrame$Interval_3_Fossils)/newFrame$Interval_3_Fossils)*100)

OriginErrorPlot <- ggplot(newFrame, aes(X, OriginError)) + geom_col()
FossilErrorPlot <- ggplot(newFrame, aes(X, TotalFossilError)) + geom_col()
Int1ErrorPlot <- ggplot(newFrame, aes(X, IntOneError)) + geom_col()
Int2ErrorPlot <- ggplot(newFrame, aes(X, IntTwoError)) + geom_col()
Int3ErrorPlot <- ggplot(newFrame, aes(X, IntThreeError)) + geom_col()

grid.arrange(OriginErrorPlot, FossilErrorPlot, Int1ErrorPlot, Int2ErrorPlot, Int3ErrorPlot, ncol = 3, nrow = 2)
