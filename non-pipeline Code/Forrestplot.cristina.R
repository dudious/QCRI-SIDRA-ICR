# Setup environment
rm(list=ls())
required.packages <- c("forestplot")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(forestplot)

  HR.table <- read.csv (file="./3 ANALISYS/stats.KLRK1.Tertile.csv")
  HR.table <- HR.table[-which(is.na(HR.table$p.value)),]
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <- 
  structure(list(
    mean  = c(NA,HR.table$HR, NA), 
    lower = c(NA,HR.table$lower, NA),
    upper = c(NA,HR.table$upper, NA)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -22L), 
    class = "data.frame")

tabletext<-cbind(
  c("Variable",    as.character(HR.table$Variable) , NA),
  #c("x", rep(100,20), NA, NA),
  #c("y", rep(100,20), NA, NA))
  c("p-value", HR.table$p.value, NA),
  c("HR",      HR.table$HR, NA))

dev.new()
forestplot(tabletext, 
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,TRUE,rep(FALSE,21),TRUE,rep(FALSE,21),TRUE,FALSE),
           #clip=c(0.1,2.5), 
           xlog=TRUE, 
           boxsize = .25,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
