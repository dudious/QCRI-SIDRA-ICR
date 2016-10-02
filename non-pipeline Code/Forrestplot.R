required.packages <- c("forestplot")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(forestplot)

HR.table <- read.csv (file="./6 Project details/Paper/OncoImmunology submission/Revised/Forrest.plot.csv")

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
  c("Gene",    as.character(HR.table$Gene) , NA),
  #c("x", rep(100,20), NA, NA),
  #c("y", rep(100,20), NA, NA))
  c("p-value", HR.table$p.value, NA),
  c("HR",      HR.table$HR, NA))

dev.new()
forestplot(tabletext, 
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,22),TRUE),
           #clip=c(0.1,2.5), 
           xlog=TRUE, 
           boxsize = .25,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
