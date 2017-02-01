required.packages <- c("forestplot")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(forestplot)

  HR.table <- read.csv (file="./6 Project details/Paper/OncoImmunology submission/Revised/new images/Forrest.MUTvsWT.plot.csv")

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <- 
  structure(list(
    mean  = c(NA,NA,1.120,1,0.390,0.200,NA,0.710,1,NA,0.290), 
    lower = c(NA,NA,0.390,1,0.150,0.050,NA,0.250,1,NA,0.130),
    upper = c(NA,NA,3.210,1,1.020,0.860,NA,2.050,1,NA,0.670)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")

tabletext<-cbind(
  c("Group",    as.character(HR.table$Group)),
  #c("x", rep(100,20), NA, NA),
  #c("y", rep(100,20), NA, NA))
  c("p-value", rep(NA,6)),
  c("HR",      HR.table$HR))

dev.new()
forestplot(tabletext, 
           cochrane_from_rmeta,new_page = TRUE,
           #is.summary=c(TRUE,  TRUE ,FALSE, FALSE, FALSE, TRUE, FALSE,FALSE,FALSE  TRUE),
           #clip=c(0.1,2.5), 
           xlog=TRUE, 
           boxsize = .25,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))



