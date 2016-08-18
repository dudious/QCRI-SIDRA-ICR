#################################################################
###
### This Script generates a pancancer correlation matrix between
### the DBGS3 immne gene signature and other selected pathways
###
#################################################################

# Setup environment
rm(list=ls())
setwd("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
#setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")
#Dependencies
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (gplots)

#load the data
file.list  <- list.files("./3 ANALISYS/Pathway correlation/",full.names = TRUE)
#extract cancer sets
cancer.sets <- gsub("./3 ANALISYS/Pathway correlation/Immune.correlation.","",file.list)
cancer.sets <- gsub(".DBGS3.RData","",cancer.sets)
#correlation matrix
load (file.list[1])
pathways.correlation.matrix <- data.frame(matrix(ncol=nrow(immune.data), nrow=length(cancer.sets)))
rownames(pathways.correlation.matrix) = cancer.sets
colnames(pathways.correlation.matrix) = rownames(immune.data)
#p.value matrix
p.value.matrix <- pathways.correlation.matrix

#loading loop
for (i in 1:length(cancer.sets)) {
  load (file.list[i])
  pathways.correlation.matrix[i,] <- immune.data$correlation
  p.value.matrix[i,] <- immune.data$p.value
}
pathways.correlation.matrix <- as.matrix(pathways.correlation.matrix)
p.value.matrix <- as.matrix(p.value.matrix)

my.palette <- colorRampPalette(c("blue", "white", "darkgreen"))(n = 297)
dev.new()
heatmap.2(pathways.correlation.matrix,
          col=my.palette,
          density.info="none",
          trace="none",
          cexRow=1.3,cexCol=1.3,margins=c(30,7)
          )

