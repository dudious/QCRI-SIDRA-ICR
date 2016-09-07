#################################################################
###
### This Script generates a pancancer correlation matrix between
### the DBGS3 immne gene signature and other selected pathways
###
#################################################################

# Setup environment
rm(list=ls())
#setwd("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
#setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")
setwd("/mnt3/wouter/BREAST-QATAR/")
#Dependencies
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (gplots)

# Set Parameters
DL.Method    = "BIOLINKS"     #Choose "ASSEMBLER" or "BIOLINKS"
sample.types = "Selected"     #Alternatives TP , TP_TM , Selected

#load the data
file.list  <- list.files(paste0("./3 ANALISYS/Pathway correlation/",DL.Method),full.names = TRUE)
#extract cancer sets
cancer.sets <- gsub(paste0("./3 ANALISYS/Pathway correlation/",DL.Method,"/",DL.Method,".",sample.types,"Immune.correlation."),"",file.list)
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
#dev.new()
png(paste0("./4 FIGURES/Heatmaps/selected_pathways/Pancancer.Selected.Pathways.",DL.Method,".",sample.types,".png"),res=600,height=15,width=10,unit="in")     # set filename
heatmap.2(pathways.correlation.matrix,
          col=my.palette,
          density.info="none",
          trace="none",
          cexRow=0.8,cexCol=0.8,margins=c(30,7)
          )
dev.off()
