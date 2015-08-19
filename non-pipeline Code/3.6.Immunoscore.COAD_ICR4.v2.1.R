#################################################################
###
### This script calculated the Immunoscore 
###
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
 #Dependencies
 required.packages <- c("plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
#estimate
library(utils)
required.packages <- c("estimate")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages("./1 CODE/R tools/estimate_1.0.11.zip", repos=NULL)
library("plyr")
library("estimate")
source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/read.gct.R")

dir.create("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",showWarnings=FALSE)

# Set Parameters
Cancerset <- "COAD.ICR4"
Geneset <- "ISGS2.IMSUG"       

## Load Data
RNASeq.data <- read.csv ("./3 ANALISYS/TCGA.COAD-ICR4_DBGS3_RNASeq.Data.csv")
inhib.genes <- c("PDCD1","CTLA4","CD274","FOXP3","IDO1","X")
RNASeq.subset <- as.matrix(RNASeq.data[,-which(colnames(RNASeq.data) %in% inhib.genes)])
ISGS2 <- colnames(RNASeq.subset)
IMSUG <- c("PDCD1","CTLA4","CD274","FOXP3","IDO1")
RNASeq.data <- read.csv ("TCGA.RNASeq.wCluster.csv") # Select subset here !!!!!
rownames(RNASeq.data) <- RNASeq.data$Row.names

# Calculate Immunoscore IMSUG
RNASeq.subset <- as.matrix(RNASeq.data[RNASeq.data$Group=="ICR4",which(colnames(RNASeq.data) %in% ISGS2)])
RNASeq.subset.scaled <- scale (RNASeq.subset,scale=FALSE)
RNASeq.subset.scaled <- cbind(RNASeq.subset.scaled,rowMeans(RNASeq.subset.scaled[,-ncol(RNASeq.subset.scaled)]))
colnames(RNASeq.subset.scaled)[ncol(RNASeq.subset.scaled)] <- c("scaled.ISGS2.IS")
immunoscore <- RNASeq.subset.scaled[,c("scaled.ISGS2.IS"),drop=FALSE]
immunoscore <- cbind(immunoscore,rowMeans(RNASeq.subset[,-ncol(RNASeq.subset)]))
colnames(immunoscore) <- c("scaled.ISGS2.IS","unscaled.ISGS2.IS")

# Calculate Immunoscore ISGS2
RNASeq.subset <- as.matrix(RNASeq.data[RNASeq.data$Group=="ICR4",which(colnames(RNASeq.data) %in% IMSUG)])
RNASeq.subset.scaled <- scale (RNASeq.subset,scale=FALSE)
RNASeq.subset.scaled <- cbind(RNASeq.subset.scaled,rowMeans(RNASeq.subset.scaled[,-ncol(RNASeq.subset.scaled)]))
colnames(RNASeq.subset.scaled)[ncol(RNASeq.subset.scaled)] <- c("scaled.IMSUG.IS")
immunoscore <- cbind(immunoscore,RNASeq.subset.scaled[,c("scaled.IMSUG.IS"),drop=FALSE])
immunoscore <- cbind(immunoscore,rowMeans(RNASeq.subset[,-ncol(RNASeq.subset)]))
colnames(immunoscore) <- c("scaled.ISGS2.IS","unscaled.ISGS2.IS","scaled.IMSUG.IS","unscaled.IMSUG.IS")

#save
write.csv (immunoscore,file=(paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.",Cancerset,".",Geneset,".csv")))


ICR.Master.data <- read.csv ("./3 ANALISYS/MASTER FILES/TCGA.COAD-merged.RNASeq_subset_DBGS3.FLTR.Master_ICR4.csv")
ICR.Master.data <- ICR.Master.data[-29,]
rownames(ICR.Master.data) <- ICR.Master.data$X
ICR.Master.data.2 <- merge (ICR.Master.data,immunoscore,by="row.names",all.x=TRUE,all.y=FALSE)
 