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
 required.packages <- c("ggplot2", "plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)

#estimate
library(utils)
required.packages <- c("estimate")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages("./1 CODE/R tools/estimate_1.0.11.zip", repos=NULL)
library("ggplot2")
library("plyr")
library("estimate")
source ("./1 CODE/R tools/read.gct.R")

#dir.create("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",showWarnings=FALSE)

## Load Data
load ("./2 DATA/SUBSETS/GBM//MA.subset.ISGS.RData")
load ("./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.RData")
dim(AffyData.NORM)#529 X 12042 I need of genes x samples
AffyData.NORM<-t(AffyData.NORM)
dim(AffyData.NORM)#12042x 529
write.table(AffyData.NORM,file="./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.txt",sep="\t",quote=FALSE)
verifica<-read.delim("./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.txt")
rm(verifica)
# Calculate estimate score
filterCommonGenes(input.f="./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.txt",
                  output.f="./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.MA.Affy.GBM.estimate.input.gct",
                  id=c("GeneSymbol","EntrezID"))
estimateScore("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.MA.Affy.GBM.estimate.input.gct",
              "./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.MA.Affy.GBM.estimate.score.gct",
              platform= "affymetrix")
estimate.gct<-read.table("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.MA.Affy.GBM.estimate.score.gct",skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe
#write.table(estimate.gct,file="./3 ANALISYS/IMMUNOSCORE/ESTIMATE/estimate.score.GBM.Affy.txt",quote=F,sep="\t")

# Calculate Immunoscore
MA.subset.scaled <- scale (MA.subset,scale=FALSE)
MA.subset.scaled <- cbind(MA.subset.scaled,rowMeans(MA.subset.scaled[,-ncol(MA.subset.scaled)]))
colnames(MA.subset.scaled)[ncol(MA.subset.scaled)] <- c("scaled.IS")
immunoscore <- MA.subset.scaled[,c("scaled.IS"),drop=FALSE]
immunoscore <- cbind(immunoscore,rowMeans(MA.subset[,-ncol(MA.subset)]))
colnames(immunoscore) <- c("scaled.IS","unscaled.IS")
estimate.gct<-t(estimate.gct[,-1])
colnames(estimate.gct) <- estimate.gct[1,]
estimate.gct<-estimate.gct[-1,]
rownames(estimate.gct) <- gsub("\\.","-",rownames(estimate.gct))
rownames(immunoscore) <- gsub("\\.","-",rownames(immunoscore))
mode(estimate.gct) <- "numeric"
estimate.gct <- estimate.gct[rownames(immunoscore),]
immunoscore <- cbind (immunoscore,estimate.gct)
write.csv (immunoscore,file=("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.MA.Affy.GBM.ISGS.csv"))


 