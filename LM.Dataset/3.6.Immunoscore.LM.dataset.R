#################################################################
###
### This script calculated the Immunoscore 
###
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")
#setwd("~/Dropbox/BREAST_QATAR/")
 #Dependencies
 required.packages <- c("plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
#estimate
library(utils)
required.packages <- c("estimate")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages("~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/R tools/estimate_1.0.11.zip", repos=NULL)
library("plyr")
library("estimate")
#source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/read.gct.R")

dir.create("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",showWarnings=FALSE)
 
## Load Data
load (paste0("./2 DATA/SUBSETS/LM.Dataset/LM.Dataset.MA.subset.DBGS3.RData"))
MA.subset <- as.matrix(MA.subset)
dir.create(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/LM.DATA"),showWarnings=FALSE)

# Calculate estimate score
#filterCommonGenes(input.f=paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt"),
#                  output.f=paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.input.gct"),
#                  id=c("GeneSymbol","EntrezID"))
#estimateScore(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.input.gct"),
#              paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.score.gct"),
#              platform= "illumina")
#estimate.gct<-read.table(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.score.gct"),skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe

# Calculate Immunoscore
MA.subset.scaled <- scale(MA.subset,scale=FALSE)
MA.subset.scaled <- cbind(MA.subset.scaled,rowMeans(MA.subset.scaled[,-ncol(MA.subset.scaled)]))
colnames(MA.subset.scaled)[ncol(MA.subset.scaled)] <- c("scaled.IS")
immunoscore <- MA.subset.scaled[,c("scaled.IS"),drop=FALSE]
immunoscore <- cbind(immunoscore,rowMeans(MA.subset[,-ncol(MA.subset)]))
colnames(immunoscore) <- c("scaled.IS","unscaled.IS")
#estimate.gct<-t(estimate.gct[,-1])
#colnames(estimate.gct) <- estimate.gct[1,]
#estimate.gct<-estimate.gct[-1,]
##rownames(estimate.gct) <- gsub("\\.","-",rownames(estimate.gct))
#rownames(immunoscore) <- gsub("\\.","-",rownames(immunoscore))
#mode(estimate.gct) <- "numeric"
#estimate.gct <- estimate.gct[rownames(immunoscore),]
#immunoscore <- cbind (immunoscore,estimate.gct)
write.csv (immunoscore,file=(paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.LM.Data.csv")))


 