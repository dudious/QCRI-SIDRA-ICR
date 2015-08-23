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
Cancerset <- "COAD-GA"
Geneset <- "DBGS3"       

## Load Data
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)
IMSUG <- c("PDCD1","CTLA4","CD274","FOXP3","IDO1")
RNASeq.subset.IMSUG <- as.matrix(RNASeq.subset[,which(colnames(RNASeq.subset) %in% IMSUG)])
RNASeq.subset.ISGS2 <- as.matrix(RNASeq.subset[,-which(colnames(RNASeq.subset) %in% IMSUG)])

load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
write.table(RNASeq.NORM_Log2,file=paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt"),sep="\t",quote=FALSE)
dir.create(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset),showWarnings=FALSE)

# Calculate estimate score
filterCommonGenes(input.f=paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt"),
                  output.f=paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.input.gct"),
                  id=c("GeneSymbol","EntrezID"))
estimateScore(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.input.gct"),
              paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.score.gct"),
              platform= "illumina")
estimate.gct<-read.table(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.score.gct"),skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe

# Calculate Immunoscore
RNASeq.subset.scaled <- scale (RNASeq.subset,scale=FALSE)
RNASeq.subset.scaled <- cbind(RNASeq.subset.scaled,rowMeans(RNASeq.subset.scaled[,-ncol(RNASeq.subset.scaled)]))
colnames(RNASeq.subset.scaled)[ncol(RNASeq.subset.scaled)] <- c("scaled.IS")
immunoscore <- RNASeq.subset.scaled[,c("scaled.IS"),drop=FALSE]
immunoscore <- cbind(immunoscore,rowMeans(RNASeq.subset[,-ncol(RNASeq.subset)]))
colnames(immunoscore) <- c("scaled.IS","unscaled.IS")
estimate.gct<-t(estimate.gct[,-1])
colnames(estimate.gct) <- estimate.gct[1,]
estimate.gct<-estimate.gct[-1,]
rownames(estimate.gct) <- gsub("\\.","-",rownames(estimate.gct))
rownames(immunoscore) <- gsub("\\.","-",rownames(immunoscore))
mode(estimate.gct) <- "numeric"
estimate.gct <- estimate.gct[rownames(immunoscore),]
immunoscore <- cbind (immunoscore,estimate.gct)
immunoscore <- cbind (immunoscore,rowMeans(RNASeq.subset.IMSUG))
immunoscore <- cbind (immunoscore,rowMeans(RNASeq.subset.ISGS2))
colnames(immunoscore) <- c(colnames(immunoscore[,1:5]),c("unscaled.IS.IMSUG","unscaled.IS.ISGS2"))


write.csv (immunoscore,file=(paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.",Cancerset,".",Geneset,".csv")))


 