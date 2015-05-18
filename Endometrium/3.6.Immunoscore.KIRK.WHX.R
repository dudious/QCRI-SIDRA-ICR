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

dir.create("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",showWarnings=FALSE)

## Load Data
load ("./2 DATA/SUBSETS/KIRC/TCGA.KIRC.RNASeq.subset.ISGS.RData")
RNASeq.subset <- as.matrix(RNASeq.subset)

load ("./2 DATA/TCGA RNAseq/RNASeq_KIRC_EDASeq/KIRC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
write.table(RNASeq.NORM_Log2,file="./2 DATA/TCGA RNAseq/RNASeq_KIRC_EDASeq/KIRC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt",sep="\t",quote=FALSE)

# Calculate estimate score
filterCommonGenes(input.f="./2 DATA/TCGA RNAseq/RNASeq_KIRC_EDASeq/KIRC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt",
                  output.f="./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.KIRC.estimate.input.gct",
                  id=c("GeneSymbol","EntrezID"))
estimateScore("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.KIRC.estimate.input.gct",
              "./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.KIRC.estimate.score.gct",
              platform= "illumina")
estimate.gct<-read.table("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.KIRC.estimate.score.gct",skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe

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
write.csv (immunoscore,file=("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.KIRC.ISGS.csv"))


 