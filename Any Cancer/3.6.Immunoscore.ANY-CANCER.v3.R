#################################################################
###
### This script calculated the Immunoscore 
###
### ESTIMATE IS NOT WORKING ANYMORE oct/2016
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")
#setwd("/mnt3/wouter/BREAST-QATAR/")
#Dependencies
 required.packages <- c("plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
#estimate
library(utils)
required.packages <- c("estimate")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
#if(length(missing.packages)) install.packages("./1 CODE/R tools/estimate_1.0.11.zip", repos=NULL)
#if(length(missing.packages)) install.packages("/mnt3/wouter/QCRI-SIDRA-ICR/R tools/estimate_1.0.11.zip", repos=NULL)
library("plyr")
#library("estimate")
#source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/read.gct.R")
#source ("/mnt3/wouter/QCRI-SIDRA-ICR/R tools/read.gct.R")

dir.create("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",showWarnings=FALSE)

# Set Parameters
Cancersets        = "ALL"
Geneset           = "DBGS3"   
DL.Method         = "PANCANCER.CLEAN.v2"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Split"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer
Genedatabase      = "Gene_selection_v2.6.txt"

#Load data
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
gene.list <- read.csv (paste0("./2 DATA/SUBSETS/",Genedatabase))                                 # Select subset here !!!!! and change filename below !!!!
gene.list.ALL <- as.character(gene.list[which(gene.list[,"DBGS3"]==1),1])
gene.list.INH <- as.character(gene.list[which(gene.list[,"ImSuGS"]==1),1])
a<-which(gene.list[,"DBGS3"]==1)
b<-which(gene.list[,"ImSuGS"]==1)
gene.list.ACT <- as.character(gene.list[a[-which(a%in%b)],1])

#test
#Cancerset = "BLCA"

# DO ALL
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}

## Load Data
load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)

#load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".NORMALIZED.LOG2.RData"))
#write.table(RNASeq.NORM_Log2,file=paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt"),sep="\t",quote=FALSE)
#dir.create(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset),showWarnings=FALSE)

# Calculate estimate score
#filterCommonGenes(input.f=paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt"),
#                  output.f=paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.input.gct"),
#                  id=c("GeneSymbol","EntrezID"))
#estimateScore(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.input.gct"),
#              paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.score.gct"),
#              platform= "illumina")
#estimate.gct<-read.table(paste0("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/",Cancerset,"/TCGA.",Cancerset,".",Geneset,".estimate.score.gct"),skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe

# Calculate Immunoscore
RNASeq.subset.scaled <- scale (RNASeq.subset,scale=FALSE)
RNASeq.subset.scaled <- cbind(RNASeq.subset.scaled,rowMeans(RNASeq.subset.scaled))
colnames(RNASeq.subset.scaled)[ncol(RNASeq.subset.scaled)] <- c("scaled.IS")
immunoscore <- RNASeq.subset.scaled[,c("scaled.IS"),drop=FALSE]

# Culculate subscore ACT/INH
#Subset the Pancancer Matrix

RNAseq.INH.scaled <- scale(RNASeq.subset[,gene.list.INH],scale=FALSE)
RNAseq.ACT.scaled <- scale(RNASeq.subset[,gene.list.ACT],scale=FALSE)
immunoscore <- cbind(immunoscore,rowMeans(RNASeq.subset),rowMeans(RNAseq.ACT.scaled),rowMeans(RNAseq.INH.scaled))
colnames(immunoscore) <- c("scaled.IS","unscaled.IS","scaled.IS.ACT","scaled.IS.INH")
#estimate.gct<-t(estimate.gct[,-1])
#colnames(estimate.gct) <- estimate.gct[1,]
#estimate.gct<-estimate.gct[-1,]
#rownames(estimate.gct) <- gsub("\\.","-",rownames(estimate.gct))
rownames(immunoscore) <- gsub("\\.","-",rownames(immunoscore))
#mode(estimate.gct) <- "numeric"
#estimate.gct <- estimate.gct[rownames(immunoscore),]
#immunoscore <- cbind (immunoscore,estimate.gct)
write.csv (immunoscore,file=(paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.",DL.Method,".",Cancerset,".",Geneset,".csv")))

}

 