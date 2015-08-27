#########################
## Script to perform pecific gene barplot
## Input: Mutation .Frequencies.RDATA file, and the cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
##         Number of clusters (num.clusters)
##         Paths to mutation file, cluster assignment file, and output filename
## 
######


## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("ggplot2", "plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")

## Parameters
Cancerset <- "LGG"           
BRCA.Filter <- "BSF2"          # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
GOF = "HLA-G"

## Load log2 transformed normalised RNAseq data
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))

#cluster assignment
if (Cancerset == "BRCA"){
  if (substring(Geneset,7,10)=="FLTR"){
    Cancerset <- paste0(Cancerset,".",BRCA.Filter)
  }
}
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#Select gene
RNASeq.Data <- t(RNASeq.NORM_Log2)
RNASeq.Data <- as.data.frame(RNASeq.Data[,GOF,drop=FALSE])

#Add Class to RNAseq data
RNASeq.Data$Cluster <- Consensus.class$Cluster[match(rownames(RNASeq.Data),rownames(Consensus.class))]
RNASeq.Data <-  RNASeq.Data[-which(is.na(RNASeq.Data$Cluster)),]
colnames (RNASeq.Data) <- c("Gene","Cluster")

#Anova
test.anova      = aov(Gene~Cluster,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
    p.value.rounded = paste0("p = ",round(p.value,4))
}

#blot
dir.create (paste0("./4 FIGURES/expession geneBYCluster/",GOF,"/"),showWarnings=FALSE)
png(paste0("./4 FIGURES/expession geneBYCluster/",GOF,"/",GOF,".",Cancerset,".",Geneset,".png", sep=""), height = 600, width = 600)
 gg = ggplot(RNASeq.Data, aes(x = factor(Cluster), y = Gene  )) +
              geom_boxplot() +
              ggtitle (paste0(GOF," in ",Cancerset," (",p.value.rounded,")"))
              
              
 print(gg)
dev.off()





