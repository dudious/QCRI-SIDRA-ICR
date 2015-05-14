#################################################################
###
### This script calculated the Immunoscore for TP53 mutated patients
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

dir.create("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/")

## Load Data

load ("./3 ANALISYS/Mutations/LIHC/Mutation.Data.split.RDATA")

Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LIHC/LIHC.TCGA.EDASeq.k7.ISGS.reps5000/LIHC.TCGA.EDASeq.k7.ISGS.reps5000.k=4.consensusClass.ICR.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]

load ("./2 DATA/SUBSETS/LIHC/TCGA.LIHC.RNASeq.subset.ISGS.RData")
RNASeq.subset <- as.matrix(RNASeq.subset)

load ("./2 DATA/TCGA RNAseq/RNASeq_LIHC_EDASeq/LIHC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
write.table(RNASeq.NORM_Log2,file="./2 DATA/TCGA RNAseq/RNASeq_LIHC_EDASeq/LIHC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt",sep="\t",quote=FALSE)

# Calculate estimate score
filterCommonGenes(input.f="./2 DATA/TCGA RNAseq/RNASeq_LIHC_EDASeq/LIHC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt",
                  output.f="./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.LIHC.estimate.input.gct",
                  id=c("GeneSymbol","EntrezID"))
estimateScore("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.LIHC.estimate.input.gct",
              "./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.LIHC.estimate.score.gct",
              platform= "illumina")
estimate.gct<-read.table("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.LIHC.estimate.score.gct",skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe

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
mode(estimate.gct) <- "numeric"
estimate.gct <- estimate.gct[rownames(immunoscore),]
immunoscore <- cbind (immunoscore,estimate.gct)
write.csv (immunoscore,file=("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.LIHC.ISGS.csv"))

# add immunoscore to mutation data
Mutation.All <- merge(Mutation.All,immunoscore[,"scaled.IS",drop=FALSE],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.All) <- Mutation.All$Row.names
Mutation.All$Row.names <- NULL

# order by ICR group and immunoscore 
Mutation.All <- Mutation.All[order(factor(Mutation.All$Group,levels = c("ICR4","ICR3","ICR2","ICR1")),-Mutation.All$scaled.IS),]    # if sorting withing cluster add ,RNASeq.subset$scaled.IS 

# Select TP53 mutations
Mutation.All.TP53 <- Mutation.All[Mutation.All$Hugo_Symbol == "TP53",]
dim (Mutation.All.TP53) #64 TP53 mutations
Mutation.All.TP53 <- unique (Mutation.All.TP53) #64 mutations

Mutation.Any.TP53 <- Mutation.Any[Mutation.Any$Hugo_Symbol == "TP53",]
dim (Mutation.Any.TP53) #64 TP53 mutations
Mutation.Any.TP53 <- unique (Mutation.Any.TP53) #64 mutations

Mutation.Missense.TP53 <- Mutation.Missense[Mutation.Missense$Hugo_Symbol == "TP53",]
dim (Mutation.Missense.TP53) #64 TP53 mutations
Mutation.Missense.TP53 <- unique (Mutation.Missense.TP53) #64 mutations

Mutation.Nonsense.TP53 <- Mutation.Nonsense[Mutation.Nonsense$Hugo_Symbol == "TP53",]
dim (Mutation.Nonsense.TP53) #64 TP53 mutations
Mutation.Nonsense.TP53 <- unique (Mutation.Nonsense.TP53) #64 mutations

#TP53 imunoscore by cluster and by type of mutation

TP53.IS.All <- aggregate (Mutation.All.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.All.TP53$Cluster),FUN=mean)
colnames (TP53.IS.any) <- c("Cluster","TP53.All")
TP53.IS.Any <- aggregate (Mutation.Any.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Any.TP53$Cluster),FUN=mean)
colnames (TP53.IS.any) <- c("Cluster","TP53.Any")
TP53.IS.Missense <- aggregate (Mutation.Missense.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Missense.TP53$Cluster),FUN=mean)
colnames (TP53.IS.any) <- c("Cluster","TP53.Missense")
TP53.IS.Nonsense <- aggregate (Mutation.Nonsense.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Nonsense.TP53$Cluster),FUN=mean)
colnames (TP53.IS.any) <- c("Cluster","TP53.Nonsense")

#TP53.IS.Silent <- aggregate (Mutation.Silent[,c("scaled.IS"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean) #(only 1 silent mutations)
#colnames (TP53.IS.Silent) <- c("Cluster","TP53.Silent")
TP53.IS.other <- aggregate (Mutation.other[,c("scaled.IS"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean)
colnames (TP53.IS.other) <- c("Cluster","TP53.other")

TP53.IS <- cbind (TP53.IS.any,TP53.IS.Missense,TP53.IS.Nonsense,TP53.IS.other)
TP53.IS <- TP53.IS[,-c(3,5,7)]
 