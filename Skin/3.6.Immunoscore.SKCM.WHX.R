#################################################################
###
### This script calculated the Immunoscore and 
### and TP53 frequency vs immunoscore
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

load ("./3 ANALISYS/Mutations/SKCM/Mutation.Data.split.RDATA")

Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/SKCM/SKCM.TCGA.EDASeq.k7.ISGS.reps5000/SKCM.TCGA.EDASeq.k7.ISGS.reps5000.k=4.consensusClass.ICR.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]

load ("./2 DATA/SUBSETS/SKCM/TCGA.SKCM.RNASeq.subset.ISGS.RData")
RNASeq.subset <- as.matrix(RNASeq.subset)

load ("./2 DATA/TCGA RNAseq/RNASeq_SKCM_EDASeq/SKCM.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
write.table(RNASeq.NORM_Log2,file="./2 DATA/TCGA RNAseq/RNASeq_SKCM_EDASeq/SKCM.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt",sep="\t",quote=FALSE)

# Calculate estimate score
filterCommonGenes(input.f="./2 DATA/TCGA RNAseq/RNASeq_SKCM_EDASeq/SKCM.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.txt",
                  output.f="./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.SKCM.estimate.input.gct",
                  id=c("GeneSymbol","EntrezID"))
estimateScore("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.SKCM.estimate.input.gct",
              "./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.SKCM.estimate.score.gct",
              platform= "illumina")
estimate.gct<-read.table("./3 ANALISYS/IMMUNOSCORE/ESTIMATE/TCGA.SKCM.estimate.score.gct",skip=2 , header = TRUE)#skip=2 tolgo le prime 2 righe

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
write.csv (immunoscore,file=("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.SKCM.ISGS.csv"))

# add immunoscore to mutation data
Mutation.All <- merge(Mutation.All,immunoscore[,"scaled.IS",drop=FALSE],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.All) <- Mutation.All$Row.names
Mutation.All$Row.names <- NULL
Mutation.Any <- merge(Mutation.Any,immunoscore[,"scaled.IS",drop=FALSE],by.x="Patient_ID",by.y="row.names",Any.x=TRUE, Any.y=FALSE)
row.names(Mutation.Any) <- Mutation.Any$Row.names
Mutation.Any$Row.names <- NULL
Mutation.Missense <- merge(Mutation.Missense,immunoscore[,"scaled.IS",drop=FALSE],by.x="Patient_ID",by.y="row.names",Missense.x=TRUE, Missense.y=FALSE)
row.names(Mutation.Missense) <- Mutation.Missense$Row.names
Mutation.Missense$Row.names <- NULL
Mutation.Nonsense <- merge(Mutation.Nonsense,immunoscore[,"scaled.IS",drop=FALSE],by.x="Patient_ID",by.y="row.names",Nonsense.x=TRUE, Nonsense.y=FALSE)
row.names(Mutation.Nonsense) <- Mutation.Nonsense$Row.names
Mutation.Nonsense$Row.names <- NULL
Mutation.Other <- merge(Mutation.Other,immunoscore[,"scaled.IS",drop=FALSE],by.x="Patient_ID",by.y="row.names",Other.x=TRUE, Other.y=FALSE)
row.names(Mutation.Other) <- Mutation.Other$Row.names
Mutation.Other$Row.names <- NULL


# order by ICR group and immunoscore 
Mutation.All <- Mutation.All[order(factor(Mutation.All$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting withing cluster add ,RNASeq.subset$scaled.IS ,-Mutation.All$scaled.IS

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

Mutation.Nonsense.TP53 <- Mutation.Nonsense[Mutation.Nonsense$Hugo_Symbol == "TP53",]
dim (Mutation.Nonsense.TP53) #64 TP53 mutations
Mutation.Nonsense.TP53 <- unique (Mutation.Nonsense.TP53) #64 mutations

Mutation.Other.TP53 <- Mutation.Other[Mutation.Other$Hugo_Symbol == "TP53",]
dim (Mutation.Other.TP53) #64 TP53 mutations
Mutation.Other.TP53 <- unique (Mutation.Other.TP53) #64 mutations

#TP53 imunoscore by cluster and by type of mutation

TP53.IS.All <- aggregate (Mutation.All.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.All.TP53$Cluster),FUN=mean)
colnames (TP53.IS.All) <- c("Cluster","TP53.All")
TP53.IS.Any <- aggregate (Mutation.Any.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Any.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Any) <- c("Cluster","TP53.Any")
TP53.IS.Missense <- aggregate (Mutation.Missense.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Missense.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Missense) <- c("Cluster","TP53.Missense")
TP53.IS.Nonsense <- aggregate (Mutation.Nonsense.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Nonsense.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Nonsense) <- c("Cluster","TP53.Nonsense")

#TP53.IS.Silent <- aggregate (Mutation.Silent[,c("scaled.IS"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean) #(only 1 silent mutations)
#colnames (TP53.IS.Silent) <- c("Cluster","TP53.Silent")
TP53.IS.Other <- aggregate (Mutation.Other[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Other$Cluster),FUN=mean)
colnames (TP53.IS.other) <- c("Cluster","TP53.Other")

TP53.IS <- cbind (TP53.IS.Any,TP53.IS.Missense,TP53.IS.Nonsense,TP53.IS.Other)
TP53.IS <- TP53.IS[,-c(3,5,7)]

print (TP53.IS)
write.csv (TP53.IS,file="./3 ANALISYS/Mutations/TP53mut_IS.scaled.csv")
 