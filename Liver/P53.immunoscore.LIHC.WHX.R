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
library("ggplot2")
library("plyr")

## Load Data

Mutation.data <- read.csv ("./2 DATA/TCGA Mutations WUSM/LIHC/LIHC.mutation.TCGA.txt",sep ="\t")
Patient_ID <- substr(Mutation.data$Tumor_Sample_Barcode,1,12)
Mutation.selected.data <- cbind(Patient_ID,Mutation.data[,c("Hugo_Symbol","Variant_Classification","Variant_Type")])
Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Hugo_Symbol == "TP53",]
dim (Mutation.selected.data) #64 TP53 mutations
Mutation.selected.data <- unique (Mutation.selected.data) #64 mutations
  
Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LIHC.TCGA.EDASeq.k4.ISGS.reps2000/LIHC.TCGA.EDASeq.k4.ISGS.reps2000.k=4.consensusClass.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]

load (paste0("./2 DATA/SUBSETS/LIHC/TCGA.LIHC.RNASeq.subset.ISGS.RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)

# Calculate Immunposcore
RNASeq.subset.scaled <- scale (RNASeq.subset,scale=FALSE)
RNASeq.subset.scaled <- cbind(RNASeq.subset.scaled,rowMeans(RNASeq.subset.scaled[,-ncol(RNASeq.subset.scaled)]))
colnames(RNASeq.subset.scaled)[ncol(RNASeq.subset.scaled)] <- c("avg")
immunoscore <- RNASeq.subset.scaled[,c("avg"),drop=FALSE]
write.csv (immunoscore,file=("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.LIHC.ISGS.csv"))

# Add Class to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,Consensus.class["Group"],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.selected.data) <- Mutation.selected.data$Row.names
Mutation.selected.data$Row.names <- NULL

# add immunoscore to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,immunoscore[,"avg",drop=FALSE],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.selected.data) <- Mutation.selected.data$Row.names
Mutation.selected.data$Row.names <- NULL

# order by ICR group and immunoscore 
Mutation.selected.data <- Mutation.selected.data[order(factor(Mutation.selected.data$Group,levels = c("ICR4","ICR3","ICR2","ICR1")),-Mutation.selected.data$avg),]    # if sorting withing cluster add ,RNASeq.subset$avg 

# split mutation types
Mutation.Missense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Missense_Mutation",]
rownames(Mutation.Missense) <- Mutation.Missense$Mutation_ID
Mutation.Missense$Mutation_ID <- NULL
Mutation.Missense$Variant_Classification <- NULL
Mutation.Missense$Variant_Type <- NULL

Mutation.Silent <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Silent",]
rownames(Mutation.Silent) <- Mutation.Silent$Mutation_ID
Mutation.Silent$Mutation_ID <- NULL
Mutation.Silent$Variant_Classification <- NULL
Mutation.Silent$Variant_Type <- NULL

Mutation.Nonsense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Nonsense_Mutation",]
rownames(Mutation.Nonsense) <- Mutation.Nonsense$Mutation_ID
Mutation.Nonsense$Mutation_ID <- NULL
Mutation.Nonsense$Variant_Classification <- NULL
Mutation.Nonsense$Variant_Type <- NULL

Mutation.other <- Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","RNA","Splice_Site","Silent"),]
rownames(Mutation.other) <- Mutation.other$Mutation_ID
Mutation.other$Variant_Classification <- NULL
Mutation.other$Variant_Type <- NULL
Mutation.other$Mutation_ID <- NULL

Mutation.any <- Mutation.selected.data
rownames(Mutation.any) <- Mutation.any$Mutation_ID
Mutation.any$Variant_Classification <- NULL
Mutation.any$Variant_Type <- NULL
Mutation.any$Mutation_ID <- NULL

#TP53 imunoscore by cluster and by type of mutation

TP53.IS.any <- aggregate (Mutation.any[,c("avg"),drop=FALSE],by=list(Mutation.any$Group),FUN=mean)
colnames (TP53.IS.any) <- c("Cluster","TP53.Any")
TP53.IS.Missense <- aggregate (Mutation.Missense[,c("avg"),drop=FALSE],by=list(Mutation.Missense$Group),FUN=mean)
colnames (TP53.IS.Missense) <- c("Cluster","TP53.Missense")
TP53.IS.Nonsense <- aggregate (Mutation.Nonsense[,c("avg"),drop=FALSE],by=list(Mutation.Nonsense$Group),FUN=mean)
colnames (TP53.IS.Nonsense) <- c("Cluster","TP53.Nonsense")
#TP53.IS.Silent <- aggregate (Mutation.Silent[,c("avg"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean) #(only 1 silent mutations)
#colnames (TP53.IS.Silent) <- c("Cluster","TP53.Silent")
TP53.IS.other <- aggregate (Mutation.other[,c("avg"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean)
colnames (TP53.IS.other) <- c("Cluster","TP53.other")

TP53.IS <- cbind (TP53.IS.any,TP53.IS.Missense,TP53.IS.Nonsense,TP53.IS.other)
TP53.IS <- TP53.IS[,-c(3,5,7)]
