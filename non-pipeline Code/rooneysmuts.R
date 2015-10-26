# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

# Read Data
Rooney.data <- read.csv ("./2 DATA/Rooneys.mutation.data.txt",sep="")
cluster.assignment <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv")

Rooney.data$Cluster <- cluster.assignment$Group[match(Rooney.data$sample,cluster.assignment$PatientID)]
Rooney.data.BRCA <- Rooney.data[-which(is.na(Rooney.data$Cluster)),]
Rooney.data.BRCA.HLA <- unique(Rooney.data.BRCA[,c("sample","hla")])
Rooney.data.BRCA.samples <- unique(Rooney.data.BRCA[,c("sample"),drop=FALSE])
