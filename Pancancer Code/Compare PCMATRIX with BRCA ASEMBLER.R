# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")

Consensus.class.ASSEMBLER <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv",header=TRUE) # select source data
Consensus.class.ASSEMBLER$X <- NULL

Consensus.class.PANCANCER <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/BRCA/BRCA.TCGA.PANCANCER.CLEAN.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.TCGA.PANCANCER.CLEAN.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv",header=TRUE) # select source data
Consensus.class.PANCANCER$X <- NULL
colnames(Consensus.class.PANCANCER) <- c("SampleID","ICR_CLASS")
Consensus.class.PANCANCER$PatientID <- substring(Consensus.class.PANCANCER$SampleID,1,12)

Consensus.class.PANCANCER$ICR_CLASS_PUB <- Consensus.class.ASSEMBLER$Group[match(Consensus.class.PANCANCER$PatientID,Consensus.class.ASSEMBLER$PatientID)]
Consensus.class.PANCANCER <- Consensus.class.PANCANCER[-which(is.na(Consensus.class.PANCANCER$ICR_CLASS_PUB)),]
Consensus.class.PANCANCER$overlap <- paste0(Consensus.class.PANCANCER$ICR_CLASS,"-",Consensus.class.PANCANCER$ICR_CLASS_PUB)
summary(as.factor(Consensus.class.PANCANCER$overlap))
dev.new()
barplot (summary(as.factor(Consensus.class.PANCANCER$overlap)))
