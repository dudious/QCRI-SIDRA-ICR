#################################################################
###
### This script compares BRCA clusterings using diffenfilters
###
### 
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")

# Set Parameters
Cancersets        = "ALL"
Geneset           = "DBGS3"   
DL.Method         = "PANCANCER.CLEAN"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Split"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer

ICR.K4.Assembler.BSF2 <- read.csv ("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv")
ICR.K4.Assembler.PCF <- read.csv ("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.PCF/BRCA.PCF.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.PCF.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv")
ICR.K4.Biolinks <-read.csv ("./3 ANALISYS/CLUSTERING/RNAseq/BRCA/BRCA.TCGA.BIOLINKS.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.TCGA.BIOLINKS.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv")
ICR.K4.Pancancer.clean <-read.csv ("./3 ANALISYS/CLUSTERING/RNAseq/BRCA/BRCA.TCGA.PANCANCER.CLEAN.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.TCGA.PANCANCER.CLEAN.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv")

ICR.matching <- ICR.K4.Pancancer.clean
rownames(ICR.matching) <- substring(ICR.matching$PatientID,1,12)
ICR.matching$X <-NULL
colnames(ICR.matching) <- c("Sample_ID","Pancancer.clean.k4")
ICR.matching$Biolinks.k4 <- ICR.K4.Biolinks$Group[match(rownames(ICR.matching),ICR.K4.Biolinks$PatientID)]
ICR.matching$Assembler.PCF.k4 <- ICR.K4.Assembler.PCF$Group[match(rownames(ICR.matching),ICR.K4.Assembler.PCF$PatientID)]
ICR.matching$Assembler.BSF2.k4 <- ICR.K4.Assembler.BSF2$Group[match(rownames(ICR.matching),ICR.K4.Assembler.BSF2$PatientID)]

not.matching <- ICR.matching[which (ICR.matching$Pancancer.clean.k4 != ICR.matching$Assembler.BSF2.k4),]
not.matching.1.4 <- ICR.matching[which (ICR.matching$Pancancer.clean.k4 != ICR.matching$Assembler.BSF2.k4 & ICR.matching$Pancancer.clean.k4 %in% c("ICR1","ICR4")),]

write.csv (not.matching,"./3 ANALISYS/CLUSTERING/RNAseq/BRCA/matching_between_clustering.csv")
