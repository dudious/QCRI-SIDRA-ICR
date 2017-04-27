# Setup environment
rm(list=ls())
#setwd("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
#setwd("/mnt3/wouter/BREAST-QATAR/")
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")

#load data
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.BIOLINKS.Selected.NORMALIZED.TP_FILTERED_LOG2.RData")

write.csv(RNASeq.NORM.TP_Log2,file="./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.BIOLINKS.Selected.NORMALIZED.TP_FILTERED_LOG2.csv",quote = FALSE)
