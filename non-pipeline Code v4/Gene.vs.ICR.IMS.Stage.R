#########################
## Script to perform pecific gene barplot
## Input: cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
##         Number of clusters (num.clusters)
##         Paths to mutation file, cluster assignment file, and output filename
## 
######


## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")
#Dependencies
required.packages <- c("ggplot2", "plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")

## Parameters
Cancerset <- "BRCA"           
GOF = "CSPG4"

## Load Data
# Load normalised RNAseq data
load (paste0("./3_DataProcessing/TCGA_Assembler/",Cancerset,"/RNASeqData/BRCA_gene_RNAseq_normalized_TP_filtered.Rdata"))
# cluster assignment
load(paste0("./4_Analysis/TCGA_Assembler/",Cancerset,"/Clustering/",Cancerset,".TCGA_Assembler.EDASeq.ICR.reps5000/",Cancerset,"_ICR_cluster_assignment_k2-6.Rdata"))# select source data
# cliniacal data
# From RAW
Raw.Clinical <- read.csv (paste0("./2_Data/TCGA_Assembler/",Cancerset,"/BiospecimenClinicalData/patient.txt"),sep = "\t")
Raw.Clinical <- Raw.Clinical[-c(1,2),]
colnames(Raw.Clinical)
# From biolinks
IMS.data <- read.csv (paste0("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Subtypes.clinicaldata.csv"))

#Select gene
RNASeq.Data <- as.data.frame(t(log(filtered.norm.RNAseqData[GOF,,drop=FALSE],2)))

#Add Class to RNAseq data
RNASeq.Data$HML_Cluster <- factor(table_cluster_assignment$HML_cluster[match(rownames(RNASeq.Data),rownames(table_cluster_assignment))],levels = c("ICR Low","ICR Medium","ICR High"))
if (any(is.na(RNASeq.Data$HML_Cluster))){RNASeq.Data <-  RNASeq.Data[-which(is.na(RNASeq.Data$HML_Cluster)),]}
RNASeq.Data$Stage <- Raw.Clinical$ajcc_pathologic_tumor_stage[match(substr(rownames(RNASeq.Data),1,12),Raw.Clinical$bcr_patient_barcode)]
RNASeq.Data$IMS <- IMS.data$PAM50.mRNA[match(substr(rownames(RNASeq.Data),1,12),IMS.data$patient)]
RNASeq.Data$Race <- Raw.Clinical$race[match(substr(rownames(RNASeq.Data),1,12),Raw.Clinical$bcr_patient_barcode)]
colnames (RNASeq.Data) <- c("Gene","Cluster","Stage","IMS","Race")


#add IMS to RNAseq data

#blot
dev.new()
#dir.create (paste0("./4 FIGURES/expession geneBYCluster/",GOF,"/"),showWarnings=FALSE)
#png(paste0("./4 FIGURES/expession geneBYCluster/",GOF,"/",GOF,".",Cancerset,".",Geneset,".png", sep=""), height = 600, width = 600)
# ICR clsuter
#Anova
test.anova      = aov(Gene~Cluster,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = Cluster, y = Gene  )) +
       geom_boxplot() +
       labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "Cluster", y = "Log 2 Expression") +
       theme_classic()
#dev.off()


# Stage
dev.new()
#Anova
test.anova      = aov(Gene~Stage,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = Stage, y = Gene  )) +
  geom_boxplot() +
  labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "Stage", y = "Log 2 Expression") +
  theme_classic()


# IMS
dev.new()
# ICR filter
RNASeq.Data <- RNASeq.Data[RNASeq.Data$Cluster=="ICR High",]
#Anova
test.anova      = aov(Gene~IMS,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = IMS, y = Gene  )) +
  geom_boxplot() +
  labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "IMS", y = "Log 2 Expression") +
  theme_classic()

## histogram
ggplot(RNASeq.Data, aes( x=rownames(RNASeq.Data), y = Gene  )) +
  geom_histogram(stat = "Identity") +
  labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "IMS", y = "Log 2 Expression") +
  theme_classic()

# Race
dev.new()
# ICR filter
#RNASeq.Data <- RNASeq.Data[RNASeq.Data$Cluster=="ICR High",]
#Anova
test.anova      = aov(Gene~Race,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = Race, y = Gene  )) +
  geom_boxplot() +
  labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "Race", y = "Log 2 Expression") +
  theme_classic()

