#################################################################
###
### Create RNASEQ based Deconvolution
### and correlation with ICR
###
### Input files:
### ./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData
### "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"
### Output files:
### 
###
#################################################################

# Setup environment
rm(list=ls())

#setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    
#code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 

library(GSVA)
library(heatmap3)
library(corrplot)

# Set Parameters
Cancerset = "BRCA"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
ALL.HLA = c("HHLA1","HHLA2","HHLA3","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2",
            "HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DRB6","HLA-E","HLA-F","HLA-F-AS1","HLA-G","HLA-L")
## Normalised RNAseq data
Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancerset,"/RNASeqData")
if(Cancerset == "SKCM"){
  load(paste0(Cancer_path, "/", Cancerset, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
} else{
  load(paste0(Cancer_path, "/", Cancerset, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
}
## cluster assignment
load(paste0("./4_Analysis/TCGA_Assembler/",Cancerset,"/Clustering/",Cancerset,".TCGA_Assembler.EDASeq.ICR.reps5000/",Cancerset,"_ICR_cluster_assignment_k2-6.Rdata"))# select source data

# subset HLA genes
expression.subset <- filtered.norm.RNAseqData[ALL.HLA,]
table_cluster_assignment <- table_cluster_assignment[colnames(expression.subset),]
expression.subset <- rbind(expression.subset,ICR = table_cluster_assignment$Scaled_ICRscore)

correlations <- cor(t(expression.subset),method = "pearson")
dev.new()
corrplot(correlations,
         order="FPC")
