#################################################################
###
### This script downloads :
### the ANY-CANCER Mutation Data (MAF)
### from the GDC database Using TCGA-BIOLINKS 
### It will process the data into 1 table. 
### Data is saved :
### ../2 DATA/GDCdata/   (RAW DOWNLOAD)
### ../2 DATA/TCGA Mutations/BIOLINKS/RNASeq_"Cancerset"_BIOLINKS/... (PROCESSED)
### File to use :
### 
### Parameters to set : Cancerset, Download-source
###
#################################################################

# Setup environment
rm(list=ls())
#setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR")
#setwd("/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR")
setwd("/mnt3/wouter/BREAST-QATAR/")
## dependencies
#required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
#missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
#if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("TCGAbiolinks","SummarizedExperiment")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

library("curl")
library("SummarizedExperiment")
library("TCGAbiolinks")

# Set Parameters
download.source     <- "TCGA"
Cancersets          <- "ALL"

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}
  Parent.Cancerset <- substring(Cancerset,1,4)

# Download the MAF file
mut <- GDCquery_Maf(tumor = Cancerset)
save (mut,file=paste0("./2 DATA/TCGA Mutations/BIOLINKS/",Cancerset,".MAF.RData"))

# Load clinical data
clincal <- read.csv(file = paste0("./3 ANALISYS/CLINICAL DATA/",download.source,".",Cancerset ,".RNASeq_BIOLINKS_subset_clinicaldata.csv"))

# Load subtype data
Subtype.available <- c("ACC","BRCA", "COAD","GBM","HNSC","KICH","KIRP","KIRC","LGG","LUAD", "LUSC", "PRAD", "PANCAN","READ","SKCM","STAD","THCA","UCEC")
if (Cancerset %in% Subtype.available) {
  subtype <- read.csv(file = paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Subtypes.clinicaldata.csv"))
}

# load clustering Data
cluster <- read.csv(file = paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".",download.source,".BIOLINKS.EDASeq.k7.DBGS3.FLTR.reps5000/",
                                  Cancerset,".",download.source,".BIOLINKS.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv"))
rownames(cluster) <- cluster$PatientID
cluster$X <- NULL
colnames(cluster) <- c("bcr_patient_barcode","ICR.Cluster")

# Oncoprint c("MAP2K4","MAP3K1","TP53")
png(paste0("./4 FIGURES/Oncoprint/",Cancerset,".BIOLINKS.png"),res=600,height=6,width=6,unit="in")  # set filename
#dev.new()
TCGAvisualize_oncoprint(mut = mut, genes = c("MAP2K4","MAP3K1","TP53"),
                        filename = NULL,
                        annotation = cluster,
                        color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
                        #column_order = rownames(cluster),
                        rows.font.size=10,
                        heatmap.legend.side = "right",
                        dist.col = 0,
                        label.font.size = 10)
dev.off()
print(paste0(Cancerset," Done..."))
Cancersets = Cancersets[-i,]
}
