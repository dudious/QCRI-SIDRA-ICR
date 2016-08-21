#################################################################
###
### This script downloads :
### the ANY-CANCER RNASeq Data 
### from the GDC database (L3 RNAseqV2 HiSeq)
### It will process the data into 1 table. 
### Data is saved :
### ../2 DATA/TCGA RNAseq/RNASeq_"Cancerset"_BIOLINKS/...
### File to use :
### NOT NORMALIZED 
### "Cancerset".RNASeq.TCGA.BIOLINKS.DATA.txt"
### NORMALIZED BY TCGA BIOLINKS
### "Cancerset".RNASeq.TCGA.BIOLINKS.DATA.GeneExp.rda"
###
### Parameters to set : Cancerset
###
#################################################################

# Setup environment
rm(list=ls())
#setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR")
setwd("/mnt3/wouter/BREAST-QATAR/")
## dependencies
#required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
#missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
#if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("TCGAbiolinks","SummarizedExperiment")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

library("TCGAbiolinks")
library("SummarizedExperiment")

# Set Parameters
download.source     <- "TCGA"
Cancerset           <- "BLCA"

# Paths and flies
Download.path       <- paste0("./2 DATA/",download.source," RNAseq/RNASeq_",Cancerset,"_BIOLINKS/")
Download.file       <- paste0(Cancerset,".RNASeq.",download.source,".BIOLINKS.DATA")  
  
# Download RNASeq data
start.time <- Sys.time ()

## LEGACY
#query <- GDCquery(project = paste0(download.source,"-",Cancerset),
#                  data.category = "Gene expression",
#                  data.type = "Gene expression quantification",
#                  legacy = TRUE,
#                  platform = "Illumina HiSeq",
#                  sample.type = "Primary solid Tumor")

## HARMONIZED
query <- GDCquery(project = paste0(download.source,"-",Cancerset),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Primary solid Tumor")

GDCdownload(query)
end.time <- Sys.time ()
time <- end.time - start.time #4.850641 mins
print (time)

## prepare and save data into a SummarizedExperiment
dir.create(Download.path, showWarnings = FALSE)
data <- GDCprepare(query,
                   save = TRUE,
                   save.filename = paste0(Download.path,Download.file,".RDA"))
#load(paste0("./2 DATA/TCGA RNAseq/RNASeq_OV_BIOLINKS/",Cancerset,".RNASeq.TCGA.BIOLINKS.DATA.RDA"))

## Extract the Expression matrix
Matrix <- assay(data,"raw_counts")

## EDAseq normalisation
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")
# Clinical Data
dir.create(paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/"), showWarnings = FALSE)
## indexed GDC
clin <- GDCquery_clinic(paste0(download.source,"-",Cancerset), type = "clinical", save.csv = FALSE)
write.csv (clin,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/indexed.GDC.clinicaldata.csv"))
## XML parsed
clin.query <- GDCquery(project = paste0(download.source,"-",Cancerset), 
                       data.category = "Clinical")
GDCdownload(clin.query)
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
clinical.follup_up <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
write.csv (clinical.patient,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Patient.XML.clinicaldata.csv"))
write.csv (clinical.follup_up,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Follow_up.XML.clinicaldata.csv"))

## subtype data
Subtype.available <- c("ACC","BRCA", "COAD","GBM","HNSC","KICH","KIRP","KIRC","LGG","LUAD", "LUSC", "PRAD", "PANCAN","READ","SKCM","STAD","THCA","UCEC")
if (Cancerset %in% Subtype.available) {subtypes <- TCGAquery_subtype(tumor = Cancerset)}



