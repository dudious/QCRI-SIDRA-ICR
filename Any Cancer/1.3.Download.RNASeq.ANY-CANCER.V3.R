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
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
#required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
#missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
#if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("TCGAbiolinks")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

library("TCGAbiolinks")

# Set Parameters
download.source     <- "TCGA"
Cancerset           <- "OV"

# Paths and flies
Download.path       <- paste0("./2 DATA/",download.source," RNAseq/RNASeq_",Cancerset,"_BIOLINKS/")
Download.file       <- paste0(Cancerset,".RNASeq.",download.source,".BIOLINKS.DATA")  
  
# Download data
start.time <- Sys.time ()
query <- GDCquery(project = paste0(download.source,"-",Cancerset),
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  legacy = TRUE,
                  platform = "Illumina HiSeq",
                  sample.type = "Primary solid Tumor")
GDCdownload(query)
end.time <- Sys.time ()
time <- end.time - start.time #4.850641 mins
print (time)

# prepare and save data

data <- GDCprepare(query,
                   save = TRUE,
                   save.filename = paste0(Download.path,Download.file,".RDA"))

Matrix <- assay(data,"raw_counts")
