#################################################################
###
### This script downloads :
### the ANY-CANCER RNASeq Data and Clinical data
### from the GDC database Using TCGA-BIOLINKS (L3 RNAseqV2 HiSeq)
### It will process the data into 1 table. 
### Data is saved :
### ../2 DATA/GDCdata/   (RAW DOWNLOAD)
### ../2 DATA/TCGA RNAseq/RNASeq_"Cancerset"_BIOLINKS/... (PROCESSED)
### ../2 DATA//Clinical Information/",Cancerset,"/BIOLINKS/... (indexed and XML versions)
### File to use :
### ../2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_BIOLINKS/Cancerset,".RNASeq.",download.source,".BIOLINKS.",sample.types,".DATA"
### Parameters to set : Cancerset, Download-source, sample.types
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
Cancerset           <- "GBM"
Profiling.method    <- "RNASeq"       #MA or RNASeq
sample.types        <- "Selected" #Alternatives TP , TP_TM , Selected

print (paste0("Downloading ",Cancerset," Data from ",download.source," using :"))
print (paste0("TCGAbiolinks version : ",packageVersion("TCGAbiolinks")))
print (paste0("SummarizedExperiment : ",packageVersion("SummarizedExperiment")))

# Paths and flies
Download.path       <- paste0("./2 DATA/",download.source," RNAseq/RNASeq_",Cancerset,"_BIOLINKS/")
Download.file       <- paste0(Cancerset,".RNASeq.",download.source,".BIOLINKS.",sample.types,".DATA")  

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
if (sample.types=="TP"){
  query <- GDCquery(project = paste0(download.source,"-",Cancerset),
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "HTSeq - Counts",
                    sample.type = c("Primary solid Tumor"))
}
if (sample.types=="TP_TM"){
  query <- GDCquery(project = paste0(download.source,"-",Cancerset),
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts",
                    sample.type = c("Primary solid Tumor","Metastatic"))
}
if (sample.types=="Selected"){
  query <- GDCquery(project = paste0(download.source,"-",Cancerset),
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts",
                    sample.type = c("Primary solid Tumor","Recurrent Solid Tumor","Additional - New Primary",
                                    "Metastatic","Additional Metastatic","Solid Tissue Normal"))
}

GDCdownload(query,
            directory="./2 DATA/GDCdata/",
            method = "api")
end.time <- Sys.time ()
time <- end.time - start.time 
print (time)

## prepare and save data into a SummarizedExperiment
dir.create(Download.path, showWarnings = FALSE)
data <- GDCprepare(query,
                   directory="./2 DATA/GDCdata/",
                   save = TRUE,
                   save.filename = paste0(Download.path,Download.file,".RDA"))
end.time <- Sys.time ()
time <- end.time - start.time 
print ("SummarizedExperiment created")
print (time)
#load(paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_BIOLINKS/",Cancerset,".RNASeq.TCGA.BIOLINKS.DATA.RDA"))

## Extract the Expression matrix
#Matrix <- assay(data)

## EDAseq normalisation using BIOLINKS
#dataPrep <- TCGAanalyze_Preprocessing(object = data,cor.cut = 0.6)
#dataNorm <- TCGAanalyze_Normalization(tabDF = Matrix,
#                                      geneInfo = geneInfo,
#                                      method = "gcContent")

# Clinical Data
dir.create(paste0("./2 DATA/Clinical Information/",Cancerset), showWarnings = FALSE)
dir.create(paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/"), showWarnings = FALSE)
## indexed GDC
clin <- GDCquery_clinic(paste0(download.source,"-",Cancerset), type = "clinical", save.csv = FALSE)
write.csv (clin,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/indexed.GDC.clinicaldata.csv"))
print ("Indexed Clinical ready")
## XML parsed
clin.query <- GDCquery(project = paste0(download.source,"-",Cancerset), 
                       data.category = "Clinical")
GDCdownload(clin.query,directory="./2 DATA/GDCdata/")
clinical.patient <- GDCprepare_clinic(clin.query,directory="./2 DATA/GDCdata/", clinical.info = "patient")
clinical.follup_up <- GDCprepare_clinic(clin.query,directory="./2 DATA/GDCdata/", clinical.info = "follow_up")
write.csv (clinical.patient,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Patient.XML.clinicaldata.csv"))
write.csv (clinical.follup_up,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Follow_up.XML.clinicaldata.csv"))
print ("XML Clinical ready")
## subtype data
Subtype.available <- c("ACC","BRCA", "COAD","GBM","HNSC","KICH","KIRP","KIRC","LGG","LUAD", "LUSC", "PRAD", "PANCAN","READ","SKCM","STAD","THCA","UCEC")
if (Cancerset %in% Subtype.available) {subtypes <- TCGAquery_subtype(tumor = Cancerset)
write.csv (subtypes,file=paste0("./2 DATA/Clinical Information/",Cancerset,"/BIOLINKS/Subtypes.clinicaldata.csv"))
print ("Subtype ready")
} else {print("No subtype data available")}
