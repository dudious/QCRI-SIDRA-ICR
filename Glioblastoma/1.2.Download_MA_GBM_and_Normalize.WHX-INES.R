#################################################################
###
### This script downloads the Glioblastoma cancer Micro Array Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA BC MA/...
### File to use :
### "GBM.MA.TCGA.ASSEMBLER.CLEANED.Rdata"
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR")
  ## dependencies
  ## install java for xlsx export
  ## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
     required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
     missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
     if(length(missing.packages)) install.packages(missing.packages)
  source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
  source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")

# Download microarray data
  RNASeqRawData=DownloadRNASeqData(traverseResultFile="./2 DATA/DirectoryTraverseResult_May-06-2015.rda",
                                 saveFolderName="./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/",
                                 cancerType= "GBM",
                                 assayPlatform="Microarray",
                                 outputFileName="GBM.MA.TCGA.ASSEMBLER.DATA");
  ### manually rename file (remove _GBM_enc.edu.....)
#   cat ("Micro Array downloaded,...
#         Please rename the file to : GBM.MA.TCGA.Affy.ASSEMBLER.DATA.txt - remove _GBM_enc.edu -,...
#         Press ENTER to continue... ")
  #line <- readline()

# Process microarray data
  GeneExpData <- ProcessRNASeqData(inputFilePath ="./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.DATA.txt.txt",
                                   outputFileName = "GBM.MA.Affy.TCGA.ASSEMBLER.PROCESSED",
                                   outputFileFolder = "./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/",
                                   dataType = "GeneExp",
                                   verType = "Microarray");
  
# transpose and cleanup 
  load("./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.PROCESSED.rda")
  rownames(Data) <- Des[,1]
  #agilentData <- data.frame(t(Data))
  View(head(Data))
  AffyData<-data.frame(t(Data))
  dim(AffyData)#558 12042
  rm("Data","Des")
  ## extract patientID from sampleID and remove normal control samples if present
  length(unique(rownames(AffyData))) #558 samples
  length(unique(rownames(AffyData))) - length(unique(substring(rownames(AffyData),1,12))) #19 patient with 2 samples 
  AffyData <- ExtractTissueSpecificSamples(inputData = t(AffyData),						#remove non tumour tissue data
                                                         tissueType = "TP",
                                                         singleSampleFlag = TRUE,
                                                         sampleTypeFile="./1 CODE/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt");# #529 samples of TP are extracted.
  AffyData <- t(AffyData)
  dim(AffyData)#529 12042
  length(unique(rownames(AffyData))) - length(unique(substring(rownames(AffyData),1,12))) #0 patient with 2 samples 
  ## extract patientID from sampleID
  AffyData<-cbind(substr(rownames(AffyData),1,12),AffyData)
  rownames(AffyData) <- AffyData[,1]
  AffyData <- AffyData[,-1]
  mode (AffyData) <- "numeric"   

# Save
  save (AffyData,file = "./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.CLEANED.Rdata")#cleaned and preQN
 ###########################################################################
#Normaization MAData
load("./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.CLEANED.Rdata")
source("./1 CODE/R tools/stefanofunctions.R")
dim(AffyData)#campione x geni
png("./4 FIGURES/Boxplot-Exp-Matrix/boxplot.MA.Affy.GBM.preQNormalization.png")
boxplot(t(AffyData[1:15,]))
dev.off()

AffyData.quantiles<- quantileNormalization(AffyData)
png("./4 FIGURES/Boxplot-Exp-Matrix/boxplot.MA.Affy.GBM.postQNormalization.png")
boxplot(t(AffyData.quantiles[1:15,]))
dev.off()
AffyData.NORM<-AffyData.quantiles
dim(AffyData.NORM)#529 12042

save(AffyData.NORM,file="./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.RData")
