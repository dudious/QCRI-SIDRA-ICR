#################################################################
###
### This script downloads the Breast cancer Micro Array Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA BC MA/...
### File to use :
### "BRCA.MA.TCGA.ASSEMBLER.CLEANED.Rdata"
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
  RNASeqRawData=DownloadRNASeqData(traverseResultFile="./2 DATA/DirectoryTraverseResult_Jan-07-2015.rda",
                                 saveFolderName="./2 DATA/TCGA BC MA/",
                                 cancerType= "BRCA",
                                 assayPlatform="Microarray",
                                 outputFileName="BRCA.MA.TCGA.ASSEMBLER.DATA");
  ### manually rename file (remove _BRCA_enc.edu.....)
  cat ("Micro Array downloaded,...
        Please rename the file to : BRCA.MA.TCGA.ASSEMBLER.DATA.txt - remove _BRCA_enc.edu -,...
        Press ENTER to continue... ")
  line <- readline()

# Process microarray data
  GeneExpData <- ProcessRNASeqData(inputFilePath ="./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.DATA.txt",
                                   outputFileName = "BRCA.MA.TCGA.ASSEMBLER.PROCESSED",
                                   outputFileFolder = "./2 DATA/TCGA BC MA/",
                                   dataType = "GeneExp",
                                   verType = "Microarray");
  
# transpose and cleanup 
  load("./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.PROCESSED.rda")
  rownames(Data) <- Des[,1]
  agilentData <- data.frame(t(Data))
  rm("Data","Des")
  ## extract patientID from sampleID and remove normal control samples if present
  length(unique(rownames(agilentData))) #604 samples
  length(unique(rownames(agilentData))) - length(unique(substring(rownames(agilentData),1,12))) #71 patient with 2 samples (normal control tissue)
  agilentData <- ExtractTissueSpecificSamples(inputData = t(agilentData),						#remove non tumour tissue data
                                                         tissueType = "TP",
                                                         singleSampleFlag = TRUE,
                                                         sampleTypeFile="./1 CODE/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt");
  agilentData <- t(agilentData)
  length(unique(rownames(agilentData))) - length(unique(substring(rownames(agilentData),1,12))) #0 patient with 2 samples 531 samples left
  ## extract patientID from sampleID
  agilentData<-cbind(substr(rownames(agilentData),1,12),agilentData)
  rownames(agilentData) <- agilentData[,1]
  agilentData <- agilentData[,-1]
  mode (agilentData) <- "numeric"   

# Save
  save (agilentData,file = "./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.Rdata")
  