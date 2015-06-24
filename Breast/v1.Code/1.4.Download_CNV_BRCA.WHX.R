#################################################################
###
### This script downloads the Breast cancer copy number Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA TCGA CNV/TCGA CNV_BRCA_ASSEMBLER/...
### File to use :
### BRCA.CN.TCGA.ASSEMBLER.DATA.hg19.txt
### BRCA.GeneLevel.CNA.hg19.txt
###
#################################################################

#Download the CNA?CNV data using TCGA Assembler
# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")

getwd()

#Download data
start.time <- Sys.time ()
CNVRawData <- DownloadCNAData(traverseResultFile="./2 DATA/DirectoryTraverseResult_Jan-07-2015.rda",
                              saveFolderName="./2 DATA/TCGA CNV/TCGA CNV_BRCA_ASSEMBLER/",
                              cancerType= "BRCA",
                              assayPlatform="genome_wide_snp_6",
                              outputFileName="BRCA.CN.TCGA.ASSEMBLER.DATA")
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

### manually rename file (remove _BRCA_enc.edu.....)
cat ("Copy Number Data downloaded,...
        Please rename the file to : BRCA.CN.TCGA.ASSEMBLER.DATA.txt - remove _BRCA_enc.edu -,...
        Press ENTER when ready... ")
line <- readline()

# Process CNV data
start.time <- Sys.time ()
BRCA.GeneLevel.CNA = ProcessCNAData(inputFilePath = "./2 DATA/TCGA CNV/TCGA CNV_BRCA_ASSEMBLER/BRCA.CN.TCGA.ASSEMBLER.DATA.hg19.txt",
                                    outputFileName = "BRCA.GeneLevel.CNA.hg19",
                                    outputFileFolder = "./2 DATA/TCGA CNV/TCGA CNV_BRCA_ASSEMBLER/",
                                    refGenomeFile = "./1 CODE/R tools/TCGA-Assembler/SupportingFiles/Hg19GenePosition.txt");
end.time <- Sys.time ()
time <- end.time - start.time 
print (time) # 4.851125 hours
