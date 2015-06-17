#################################################################
###
### This script downloads ANY-CANCER copy number Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA TCGA CNV/TCGA CNV_"Cancerset"_ASSEMBLER/...
### File to use :
### "Cancerset".CN.TCGA.ASSEMBLER.DATA.hg19.txt
### "Cancerset".CN.TCGA.ASSEMBLER.DATA..GeneLevel.hg19.txt
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

# Set Parameters
Cancerset           <- "LGG"
TCGA.structure.file <- "./2 DATA/DirectoryTraverseResult_May-06-2015.rda"
GenomeFile = "./1 CODE/R tools/TCGA-Assembler/SupportingFiles/Hg19GenePosition.txt"

# Paths and flies
Download.path       <- paste0("./2 DATA/TCGA CNV/TCGA CNV_",Cancerset,"_ASSEMBLER/")
Download.file       <- paste0(Cancerset,".CN.TCGA.ASSEMBLER.DATA") 

#Download data
start.time <- Sys.time ()
CNVRawData <- DownloadCNAData(traverseResultFile=TCGA.structure.file,
                              saveFolderName=Download.path,
                              cancerType= Cancerset,
                              assayPlatform="genome_wide_snp_6",
                              outputFileName=Download.file)
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

# rename files (NEEDS change for GA/Hiseq datasets(COAD,READ,UCEC))
file.list.old  <- list.files(Download.path,full.names = TRUE,pattern = paste0("_",Cancerset,"_"))
file.list.new  <- paste0(str_split(file.list.old,paste0("DATA_",Cancerset,"_"))[[1]][[1]],"DATA.hg19.txt")
file.rename (file.list.old,file.list.new)

# Process CNV data
start.time <- Sys.time ()
Input.file  <- paste0(Download.path,Download.file,".hg19.txt")
Output.file <- paste0(Download.file,".GeneLevel.hg19")
BRCA.GeneLevel.CNA = ProcessCNAData(inputFilePath = Input.file,
                                    outputFileName = Output.file,
                                    outputFileFolder = Download.path,
                                    refGenomeFile = GenomeFile)
end.time <- Sys.time ()
time <- end.time - start.time 
print (time) # 4.851125 hours



