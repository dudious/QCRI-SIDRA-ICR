#################################################################
###
### This script downloads :
### the ANY-CANCER RNASeq Data 
### from the TCGA database (L3 RNAseqV2 HiSeq)
### It will process the data into 1 table. 
### Data is saved :
### ../2 DATA/TCGA RNAseq/RNASeq_"Cancerset"_ASSEMBLER/...
### File to use :
### NOT NORMALIZED 
### "Cancerset".RNASeq.TCGA.ASSEMBLER.DATA.txt"
### NORMALIZED BY TCGA ASSEMBLER
### "Cancerset".RNASeq.TCGA.ASSEMBLER.DATA.GeneExp.rda"
###
### Parameters to set : Cancerset and TCGA.structure.file
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/Module_A.r")
source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/Module_B.r")

# Set Parameters
Cancerset           <- "KIRC"
TCGA.structure.file <- "./2 DATA/DirectoryTraverseResult_Jul-02-2015.rda"

# Paths and flies
Download.path       <- paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_ASSEMBLER/")
Download.file       <- paste0(Cancerset,".RNASeq.TCGA.ASSEMBLER.DATA")  
  
# Download data
start.time <- Sys.time ()
RNASeqRawData=DownloadRNASeqData(traverseResultFile=TCGA.structure.file,
                                 saveFolderName=Download.path,
                                 cancerType=Cancerset,
                                 assayPlatform="RNASeqV2",
                                 dataType="rsem.genes.results",
                                 outputFileName=Download.file)
end.time <- Sys.time ()
time <- end.time - start.time #4.850641 mins
print (time)

# rename files (NEEDS check for GA/Hiseq datasets(COAD,READ,UCEC))
file.hiseq <- list.files(Download.path,full.names = TRUE,pattern = "illuminahiseq" )
hiseq <- length(file.hiseq)
file.GA <- list.files(Download.path,full.names = TRUE,pattern = "illuminaga" )
GA <- length(file.GA)
if (GA==1){
  file.hiseq.new <- paste0(str_split(file.hiseq,paste0("DATA__",Cancerset,"__"))[[1]][[1]],"hiseq.DATA.txt")
  file.GA.new    <- paste0(str_split(file.hiseq,paste0("DATA__",Cancerset,"__"))[[1]][[1]],"GA.DATA.txt")
  file.list.old  <- c(file.hiseq,file.GA)
  file.list.new  <- c(file.hiseq.new,file.GA.new)
  tech="hiseq"
  Input.file     <- paste0(Download.path,Cancerset,".RNASeq.TCGA.ASSEMBLER.",tech,".DATA.txt")
  Download.file  <- paste0(Cancerset,".RNASeq.TCGA.ASSEMBLER.hiseq.DATA")  
}
if (GA==0){
  file.list.old  <- list.files(Download.path,full.names = TRUE,pattern = paste0("__",Cancerset,"__"))
  file.list.new  <- paste0(str_split(file.list.old,paste0("DATA__",Cancerset,"__"))[[1]][[1]],"DATA.txt")
  Input.file  <- paste0(Download.path,Cancerset,".RNASeq.TCGA.ASSEMBLER.DATA.txt")
}
file.rename (file.list.old,file.list.new)

# Process the downloaded normalized gene expression data and save the results
start.time  <- Sys.time ()
Output.file <- paste0(Download.file,".GeneExp")
GeneExpData <- ProcessRNASeqData(inputFilePath = Input.file,
                                outputFileName = Output.file,
                                outputFileFolder = Download.path,
                                dataType = "GeneExp",
                                verType = "RNASeqV2");
end.time <- Sys.time ()
time     <- end.time - start.time #4.850641 mins
print (time)

# Process the downloaded exon expression data and save the results
start.time <- Sys.time ()
Input.file  <- paste0(Download.path,Download.file,".GeneExp.txt")
Output.file <- paste0(Download.file,".ExonExp")
ExonExpData = ProcessRNASeqData(inputFilePath = Input.file,
                                outputFileName = Output.file,
                                outputFileFolder = Download.path,
                                dataType = "ExonExp",
                                verType = "RNASeqV2"); 
end.time <- Sys.time ()
time <- end.time - start.time #31.99029 secs
print (time)
 