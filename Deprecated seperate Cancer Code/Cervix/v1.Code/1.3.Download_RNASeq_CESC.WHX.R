#################################################################
###
### This script downloads :
### the Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] RNASeq Data 
### from the TCGA database (L3 CESC RNAseqV2 HiSeq)
### It will process the data into 1 table. 
### Data is saved :
### ../2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER/...
### File to use :
### NOT NORMALIZED 
### "CESC.RNASeq.TCGA.ASSEMBLER.DATA.txt"
### NORMALIZED BY TCGA ASSEMBLER
### "CESC.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp.rda"
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
source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")

getwd()

#Download data
RNASeqRawData=DownloadRNASeqData(traverseResultFile="./2 DATA/DirectoryTraverseResult_May-06-2015.rda",
                                 saveFolderName="./2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER/",
                                 cancerType= "CESC",
                                 assayPlatform="RNASeqV2",
                                 dataType="rsem.genes.results",
                                 outputFileName="CESC.RNASeq.TCGA.ASSEMBLER.DATA")

### manually rename file (remove _CESC_enc.edu.....)
cat ("RNAseq Data downloaded,...
        Please rename the file to : CESC.RNASeq.TCGA.ASSEMBLER.DATA.txt - remove _CESC_enc.edu -,...
        Press ENTER to continue... ")
line <- readline()

# Process the downloaded normalized gene expression data and save the results
start.time <- Sys.time ()
GeneExpData = ProcessRNASeqData(inputFilePath = "./2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER/CESC.RNASeq.TCGA.ASSEMBLER.DATA.txt",
                                outputFileName = "CESC.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp",
                                outputFileFolder = "./2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER",
                                dataType = "GeneExp",
                                verType = "RNASeqV2");
end.time <- Sys.time ()
time <- end.time - start.time #6.398157 mins
print (time)

# Process the downloaded exon expression data and save the results
start.time <- Sys.time ()
ExonExpData = ProcessRNASeqData(inputFilePath = "./2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER/CESC.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp.txt",
                                outputFileName = "CESC.RNASeq.TCGA.ASSEMBLER.DATA.ExonExp",
                                outputFileFolder = "./2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER",
                                dataType = "ExonExp",
                                verType = "RNASeqV2"); 
end.time <- Sys.time ()
time <- end.time - start.time #49.58051 secs
print (time)
 