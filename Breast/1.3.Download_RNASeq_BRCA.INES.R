#################################################################
###
### This script downloads the Breast cancer RNASeq Data 
### from the TCGA database (L3 BRCA RNAseqV2, HiSeq)
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER/...
### File to use :
### NOT NORMALIZED 
### "BRCA.RNASeq.TCGA.ASSEMBLER.DATA.txt"
### NORMALIZED BY TCGA ASSEMBLER
### "BRCA.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp.rda"
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
RNASeqRawData=DownloadRNASeqData(traverseResultFile="./2 DATA/DirectoryTraverseResult_Jan-07-2015.rda",
                                 saveFolderName="./2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER/",
                                 cancerType= "BRCA",
                                 assayPlatform="RNASeqV2",
                                 dataType="rsem.genes.results",
                                 outputFileName="BRCA.RNASeq.TCGA.ASSEMBLER.DATA")

### manually rename file (remove _BRCA_enc.edu.....)
cat ("Micro Array downloaded,...
        Please rename the file to : BRCA.RNASeq.TCGA.ASSEMBLER.DATA.txt - remove _BRCA_enc.edu -,...
        Press ENTER to continue... ")
line <- readline()

# Process the downloaded normalized gene expression data and save the results
start.time <- Sys.time ()
GeneExpData = ProcessRNASeqData(inputFilePath = "./2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER/BRCA.RNASeq.TCGA.ASSEMBLER.DATA.txt",
                                outputFileName = "BRCA.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp",
                                outputFileFolder = "./2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER",
                                dataType = "GeneExp",
                                verType = "RNASeqV2");
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

# Process the downloaded exon expression data and save the results
start.time <- Sys.time ()
ExonExpData = ProcessRNASeqData(inputFilePath = "./2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER/BRCA.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp.txt",
                                outputFileName = "BRCA.RNASeq.TCGA.ASSEMBLER.DATA.ExonExp",
                                outputFileFolder = "./2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER",
                                dataType = "ExonExp",
                                verType = "RNASeqV2"); 
end.time <- Sys.time ()
time <- end.time - start.time
print (time)
 