###
# input variables for the subtype prediction script
# Source for this script TCGA, adjusted by WHX for our Paths and data
# Alternative for Genefu
###

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
required.packages <- c("heatmap.plus")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("ctc","impute")
 missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
 source("http://bioconductor.org/biocLite.R")
 if(length(missing.packages)) biocLite(missing.packages)
library(ctc)
library(heatmap.plus)

## RNASeq DAta from METABRIC
load ("./2 DATA/METABRIC/RNASEQ.DATA1.RData")
dim (RNASEQ.DATA)                                         #-> 997 25159

# Save processed data to pass on to IMS function
write.table (RNASEQ.DATA,file="./3 ANALISYS/IMS/TCGA IMS/RNASeq/BRCA.METABRIC.RNASeq.CLEANED.txt",sep = "\t",quote=FALSE,col.names=NA)

paramDir<- "./1 CODE/PAM50.TCGA.method/bioclassifier_R" # the location of unchanging files such as the function library and main program
inputDir<- "./3 ANALISYS/IMS/TCGA IMS/RNASeq/"  # the location of the data matrix, and where output will be located

inputFile<- "BRCA.METABRIC.RNASeq.CLEANED.txt" # the input data matrix as a tab delimited text file
short<-"BRCA.METABRIC.IMS.RNASeq.OUTPUT" # short name that will be used for output files

#MPD <- read.table ("./CODE/PAM50.TCGA.method/bioclassifier_R/mediansPerDataset_v2.1.txt") ## file source Katie Hoadley
calibrationParameters<- 9 	#the column of the "mediansPerDataset.txt" file to use for calibration; 
														              #NA will force centering within the test set & -1 will not do any 
														              #adjustment (when adjustment performed by used)

hasClinical<-FALSE 	#may include tumor size as second row, with 'T' as the gene name, 
										#and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
										#set this variable to FALSE if tumor size is not available

collapseMethod<-"mean" # can be mean or iqr (probe with max iqr is selected)
											# typically, mean is preferred for long oligo and
											# iqr is preferred for short oligo platforms


####
# run the assignment algorithm
####

source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))

TCGA.IMS <- read.table ("./3 ANALISYS/IMS/TCGA IMS/RNASeq/BRCA.TCGA.IMS.RNASeq.OUTPUT_pam50scores.txt",header=TRUE)
##TCGA.IMS <- TCGA.IMS[complete.cases(TCGA.IMS),] ## remove incomlete cases
rownames (TCGA.IMS) <- TCGA.IMS[,1]
TCGA.IMS[,1] <- NULL
print (table (TCGA.IMS$Call))
