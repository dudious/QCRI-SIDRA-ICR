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
source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")

## RNASeq DAta from TCGA assembler
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_ASSEMBLER/BRCA.RNASeq.TCGA.ASSEMBLER.DATA.GeneExp.rda")
Des <- CheckGeneSymbol(Des)
RNASeq.Data <- cbind (Des,Data)
rm (Data,Des)
dim (RNASeq.Data)                                         #-> 20531  2433
RNASeq.Data <- RNASeq.Data[which(RNASeq.Data[,1]!="?"),]  #drop rows(genes) without approved genesymbol
dim (RNASeq.Data)                                         #-> 20502  2433 
row.names(RNASeq.Data) <- RNASeq.Data[,1]                 #set row.names to GeneSymbol
RNASeq.Data <- RNASeq.Data [,-c(1:3)]                     # Drop GeneSymbol, EntrezID, Hybridization REF columns
dim (RNASeq.Data)                                         #-> 20502  2430 
RNASeq.Data <- RNASeq.Data[,which(substring (colnames (RNASeq.Data),29,31) == "")] ##drop the ???? collumns
dim (RNASeq.Data)                                         #-> 20502  1215
mode (RNASeq.Data) <- "numeric"                           #convert the matrix to numeric


### extract patientID from sampleID and remove normal control samples
length(unique(colnames(RNASeq.Data))) - length(unique(substring(colnames(RNASeq.Data),1,12))) #120 patient with 2 samples (normal controll tissue)
RNASeq.Data.all <- RNASeq.Data
RNASeq.Data <- ExtractTissueSpecificSamples(inputData = RNASeq.Data.all,
                                             tissueType = "TP",
                                             singleSampleFlag = TRUE,
                                             sampleTypeFile="./1 CODE/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt");
length(unique(colnames(RNASeq.Data))) - length(unique(substring(colnames(RNASeq.Data),1,12))) #0 patient with 2 samples (normal controll tissue)
colnames(RNASeq.Data) <- substring(colnames(RNASeq.Data),1,12) #rename collumns from sample to patient ID


#rename to old genenames to fit with this method
rownames(RNASeq.Data)[rownames(RNASeq.Data) == "NDC80"] <- "KNTC2"
rownames(RNASeq.Data)[rownames(RNASeq.Data) == "NUF2"] <- "CDCA1"
rownames(RNASeq.Data)[rownames(RNASeq.Data) == "ORC6"] <- "ORC6L"

# LOG 2 transform RNASeq Data
RNASeq.Data<-log(RNASeq.Data+1,2)

# Save processed data to pass on to IMS function
write.table (RNASeq.Data,file="./3 ANALISYS/IMS/TCGA IMS/RNASeq/BRCA.TCGA.ASSEMBLER.RNASeq.CLEANED.txt",sep = "\t",quote=FALSE,col.names=NA)

paramDir<- "./1 CODE/PAM50.TCGA.method/bioclassifier_R" # the location of unchanging files such as the function library and main program
inputDir<- "./3 ANALISYS/IMS/TCGA IMS/RNASeq/"  # the location of the data matrix, and where output will be located

inputFile<- "BRCA.TCGA.ASSEMBLER.RNASeq.CLEANED.txt" # the input data matrix as a tab delimited text file
short<-"BRCA.TCGA.IMS.RNASeq.OUTPUT" # short name that will be used for output files

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
