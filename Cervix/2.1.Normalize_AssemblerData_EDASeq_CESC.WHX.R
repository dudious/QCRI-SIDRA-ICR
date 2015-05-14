#################################################################
###
### This script normalizes the Skin cancer RNASeq Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/...
### File to use :
### "CESC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
   required.packages.BioC <- c("EDASeq","base64enc","HGNChelper","RCurl","httr","stringr","digest","bitops")
   missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
   source("http://bioconductor.org/biocLite.R")
   if(length(missing.packages)) biocLite(missing.packages)

library("EDASeq","base64enc")
source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")
source("./1 CODE/R tools/stefanofunctions.R")

dir.create("./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/", showWarnings = FALSE)


# Load RAW Level 3 RNASeq Data downloaded with TCGA Assembler
RNASeq.DATA <- read.csv("./2 DATA/TCGA RNAseq/RNASeq_CESC_ASSEMBLER/CESC.RNASeq.TCGA.ASSEMBLER.DATA.txt",as.is=T,sep="\t")
#save(RNASeq.DATA,file="./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.DATA.RData") #Backup Source Data
#load("./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.DATA.RData") #Restore Source Data

# Create Gene annotaion table split column 1 in 2 (Gene Name | EntrezID)
RNASeq.DES <- as.data.frame(strsplit(RNASeq.DATA[,1], "\\|")) 
RNASeq.DES <- t(RNASeq.DES [,-1]) 
rownames(RNASeq.DES) <- NULL
colnames(RNASeq.DES) <- c("GeneSymbol","EntrezID")

# clean-up Data
RNASeq.DATA <- RNASeq.DATA[-1,-c(1,2)]                                              #drop annotaion from data
RNASeq.DATA <- RNASeq.DATA[,which(substring (colnames (RNASeq.DATA),28,31) == "7")] #remove "scaled_estimate" data
colnames(RNASeq.DATA) <- gsub("\\.","-",colnames(RNASeq.DATA))                      #replace "." with "-" in sample names
rownames(RNASeq.DATA) <- NULL

#save(RNASeq.DATA,RNASeq.DES,file="./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.CLEANED.RData") #Backup clean Source Data
#load("./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.CLEANED.RData") #Restore clean Source Data

# Prepare data for normalisation
length(unique(colnames(RNASeq.DATA))) - length(unique(substring(colnames(RNASeq.DATA),1,12))) #3 patient with 2 samples (normal controll tissue)
RNASeq.DATA <- ExtractTissueSpecificSamples(inputData = RNASeq.DATA,                          #this function will also remove anotation column
                                            tissueType = c("TP","TM"),
                                            singleSampleFlag = TRUE,
                                            sampleTypeFile="./1 CODE/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt"); #select primary tumour and metastatic only
RNASeq.DATA <- cbind (RNASeq.DES,RNASeq.DATA)
length(unique(colnames(RNASeq.DATA))) - length(unique(substring(colnames(RNASeq.DATA),1,12))) #0 patient with 2 samples (normal control tissue)
colnames(RNASeq.DATA) <- substring(colnames(RNASeq.DATA),1,12) #rename columns from sample to patient ID
dim (RNASeq.DATA)                                         #-> 20502  309
RNASeq.DATA <- RNASeq.DATA[which(RNASeq.DATA[,1]!="?"),]  #drop rows(genes) without approved gene-symbol
dim (RNASeq.DATA)    #-> 20502  309
#RNASeq.DATA[((which(duplicated(RNASeq.DATA[,1]))-1):which(duplicated(RNASeq.DATA[,1]))),1:10]
levels (RNASeq.DATA[,1]) <- c(levels (RNASeq.DATA[,1]),"SLC35E2B") # Fix 1 duplicate gene with updated name
RNASeq.DATA[RNASeq.DATA$EntrezID=="728661",1] <- "SLC35E2B"
row.names(RNASeq.DATA) <- RNASeq.DATA[,1]                 #set row.names to GeneSymbol
RNASeq.DATA <- RNASeq.DATA [,-c(1:2)]                     #Drop GeneSymbol, EntrezID, Hybridization REF columns
RNASeq.DATA <- as.matrix(RNASeq.DATA)
mode(RNASeq.DATA) <- "numeric"
RNASeq.DATA<-floor(RNASeq.DATA)                           #round off values lower integer
#save(RNASeq.DATA,RNASeq.DES,file="./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.PRENORM.RData") #Backup clean Source Data
#load("./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.PRENORM.RData") #Restore clean Source Data

# Normalisation
load ("./2 DATA/geneInfo.RData")                                                # file source : ???
geneInfo <- as.data.frame(geneInfo)
geneInfo <- geneInfo[rownames(RNASeq.DATA),]                                    # drop the genes without RNAseq.DATA
geneInfo <- geneInfo[!is.na(geneInfo[,1]),]                                     # drop the genes without information and convert to data frame
RNASeq.DATA <- RNASeq.DATA[rownames(geneInfo),]                                 # drop the genes without information from the RNAseq.DATA
RNASeq.expr.set <- newSeqExpressionSet(RNASeq.DATA, featureData = geneInfo)     # create a new SeqExpressionSet object.
dim (RNASeq.DATA)                                                               # -> 20322  371
fData(RNASeq.expr.set)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])    # make sure gcContenet is numeric
RNASeq.expr.set <- withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE) #removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set <- betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)             #removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM <-  log(RNASeq.DATA + .1) + offst(RNASeq.expr.set)                  
RNASeq.NORM <-  floor(exp(RNASeq.NORM) - .1)
#RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                       ## Alternative function using "preprocessCore" vs "stefanofunctions.R" , Process will take very long time
RNASeq.NORM.quantiles <- t(quantileNormalization(t(RNASeq.NORM)))
RNASeq.NORM <- floor(RNASeq.NORM.quantiles)
dim (RNASeq.NORM)

RNASeq.NORM_Log2<-log(RNASeq.NORM+1,2)

# Save Data
save(RNASeq.NORM,file="./2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")            #without log2 transformation
save(RNASeq.NORM_Log2,file="././2 DATA/TCGA RNAseq/RNASeq_CESC_EDASeq/CESC.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData") #with log2 transformation: the matrix to use

