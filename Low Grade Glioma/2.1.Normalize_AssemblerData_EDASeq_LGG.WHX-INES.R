#################################################################
###
### This script normalizes the Glioblastoma cancer RNASeq Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/...
### File to use :
### "LGG.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"
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


# Load RAW Level 3 RNASeq Data downloaded with TCGA Assembler
RNASeq.DATA <- read.csv("./2 DATA/TCGA RNAseq/RNASeq_LGG_ASSEMBLER/LGG.RNASeq.TCGA.ASSEMBLER.DATA.txt",as.is=T,sep="\t")
#save(RNASeq.DATA,file="./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.DATA.RData") #Backup Source Data
#load("./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.DATA.RData") #Restore Source Data
View(head(RNASeq.DATA))
# Create Gene annotaion table split column 1 in 2 (Gene Name | EntrezID)
# RNASeq.Des <- read.csv("./2 DATA/TCGA RNAseq/RNASeq_LGG_ASSEMBLER/LGG.RNASeq.TCGA.ASSEMBLER.DATA.txt",as.is=T,sep="\t")
# RNASeq.DES <- as.data.frame(strsplit(RNASeq.Des[,1], "\\|")) 
# RNASeq.DES <- t(RNASeq.DES [,-1]) 
#  rownames(RNASeq.DES) <- NULL
#  colnames(RNASeq.DES) <- c("GeneSymbol","EntrezID")

RNASeq.DATA <- RNASeq.DATA[, c(1, grep("raw_count", RNASeq.DATA[1, ]))]#num of variables:94
rowNames <- sapply(strsplit(RNASeq.DATA[, 1], "|", fixed = T), function(x) ifelse(x[1] != "?", x[1], x[2]))
RNASeq.DATA<- RNASeq.DATA[-1, -1]
RNASeq.DATA <- apply(RNASeq.DATA, 2, as.numeric)
rownames(RNASeq.DATA) <- rowNames[2:length(rowNames)]
head(colnames(RNASeq.DATA))
colnames(RNASeq.DATA) <- gsub("\\.","-",colnames(RNASeq.DATA))
#View(RNASeq.DATA)
# dim(RNASeq.DATA)#20531    534
########################################################################################
# clean-up Data
# RNASeq.DATA <- RNASeq.DATA[-1,-c(1,2)]                                              #drop annotaion from data
# RNASeq.DATA <- RNASeq.DATA[,which(substring (colnames (RNASeq.DATA),28,31) == "7")] #remove "scaled_estimate" data
# colnames(RNASeq.DATA) <- gsub("\\.","-",colnames(RNASeq.DATA))                      #replace "." with "-" in sample names
# rownames(RNASeq.DATA) <- NULL
##################################################################################################

#save(RNASeq.DATA, file="./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.CLEANED.RData") #Backup clean Source Data
#load("./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.CLEANED.RData") #Restore clean Source Data

# Prepare data for normalization
#length(unique(colnames(RNASeq.DATA))):534
#length(unique(substring(colnames(RNASeq.DATA),1,12))):516
length(unique(colnames(RNASeq.DATA))) - length(unique(substring(colnames(RNASeq.DATA),1,12))) #18 patient with 2 samples (Tumor but different vial)
# tmp<-table(substring(colnames(RNASeq.DATA),1,12))
# tmp[which(tmp > 1)]#3
#table(colnames(RNASeq.DATA))
#TCGA-06-0190# TCGA-06-0190-01A-01R-1849-01 TCGA-06-0190-02A-01R-2005-01 
#TCGA-06-0211# TCCGA-06-0211-01A-01R-1849-01 TCGA-06-0211-01B-01R-1849-01
#TCGA-14-1034# TCGA-14-1034-01A-01R-1849-01 TCGA-14-1034-02B-01R-2005-01
#rm(tmp)
RNASeq.DATA <- ExtractTissueSpecificSamples(inputData = RNASeq.DATA,                          #this function will also remove anotation column
                                            tissueType = "TP",
                                            singleSampleFlag = TRUE,
                                            sampleTypeFile="./1 CODE/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt"); #select primary tumour
#516 TP tumors. NB: There is a tumor(01) duplicated that is heve been delated from using function

View(RNASeq.DATA[1:50,1:10])
#RNASeq.DATA <- cbind (RNASeq.DES,RNASeq.DATA)
length(unique(colnames(RNASeq.DATA))) - length(unique(substring(colnames(RNASeq.DATA),1,12))) #0 patient with 2 samples 
colnames(RNASeq.DATA) <- substring(colnames(RNASeq.DATA),1,12) #rename columns from sample to patient ID
dim (RNASeq.DATA)  #   20531   516                                    #-> 20531  373
#RNASeq.DATA <- RNASeq.DATA[which(rownames(RNASeq.DATA[,1]!="?"),]  #drop rows(genes) without approved gene-symbol
tmp <- table(rownames(RNASeq.DATA))
#tmp
tmp[which(tmp > 1)]
which(rownames(RNASeq.DATA)=="SLC35E2")
RNASeq.DATA[16302: 16303,1:10]
#which(RNASeq.DES[,1]=="SLC35E2")
#RNASeq.DES[16302: 16303,]
rownames(RNASeq.DATA)[16302] <- "SLC35E2B"#the EZ=728661 has been renamed in SLCE3B
RNASeq.DATA[16302: 16303,1:10]
#dim(RNASeq.DATA)#20531    83
EZid<-1:29
RNASeq.DATA<-RNASeq.DATA[-EZid,]
dim(RNASeq.DATA)#20502    516
View(head(RNASeq.DATA))

#levels (RNASeq.DATA[,1]) <- c(levels (RNASeq.DATA[,1]),"SLC35E2B") # Fix 1 duplicate gene with updated name
#RNASeq.DATA[16273,1] <- "SLC35E2B"#EZid=728661
#row.names(RNASeq.DATA) <- RNASeq.DATA[,1]                 #set row.names to GeneSymbol
#RNASeq.DATA <- RNASeq.DATA [,-c(1:2)]                     #Drop GeneSymbol, EntrezID, Hybridization REF columns
RNASeq.DATA <- as.matrix(RNASeq.DATA)
mode(RNASeq.DATA) <- "numeric"
RNASeq.DATA<-floor(RNASeq.DATA) #round off values lower integer
#View(head(RNASeq.DATA))
save(RNASeq.DATA,file="./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.PRENORM.RData") #Backup clean Source Data
#load("./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.PRENORM.RData") #Restore clean Source Data

# Normalization
load ("./2 DATA/geneInfo.RData")
geneInfo <- as.data.frame(geneInfo)
geneInfo <- geneInfo[rownames(RNASeq.DATA),] # drop the genes without RNAseq.DATA
dim(geneInfo)#20502     3
geneInfo <- geneInfo[!is.na(geneInfo[,1]),]                                # drop the genes without information from the RNAseq.DATA
dim(geneInfo)#20322     3
RNASeq.DATA <- RNASeq.DATA[rownames(geneInfo),] 
dim (RNASeq.DATA) #20322     516
RNASeq.expr.set <- newSeqExpressionSet(RNASeq.DATA, featureData = geneInfo)     # create a new SeqExpressionSet object.
dim (RNASeq.DATA)                                                             # -> 20322  83
fData(RNASeq.expr.set)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])    # make sure gcContenet is numeric
RNASeq.expr.set <- withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE) #removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set <- betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)             #removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM <-  log(RNASeq.DATA + .1) + offst(RNASeq.expr.set)                  
RNASeq.NORM <-  floor(exp(RNASeq.NORM) - .1)
#RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                       ## Alternative function using "preprocessCore" vs "stefanofunctions.R" , Process will take very long time
RNASeq.NORM.quantiles <- t(quantileNormalization(t(RNASeq.NORM)))
RNASeq.NORM <- floor(RNASeq.NORM.quantiles)
dim (RNASeq.NORM)#20322    516
# pdf("density_plot.pdf")
# plot(density(log(RNASeq.NORM[,1])))
# dev.off()
RNASeq.NORM_Log2<-log(RNASeq.NORM+1,2)
# pdf("boxplotNormCountsLog2.pdf")
# boxplot(RNASeq.NORM_Log2[,1:10])
# dev.off()

# Save Data
save(RNASeq.NORM,file="./2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")            #without log2 transformation
save(RNASeq.NORM_Log2,file="././2 DATA/TCGA RNAseq/RNASeq_LGG_EDASeq/LGG.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData") #with log2 transformation: the matrix to use

