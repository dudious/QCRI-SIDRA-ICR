#################################################################
###
### This script downloads :
### the ANY-CANCER Micro Array Data 
### from the TCGA database (L3 RNAseqV2 HiSeq)
### It will process the data into 1 table. 
### Data is saved :
### ../2 DATA/TCGA MA/MA_"Cancerset"_ASSEMBLER/...
### 
### File to use :
### "Cancerset".MA.TCGA.ASSEMBLER.CLEANED.Rdata"
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
source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/stefanofunctions.R")

# Set Parameters
Cancerset           <- "GBM"
TCGA.structure.file <- "./2 DATA/DirectoryTraverseResult_Jul-02-2015.rda"

# Paths and flies
Download.path       <- paste0("./2 DATA/TCGA MA/",Cancerset,"/")
Download.file       <- paste0(Cancerset,".MA.TCGA.ASSEMBLER.DATA")  
  
# Download data
start.time <- Sys.time ()
GeneExpData=DownloadRNASeqData(traverseResultFile=TCGA.structure.file,
                               saveFolderName=Download.path,
                               outputFileName=Download.file,
                               cancerType=Cancerset,
                               assayPlatform="Microarray")
end.time <- Sys.time ()
time <- end.time - start.time #4.850641 mins
print (time)

# rename files
file.agilent <- list.files(Download.path,full.names = TRUE,pattern = "agilent" )
agilent <- length(file.agilent)
file.affy <- list.files(Download.path,full.names = TRUE,pattern = "u133" )
affy <- length(file.affy)
 ## 1 file
if (agilent==1){
  file.agilent.new <- paste0(str_split(file.agilent,paste0("DATA__",Cancerset,"__"))[[1]][[1]],"agilent.DATA.txt")
}
if (affy==1){
  file.affy.new <- paste0(str_split(file.affy,paste0("DATA__",Cancerset,"__"))[[1]][[1]],"affy.DATA.txt")
}
 ## more then 1 file
if (agilent>1){
  file.agilent.new <- as.character(agilent)
  for (i in 1:agilent){
    file.agilent.new <- c(file.agilent.new,paste0(str_split(file.agilent[i],paste0("DATA__",Cancerset,"__"))[[1]][[1]],"agilent.",i,".DATA.txt"))
  }
  file.agilent.new <- file.agilent.new[-1]
}
file.list.old  <- c(file.agilent,file.affy)
file.list.new  <- c(file.agilent.new,file.affy.new)
file.rename (file.list.old,file.list.new)
file.source <- data.frame(cancer=character(),institute=character(),assay=character(),data.type=character(),download.date=character(),stringsAsFactors=FALSE)
 ## save source of data
for (i in 1: length(file.list.old)){
  file.source <- rbind (file.source,t(as.data.frame(str_split(file.list.old[i],"__")[[1]][-1],stringsAsFactors=FALSE)))
}
colnames(file.source) <- c("cancer","institute","assay","data.type","download.date")
rownames(file.source) <- NULL
write.csv (file.source,file=paste0(Download.path,Cancerset,".data.source.txt"))

# Process microarray data

if (agilent>0){
  for (i in 1:agilent){
   if (agilent==1) {
      i.text=""
    } else {
      i.text <- paste0(".",i)
    }
    print(paste0("Processing agilent data : ",i," of ",agilent," files"))
    input.file  <- paste0 (Download.path,Cancerset,".MA.TCGA.ASSEMBLER.agilent",i.text,".DATA.txt")
    output.file <- paste0(Cancerset,".MA.TCGA.ASSEMBLER.agilent",i.text,".PROCESSED")
    GeneExpData <- ProcessRNASeqData(inputFilePath =input.file,
                                     outputFileName = output.file,
                                     outputFileFolder = Download.path,
                                     dataType = "GeneExp",
                                     verType = "Microarray")
  ## transpose and cleanup 
  load(paste0(Download.path,output.file,".rda"))
  rownames(Data) <- Des[,1]
  agilentData <- data.frame(t(Data))
  rm("Data","Des")
  ## extract patientID from sampleID and remove normal control samples if present
  print(paste0(length(unique(rownames(agilentData)))," samples"))
  print(paste0(length(unique(rownames(agilentData))) - length(unique(substring(rownames(agilentData),1,12)))," Doubles"))
  agilentData <- ExtractTissueSpecificSamples(inputData = t(agilentData),    				#remove non tumour tissue data
                                              tissueType = "TP",
                                              singleSampleFlag = TRUE,
                                              sampleTypeFile="~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt");
  agilentData <- t(agilentData)
  print(paste0(length(unique(rownames(agilentData))) - length(unique(substring(rownames(agilentData),1,12)))," Doubles after cleanup"))
  ## extract patientID from sampleID
  agilentData<-cbind(substr(rownames(agilentData),1,12),agilentData)
  rownames(agilentData) <- agilentData[,1]
  agilentData <- agilentData[,-1]
  mode (agilentData) <- "numeric"   
  ## Save
  save (agilentData,file = paste0(Download.path,Cancerset,".MA.TCGA.ASSEMBLER.agilent",i.text,".CLEANED.Rdata"))
  #load (paste0(Download.path,Cancerset,".MA.TCGA.ASSEMBLER.agilent",i.text,".CLEANED.Rdata"))
  ## Quantile normalisation
  png(paste0(Download.path,Cancerset,".boxplot.MA.agilent",i.text,".preQNormalization.png"))
  boxplot(t(agilentData[1:15,]))
  dev.off()
    MA.quantiles<- quantileNormalization(agilentData)
  png(paste0(Download.path,Cancerset,".boxplot.MA.agilent",i.text,".postQNormalization.png"))
  boxplot(t(MA.quantiles[1:15,]))
  dev.off()
  
  save(MA.quantiles,file=paste0(Download.path,Cancerset,".MA.TCGA.ASSEMBLER.agilent",i.text,".QN.NORMALIZED.RData"))
  }
}

if (affy==1){
  print(paste0("Processing Affymetrix data : 1 of ",affy," files"))
  input.file  <- paste0 (Download.path,Cancerset,".MA.TCGA.ASSEMBLER.affy.DATA.txt")
  output.file <- paste0(Cancerset,".MA.TCGA.ASSEMBLER.affy.PROCESSED")
  GeneExpData <- ProcessRNASeqData(inputFilePath =input.file,
                                   outputFileName = output.file,
                                   outputFileFolder = Download.path,
                                   dataType = "GeneExp",
                                   verType = "Microarray")
  ## transpose and cleanup 
  load(paste0(Download.path,output.file,".rda"))
  rownames(Data) <- Des[,1]
  AffyData <- data.frame(t(Data))
  rm("Data","Des")
  ## extract patientID from sampleID and remove normal control samples if present
  print(paste0(length(unique(rownames(AffyData)))," samples"))
  print(paste0(length(unique(rownames(AffyData))) - length(unique(substring(rownames(AffyData),1,12)))," Doubles"))
  AffyData <- ExtractTissueSpecificSamples(inputData = t(AffyData),      			#remove non tumour tissue data
                                              tissueType = "TP",
                                              singleSampleFlag = TRUE,
                                              sampleTypeFile="~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/SupportingFiles/TCGASampleType.txt");
  AffyData <- t(AffyData)
  print(paste0(length(unique(rownames(AffyData))) - length(unique(substring(rownames(AffyData),1,12)))," Doubles after cleanup"))
  ## extract patientID from sampleID
  AffyData<-cbind(substr(rownames(AffyData),1,12),AffyData)
  rownames(AffyData) <- AffyData[,1]
  AffyData <- AffyData[,-1]
  mode (AffyData) <- "numeric"   
  ## Save
  save (AffyData,file = paste0(Download.path,Cancerset,".MA.TCGA.ASSEMBLER.affy.CLEANED.Rdata"))
  load (paste0(Download.path,Cancerset,".MA.TCGA.ASSEMBLER.affy.CLEANED.Rdata"))
  ## Quantile normalisation
  png(paste0(Download.path,Cancerset,".boxplot.MA.affy.preQNormalization.png"))
  boxplot(t(AffyData[1:15,]))
  dev.off()
    MA.quantiles<- quantileNormalization(AffyData)
  png(paste0(Download.path,Cancerset,".boxplot.MA.affy.postQNormalization.png"))
  boxplot(t(MA.quantiles[1:15,]))
  dev.off()
  
  save(MA.quantiles,file=paste0(Download.path,Cancerset,".MA.TCGA.ASSEMBLER.affy.QN.NORMALIZED.RData"))
}




