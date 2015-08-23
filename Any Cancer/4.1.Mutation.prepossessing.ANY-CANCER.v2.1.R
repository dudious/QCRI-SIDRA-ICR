#################################################################
###
### This Script creates a merged myutation data file for 
### TCGA mutation data sources
### The order of selection is :
### Curated -> UCSC -> BI -> WUSM -> BCGSC ->BCM
### Input data :
### ./2 DATA/TCGA Mutations/,Cancerset,/Somatic_Mutations/.../*.maf
### Data is saved :
### ./2 DATA/TCGA Mutations/,Cancerset,/Somatic_Mutations/,Cancerset,.TCGA.combined.Mutation.Data.maf.Rdata
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

#Dependencies
required.packages <- c("Hmisc")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("Hmisc")

# Set Parameters
Cancerset <- "COAD-GA"

#read in available maf files
path <- paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations")
maf.all   <- list.files(path,recursive=TRUE,pattern=".maf",full.names=TRUE)
maf.total.number <- length(maf.all) - length(list.files(path,recursive=TRUE,pattern=".maf.Rdata",full.names=TRUE))
print (paste0("Total number of maf files available : ",maf.total.number))

# load curated data if availabe
maf.curated        <- list.files(paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations"),recursive=TRUE,pattern="curated",full.names=TRUE)
maf.curated.number <- length(maf.curated)
maf.curated.table  <- data.frame(Tumor_Sample_Barcode=character())
if (maf.curated.number==1){
  print ("1 curated maf file available")
  maf.curated.source <- strsplit (substring (maf.curated,nchar(path)+2,nchar(maf.curated)),"__")[[1]][[1]]
  print (paste0("Source : ",maf.curated.source))
  maf.curated.table <- cbind(read.csv(maf.curated,,sep ="\t",stringsAsFactors=FALSE),origin=maf.curated.source)
} else if (maf.curated.number>1){
  stop ("multiple curated maf files available")
} else if (maf.curated.number==0){
  print ("NO curated maf file available")
}

# load UCSC data if availabe
maf.UCSC        <- list.files(paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations"),recursive=TRUE,pattern="ucsc",full.names=TRUE)
maf.UCSC.number <- length(maf.UCSC)
maf.UCSC.table  <- data.frame(Tumor_Sample_Barcode=character())
if (maf.UCSC.number==1){
  print ("1 UCSC maf file available")
  maf.UCSC.source <- strsplit (substring (maf.UCSC,nchar(path)+2,nchar(maf.UCSC)),"__")[[1]][[1]]
  #print (paste0("Source : ",maf.UCSC.source))
  maf.UCSC.table <- cbind(read.csv(maf.UCSC,,sep ="\t",stringsAsFactors=FALSE),origin=maf.UCSC.source)
} else if (maf.UCSC.number>1){
  stop ("multiple UCSC maf files available")
} else if (maf.UCSC.number==0){
  print ("NO UCSC maf file available")
}

# load BI data if availabe
maf.BI        <- list.files(paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations"),recursive=TRUE,pattern="broad",full.names=TRUE)
maf.BI.number <- length(maf.BI)
maf.BI.table  <- data.frame(Tumor_Sample_Barcode=character())
if (maf.BI.number==1){
  print ("1 BI maf file available")
  maf.BI.source <- strsplit (substring (maf.BI,nchar(path)+2,nchar(maf.BI)),"__")[[1]][[1]]
  #print (paste0("Source : ",maf.BI.source))
  maf.BI.table <- cbind(read.csv(maf.BI,,sep ="\t",stringsAsFactors=FALSE),origin=maf.BI.source)
} else if (maf.BI.number>1){
  stop ("multiple BI maf files available")
} else if (maf.BI.number==0){
  print ("NO BI maf file available")
}

# load WUSM data if availabe
maf.WUSM        <- list.files(paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations"),recursive=TRUE,pattern="wust",full.names=TRUE)
maf.WUSM.number <- length(maf.WUSM)
maf.WUSM.table  <- data.frame(Tumor_Sample_Barcode=character())
if (maf.WUSM.number==1){
  print ("1 WUSM maf file available")
  maf.WUSM.source <- strsplit (substring (maf.WUSM,nchar(path)+2,nchar(maf.WUSM)),"__")[[1]][[1]]
  #print (paste0("Source : ",maf.WUSM.source))
  maf.WUSM.table <- cbind(read.csv(maf.WUSM,,sep ="\t",stringsAsFactors=FALSE),origin=maf.WUSM.source)
} else if (maf.WUSM.number>1){
  stop ("multiple WUSM maf files available")
} else if (maf.WUSM.number==0){
  print ("NO WUSM maf file available")
}

# load BCGSC data if availabe
maf.BCGSC        <- list.files(paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations"),recursive=TRUE,pattern="bcgsc",full.names=TRUE)
maf.BCGSC.number <- length(maf.BCGSC)
maf.BCGSC.table  <- data.frame(Tumor_Sample_Barcode=character())
if (maf.BCGSC.number==1){
  print ("1 BCGSC maf file available")
  maf.BCGSC.source <- strsplit (substring (maf.BCGSC,nchar(path)+2,nchar(maf.BCGSC)),"__")[[1]][[1]]
  #print (paste0("Source : ",maf.BCGSC.source))
  maf.BCGSC.table <- cbind(read.csv(maf.BCGSC,,sep ="\t",stringsAsFactors=FALSE),origin=maf.BCGSC.source)
} else if (maf.BCGSC.number>1){
  stop("multiple BCGSC maf files available")
} else if (maf.BCGSC.number==0){
  print ("NO BCGSC maf file available")
}

# load BCM data if availabe
maf.BCM        <- list.files(paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations"),recursive=TRUE,pattern="bcm",full.names=TRUE)
maf.BCM.number <- length(maf.BCM)
maf.BCM.table  <- data.frame(Tumor_Sample_Barcode=character())
if (maf.BCM.number==1){
  print ("1 BCM maf file available")
  maf.BCM.source <- strsplit (substring (maf.BCM,nchar(path)+2,nchar(maf.BCM)),"__")[[1]][[1]]
  #print (paste0("Source : ",maf.BCM.source))
  maf.BCM.table <- cbind(read.csv(maf.BCM,,sep ="\t",stringsAsFactors=FALSE),origin=maf.BCM.source)
} else if (maf.BCM.number>1){
  stop("multiple BCM maf files available")
} else if (maf.BCM.number==0){
  print ("NO BCM maf file available")
}


available.samples <- unique(c(maf.curated.table$Tumor_Sample_Barcode,
                              maf.UCSC.table$Tumor_Sample_Barcode,
                              maf.BI.table$Tumor_Sample_Barcode,
                              maf.WUSM.table$Tumor_Sample_Barcode,
                              maf.BCGSC.table$Tumor_Sample_Barcode,
                              maf.BCM.table$Tumor_Sample_Barcode))
print (paste0("Number of available patients : ",length(available.samples)))

# Add curated data
maf.merged.table <- maf.curated.table
print (paste0("Number of selected patients after curated is added : ",length(unique(maf.merged.table$Tumor_Sample_Barcode))))
# Add UCSC data
maf.merged.table<-rbind(maf.merged.table,maf.UCSC.table[substr(maf.UCSC.table$Tumor_Sample_Barcode,1,12) %nin% substr(maf.merged.table$Tumor_Sample_Barcode,1,12),])
print (paste0("Number of selected patients after UCSC is added : ",length(unique(maf.merged.table$Tumor_Sample_Barcode))))
# Add BI data
maf.merged.table<-rbind(maf.merged.table,maf.BI.table[substr(maf.BI.table$Tumor_Sample_Barcode,1,12) %nin% substr(maf.merged.table$Tumor_Sample_Barcode,1,12),])
print (paste0("Number of selected patients after BI is added : ",length(unique(maf.merged.table$Tumor_Sample_Barcode))))
# add WUSM data
maf.merged.table<-rbind(maf.merged.table,maf.WUSM.table[substr(maf.WUSM.table$Tumor_Sample_Barcode,1,12) %nin% substr(maf.merged.table$Tumor_Sample_Barcode,1,12),])
print (paste0("Number of selected patients after WUSM is added : ",length(unique(maf.merged.table$Tumor_Sample_Barcode))))
# add BCGSC data
maf.merged.table<-rbind(maf.merged.table,maf.BCGSC.table[substr(maf.BCGSC.table$Tumor_Sample_Barcode,1,12) %nin% substr(maf.merged.table$Tumor_Sample_Barcode,1,12),])
print (paste0("Number of selected patients after BCGSC is added : ",length(unique(maf.merged.table$Tumor_Sample_Barcode))))
# add BCM data
maf.merged.table<-rbind(maf.merged.table,maf.BCM.table[substr(maf.BCM.table$Tumor_Sample_Barcode,1,12) %nin% substr(maf.merged.table$Tumor_Sample_Barcode,1,12),])
print (paste0("Number of selected patients after BCM is added : ",length(unique(maf.merged.table$Tumor_Sample_Barcode))))

save (maf.merged.table,file=paste0(path,"/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))

