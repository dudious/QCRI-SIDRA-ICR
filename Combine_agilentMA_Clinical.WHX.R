#################################################################
###
### This Script combines the clinical data relating to patients 
### that have Micro Array Data available with IMS predictions based 
### based on MicroArray data 
### (Genfu.PAM50,Genefu.Sorlie,Genefu.Hu and TCGA.Method.PAM50.)
### It will aslo add exclusion flag
### a "Master file" is created merging MA and clinical data
### Data is saved in :
### "./3 ANALISYS/CLINICAL DATA/"
### File to use :
### "Agilent_subset_clinicaldata.csv"
### Master file :
### "Agilent_subset_clinical_AND_MA_DATA.csv"
###
#################################################################

# Setup environment
  rm(list=ls())
  ## dependencies
  required.packages <- c("xlsx","Hmisc")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  library (xlsx) #xlsx needs java installed
  library (Hmisc)
  setwd("~/Dropbox/BREAST_QATAR/")

# Load data files 
  ## agilent DAta from TCGA ASSEMBLER
  load ("./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.Rdata")
  
  ## clinical data retrieved through TCGA assembler
  clinicalData <- read.csv ("./2 DATA/Clinical Information/BRCA/selected_clinicaldata.txt", header = TRUE)
  ## IMS 1 ### prediction data from TCGA ( "https://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt")
  IMS.TCGA<-read.csv("./2 DATA/Clinical Information/Downloaded files/BRCA.547.PAM50.SigClust.Subtypes.txt",sep="\t",header = TRUE)
    ### extract patientID from sampleID
    IMS.TCGA<-cbind(IMS.TCGA,substr(IMS.TCGA$Sample,1,12))
    col.names <- colnames(IMS.TCGA)
    col.names[length(col.names)] <- "PatientID"
    colnames(IMS.TCGA) <- col.names
    ### remove non-tumor tissue samples (doubles)
    IMS.TCGA<-subset (IMS.TCGA,IMS.TCGA$Type == "tumor")
    ### Set rownames
    rownames(IMS.TCGA) <-IMS.TCGA$PatientID
    ### rename PAM50 names
    levels(IMS.TCGA$PAM50) <- c(levels(IMS.TCGA$PAM50), c("Basal-like","HER2-enriched","Luminal A","Luminal B","Normal-like")) 
    IMS.TCGA$PAM50[IMS.TCGA$PAM50=="Basal"]  <- "Basal-like"
    IMS.TCGA$PAM50[IMS.TCGA$PAM50=="Her2"]   <- "HER2-enriched"
    IMS.TCGA$PAM50[IMS.TCGA$PAM50=="LumA"]   <- "Luminal A"
    IMS.TCGA$PAM50[IMS.TCGA$PAM50=="LumB"]   <- "Luminal B"
    IMS.TCGA$PAM50[IMS.TCGA$PAM50=="Normal"] <- "Normal-like"
    IMS.TCGA$PAM50 <- droplevels(IMS.TCGA$PAM50)
  ## IMS 2 ### prediction data from TCGA.PAM50.METHOD using Agilent data through ASSEMBLER
  IMS.TCGA.METHOD<-read.csv("./3 ANALISYS/IMS/TCGA IMS/MA/BRCA.TCGA.MA.IMS.OUTPUT_pam50scores.txt",sep="\t",header = TRUE)
  IMS.TCGA.METHOD <- IMS.TCGA.METHOD[complete.cases(IMS.TCGA.METHOD),] ## remove incomlete cases
  rownames (IMS.TCGA.METHOD) <- IMS.TCGA.METHOD[,1]
  IMS.TCGA.METHOD[,1] <- NULL
  ### rename PAM50 names
  levels(IMS.TCGA.METHOD$Call) <- c(levels(IMS.TCGA.METHOD$Call), c("Basal-like","HER2-enriched","Luminal A","Luminal B","Normal-like")) 
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="Basal"]  <- "Basal-like"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="Her2"]   <- "HER2-enriched"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="LumA"]   <- "Luminal A"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="LumB"]   <- "Luminal B"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="Normal"] <- "Normal-like"
  IMS.TCGA.METHOD$Call <- droplevels(IMS.TCGA.METHOD$Call)
  ## IMS 3,4,5 ###prediction from Genefu using Agilent data through ASSEMBLER
  IMSData <- read.csv ("./3 ANALISYS/IMS/Genefu IMS/PredictionTable.ASSEMBLER.csv", header = TRUE)
  rownames(IMSData) <-IMSData$Sample.ID
    ### rename IMS names
      #### rename pam50.prediction.subtype names
      levels(IMSData$pam50.prediction.subtype) <- c(levels(IMSData$pam50.prediction.subtype),levels(IMS.TCGA$PAM50)) 
      IMSData$pam50.prediction.subtype[IMSData$pam50.prediction.subtype=="Basal"]  <- "Basal-like"
      IMSData$pam50.prediction.subtype[IMSData$pam50.prediction.subtype=="Her2"]   <- "HER2-enriched"
      IMSData$pam50.prediction.subtype[IMSData$pam50.prediction.subtype=="LumA"]   <- "Luminal A"
      IMSData$pam50.prediction.subtype[IMSData$pam50.prediction.subtype=="LumB"]   <- "Luminal B"
      IMSData$pam50.prediction.subtype[IMSData$pam50.prediction.subtype=="Normal"] <- "Normal-like"
      IMSData$pam50.prediction.subtype <- droplevels(IMSData$pam50.prediction.subtype)
      #### rename Sorlie.prediction.subtype names
      levels(IMSData$Sorlie.prediction.subtype) <- c(levels(IMSData$Sorlie.prediction.subtype),levels(IMS.TCGA$PAM50)) 
      IMSData$Sorlie.prediction.subtype[IMSData$Sorlie.prediction.subtype=="Basal"]  <- "Basal-like"
      IMSData$Sorlie.prediction.subtype[IMSData$Sorlie.prediction.subtype=="Her2"]   <- "HER2-enriched"
      IMSData$Sorlie.prediction.subtype[IMSData$Sorlie.prediction.subtype=="LumA"]   <- "Luminal A"
      IMSData$Sorlie.prediction.subtype[IMSData$Sorlie.prediction.subtype=="LumB"]   <- "Luminal B"
      IMSData$Sorlie.prediction.subtype[IMSData$Sorlie.prediction.subtype=="Normal"] <- "Normal-like"
      IMSData$Sorlie.prediction.subtype <- droplevels(IMSData$Sorlie.prediction.subtype)
      #### rename Hu.prediction.subtype names
      levels(IMSData$Hu.prediction.subtype) <- c(levels(IMSData$Hu.prediction.subtype),levels(IMS.TCGA$PAM50)) 
      IMSData$Hu.prediction.subtype[IMSData$Hu.prediction.subtype=="Basal"]  <- "Basal-like"
      IMSData$Hu.prediction.subtype[IMSData$Hu.prediction.subtype=="Her2"]   <- "HER2-enriched"
      IMSData$Hu.prediction.subtype[IMSData$Hu.prediction.subtype=="LumA"]   <- "Luminal A"
      IMSData$Hu.prediction.subtype[IMSData$Hu.prediction.subtype=="LumB"]   <- "Luminal B"
      IMSData$Hu.prediction.subtype[IMSData$Hu.prediction.subtype=="Normal"] <- "Normal-like"
      IMSData$Hu.prediction.subtype <- droplevels(IMSData$Hu.prediction.subtype)
  print ("Data loaded...")

# subset clinical data
  ClinicalData.subset <- subset (clinicalData,is.element (clinicalData$bcr_patient_barcode,row.names(agilentData)))
  row.names(ClinicalData.subset) <- ClinicalData.subset$bcr_patient_barcode
  ClinicalData.subset$bcr_patient_barcode <- NULL
# merge clinical data with IMS data
  ClinicalData.subset <- merge(ClinicalData.subset, IMSData, by="row.names",all.x=TRUE,all.y=FALSE)
  row.names(ClinicalData.subset) <- ClinicalData.subset$Row.names
  ClinicalData.subset$Row.names <- NULL
  ClinicalData.subset <- merge(ClinicalData.subset,IMS.TCGA.METHOD["Call"],by="row.names",all.x=TRUE, all.y=FALSE)
  row.names(ClinicalData.subset) <- ClinicalData.subset$Row.names
  ClinicalData.subset$Row.names <- NULL
  ClinicalData.subset <- merge(ClinicalData.subset,IMS.TCGA["PAM50"],by="row.names",all.x=TRUE, all.y=FALSE)
  row.names(ClinicalData.subset) <- ClinicalData.subset$Row.names
  ClinicalData.subset$Row.names <- NULL
  
  colnames(ClinicalData.subset)[colnames(ClinicalData.subset) == "Call"] <- "TCGA.PAM50.RMethod"
  ClinicalData.MA.DATA <- merge(ClinicalData.subset, agilentData,by="row.names",all.x=TRUE, all.y=FALSE)
  print ("Data merged...")
 
# test match in PAM50 IMS
  ClinicalData.subset <- cbind(ClinicalData.subset,ClinicalData.subset$pam50.prediction.subtype == ClinicalData.subset$PAM50)
  col.names <- colnames(ClinicalData.subset)
  col.names[length(col.names)] <- "TCGA.CallvsGenefu.PAM50"
  colnames(ClinicalData.subset) <- col.names
  ClinicalData.subset <- cbind(ClinicalData.subset,ClinicalData.subset$TCGA.PAM50.RMethod == ClinicalData.subset$PAM50)
  col.names <- colnames(ClinicalData.subset)
  col.names[length(col.names)] <- "TCGA.CallvsTCGA.PAM50.RMethod"
  colnames(ClinicalData.subset) <- col.names
  
  print ("TCGA.Call")
  print (table (ClinicalData.subset$PAM50))
  print ("Genefu.PAM50")
  print (table (ClinicalData.subset$pam50.prediction.subtype))
  print (table (ClinicalData.subset$TCGA.CallvsGenefu.PAM50))
  print ("TCGA.Rmethod.PAM50")
  print (table (ClinicalData.subset$TCGA.PAM50.RMethod))
  print (table (ClinicalData.subset$TCGA.CallvsTCGA.PAM50.RMethod))
  nomatch <- ClinicalData.subset[which(ClinicalData.subset$TCGA.CallvsTCGA.PAM50.RMethod == FALSE),31:32]
  nomatch <- merge(nomatch,IMS.TCGA.METHOD[,1:7],by="row.names",all.x=TRUE, all.y=FALSE)
  #write.csv (nomatch, file = "./DATA/IMS/TCGA IMS/nomatchwithpub.csv",row.names = TRUE);
  
# append exclusion parameter
  exclude.samples.males <- which(ClinicalData.subset$gender == "MALE")
  exclude.samples.nat <- which(ClinicalData.subset$history_neoadjuvant_treatment == "Yes")
  exclude.samples.histo <- which(ClinicalData.subset$histological_type %nin% c("Infiltrating Ductal Carcinoma","Infiltrating Lobular Carcinoma","Mixed Histology (please specify)"))
  exclude.samples <- unique(c(exclude.samples.males,exclude.samples.nat,exclude.samples.histo ))
  ClinicalData.subset$exclude <- "No"
  ClinicalData.subset[exclude.samples,"exclude"] <-"Yes"
  print ("exclusion parameter added...")
  
# export data to txt and excell
  write.csv (ClinicalData.subset, file = "./3 ANALISYS/CLINICAL DATA/Agilent_subset_clinicaldata.csv",row.names = TRUE);
  write.xlsx (ClinicalData.subset, file = "./3 ANALISYS/CLINICAL DATA/Agilent_subset_clinicaldata.xlsx", sheetName ="Agilent MA subset clinical data", row.names=TRUE);
  write.csv (as.data.frame(ClinicalData.MA.DATA), file = "./3 ANALISYS/CLINICAL DATA/Agilent_subset_clinical_AND_MA_DATA.csv",row.names = TRUE)
  print ("Data on all Samples are saved in Agilent_subset_clinicaldata.xlsx and Agilent_subset_clinicaldata.txt.");
 
