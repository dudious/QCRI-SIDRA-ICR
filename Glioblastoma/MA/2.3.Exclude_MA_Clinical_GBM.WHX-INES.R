#################################################################
###
### This Script adds the exlusion parameter to clinical data.
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR/")
  ## dependencies
  ## install java for xlsx export
  ## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
  required.packages <- c("xlsx","Hmisc","HGNChelper")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  library (xlsx) #xlsx needs java installed
  library (Hmisc)
 
  source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
  source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")
  
# Load data files 
  ## RNASeq DAta from TCGA assembler
  load ("./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.RData")
  AffyData.NORM<-t(AffyData.NORM) 
  dim(AffyData.NORM)#12042 x 529
  PatientIDs <- unique(substr(colnames(AffyData.NORM),1,12)) 
  
  ## clinical data retrieved through TCGA assembler
  clinicalData <- read.csv ("./2 DATA/Clinical Information/GBM/selected_clinicaldata.txt", header = TRUE)
    
# subset clinical data
# PatientIDs:83
# clinicalData:594

  ClinicalData.subset <- subset (clinicalData,is.element (clinicalData$bcr_patient_barcode,PatientIDs))#number:525
  row.names(ClinicalData.subset) <- ClinicalData.subset$bcr_patient_barcode
  ClinicalData.subset$bcr_patient_barcode <- NULL

# append exclusion parameter
  exclude.samples.nat <- which(ClinicalData.subset$history_neoadjuvant_treatment == "Yes")#21
  #exclude.samples.histo <- which(ClinicalData.subset$histological_type %nin% c("Infiltrating Ductal Carcinoma","Infiltrating Lobular Carcinoma","Mixed Histology (please specify)"))
  exclude.samples <- unique(c(exclude.samples.nat))
  ClinicalData.subset$exclude <- "No"
  ClinicalData.subset[exclude.samples,"exclude"] <-"Yes"
  print ("exclusion parameter added...")
  
# export data to txt and excel
  write.csv (ClinicalData.subset, file = "./3 ANALISYS/CLINICAL DATA/TCGA.GBM.MA.Affy_subset_clinicaldata.csv",row.names = TRUE);
  write.xlsx (ClinicalData.subset, file = "./3 ANALISYS/CLINICAL DATA/TCGA.GBM.MA.Affy_subset_clinicaldata.xlsx", sheetName ="RNASeq subset clinical data", row.names=TRUE);
  print ("Data on all Samples are saved in RNASeq_subset_clinicaldata.xlsx and RNASeq_subset_clinicaldata.txt.");
 
