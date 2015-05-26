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

# Parameters
  Cancerset <- "READ-GA"
  Parent.Cancerset <- substring(Cancerset,1,4) 
  
# Load data files 
  ## RNASeq DAta from TCGA assembler
  load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
  PatientIDs <- unique(substr(colnames(RNASeq.NORM_Log2),1,12)) 
  rm(RNASeq.NORM_Log2)
  ## clinical data retrieved through TCGA assembler
  clinicalData <- read.csv (paste0("./2 DATA/Clinical Information/",Parent.Cancerset,"/selected_clinicaldata.txt"), header = TRUE)
    
# subset clinical data
  ClinicalData.subset <- subset (clinicalData,is.element (clinicalData$bcr_patient_barcode,PatientIDs))
  row.names(ClinicalData.subset) <- ClinicalData.subset$bcr_patient_barcode
  ClinicalData.subset$bcr_patient_barcode <- NULL

# append exclusion parameter
  exclude.samples.nat <- which(ClinicalData.subset$history_neoadjuvant_treatment == "Yes")
  #exclude.samples.histo <- which(ClinicalData.subset$histological_type %nin% c("Infiltrating Ductal Carcinoma","Infiltrating Lobular Carcinoma","Mixed Histology (please specify)"))
  exclude.samples <- unique(c(exclude.samples.nat))
  ClinicalData.subset$exclude <- "No"
  ClinicalData.subset[exclude.samples,"exclude"] <-"Yes"
  print ("exclusion parameter added...")
  
# export data to txt and excel
  write.csv (ClinicalData.subset, file = paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"),row.names = TRUE);
  write.xlsx (ClinicalData.subset, file = paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.xlsx"), sheetName ="RNASeq subset clinical data", row.names=TRUE);
  print (paste0("Data on all Samples are saved in TCGA.",Cancerset,".RNASeq_subset_clinicaldata.xlsx and TCGA.",Cancerset,".RNASeq_subset_clinicaldata.txt."))
