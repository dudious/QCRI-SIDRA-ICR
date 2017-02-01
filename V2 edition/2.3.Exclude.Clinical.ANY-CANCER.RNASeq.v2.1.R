#################################################################
###
### This Script adds the exlusion parameters to clinical data.
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR/")
  ## dependencies
  ## install java for xlsx export ( in ubuntu : sudo apt-get install openjdk-7-jdk and sudo R CMD javareconf)
  ## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
  ## install libcurl on ubuntu if needed (sudo apt-get install libcurl4-gnutls-dev)
  required.packages <- c("xlsx","Hmisc","HGNChelper","RCurl","httr")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  library (xlsx) #xlsx needs java installed
  library (Hmisc)
 
  source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/Module_A.r")
  source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/Module_B.r")

# Parameters
  Cancerset    = "BLCA"
  DL.Method    = "ASSEMBLER" #Choose "ASSEMBLER" or "BIOLINKS"
  sample.types = "Selected" #Alternatives TP , TP_TM , Selected
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

# append exclusion parameters
  exclude.samples.nat <- which(ClinicalData.subset$history_neoadjuvant_treatment == "Yes")
  exclude.samples.history <- which(ClinicalData.subset$history_other_malignancy == "Yes")
  exclude.samples.preclust <- unique(exclude.samples.nat) 
  exclude.samples.postclust <- unique(exclude.samples.history)
  ClinicalData.subset$exclude.pre <- "No"
  ClinicalData.subset$exclude.post <- "No"
  ClinicalData.subset[exclude.samples.preclust,"exclude.pre"] <-"Yes"
  ClinicalData.subset[exclude.samples.postclust,"exclude.post"] <-"Yes"
  print ("exclusion parameter added...")
  
# export data to txt and excel
  write.csv (ClinicalData.subset, file = paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"),row.names = TRUE);
  write.xlsx (ClinicalData.subset, file = paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.xlsx"), sheetName ="RNASeq subset clinical data", row.names=TRUE);
  print (paste0("Data on all Samples are saved in TCGA.",Cancerset,".RNASeq_subset_clinicaldata.xlsx and TCGA.",Cancerset,".RNASeq_subset_clinicaldata.txt."))
