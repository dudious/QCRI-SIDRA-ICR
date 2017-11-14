#################################################################
###
### This Script adds the exlusion parameters to clinical data.
###
#################################################################

# Setup environment
  rm(list=ls())
  #setwd("~/Dropbox/BREAST_QATAR/")
  setwd("/mnt3/wouter/BREAST-QATAR/")
  ## dependencies
  ## install java for xlsx export ( in ubuntu : sudo apt-get install openjdk-7-jdk and sudo R CMD javareconf)
  ## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
  ## install libcurl on ubuntu if needed (sudo apt-get install libcurl4-gnutls-dev)
  required.packages <- c("xlsx","Hmisc","HGNChelper","RCurl","httr")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  #library (xlsx) #xlsx needs java installed
  library (Hmisc)
 
# Parameters
  Cancersets   = "ALL"
  DL.Method    = "BIOLINKS" #Choose "ASSEMBLER" or "BIOLINKS"
  sample.types = "Selected" #Alternatives TP , TP_TM , Selected
  
  
# DO ALL
  TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
  if (Cancersets == "ALL") { 
    Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
  }
  N.sets = length(Cancersets)
  for (i in 1:N.sets) {
    Cancerset = Cancersets[i]
    if (Cancerset %in% c("LAML","FPPP","BRCA")) {next}
    Parent.Cancerset <- substring(Cancerset,1,4)
  
if (DL.Method=="ASSEMBLER"){
  # Load data files 
   ## RNASeq DAta from TCGA assembler
   load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".NORMALIZED.LOG2.RData"))
   PatientIDs <- unique(substr(colnames(RNASeq.NORM_Log2),1,12)) 
   rm(RNASeq.NORM_Log2)
   ## clinical data retrieved through TCGA assembler
   clinicalData <- read.csv (paste0("./2 DATA/Clinical Information/",Parent.Cancerset,"/selected_clinicaldata.txt"), header = TRUE)
  # subset clinical data
  ClinicalData.subset <- subset (clinicalData,is.element (clinicalData$bcr_patient_barcode,PatientIDs))
  row.names(ClinicalData.subset) <- ClinicalData.subset$bcr_patient_barcode
  ClinicalData.subset$bcr_patient_barcode <- NULL
  exclude.samples.nat <- which(ClinicalData.subset$history_neoadjuvant_treatment == "Yes")
  exclude.samples.history <- which(ClinicalData.subset$history_other_malignancy == "Yes")
}

if (DL.Method=="BIOLINKS"){
  # Load data files 
   ## RNASeq DAta from TCGA Biolinks
  load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".",sample.types,".NORMALIZED.TP_FILTERED_LOG2.RData"))
  PatientIDs <- colnames(RNASeq.NORM.TP_Log2)
  rm(RNASeq.NORM.TP_Log2)
  ## indexed clinical data retrieved through TCGA BIOLINKS
  clinicalData.indexed <- read.csv (paste0("./2 DATA/Clinical Information/",Parent.Cancerset,"/BIOLINKS/indexed.GDC.clinicaldata.csv"), header = TRUE)
  clinicalData.XML <- read.csv (paste0("./2 DATA/Clinical Information/",Parent.Cancerset,"/BIOLINKS/Patient.XML.clinicaldata.csv"), header = TRUE)
  row.names(clinicalData.indexed) <- clinicalData.indexed$bcr_patient_barcode
  row.names(clinicalData.XML) <- clinicalData.XML$bcr_patient_barcode
  clinicalData.indexed.subset <- clinicalData.indexed[PatientIDs,]
  clinicalData.XML.subset <- clinicalData.XML[PatientIDs,]
  exclude.samples.nat <- rownames(clinicalData.XML.subset[which(clinicalData.XML.subset$history_of_neoadjuvant_treatment == "Yes"),])
  exclude.samples.history <- rownames(clinicalData.indexed.subset[which(clinicalData.indexed.subset$prior_malignancy != "not reported"),])
  ## Select Data
  ClinicalData.subset <- clinicalData.indexed.subset[,c("gender","vital_status","days_to_death","days_to_last_follow_up","prior_malignancy")]
  ClinicalData.subset$history_neoadjuvant_treatment <- clinicalData.XML.subset$history_of_neoadjuvant_treatment[match(rownames(ClinicalData.subset),rownames(clinicalData.XML.subset))]
  ## ASSEMBLER PROCESSED FILE FOR COMPARISON
  #clinicalData.ASSEMBLER.subset <- read.csv(paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))
  #rownames(clinicalData.ASSEMBLER.subset) <- clinicalData.ASSEMBLER.subset$X
  #clinicalData.ASSEMBLER.subset <- clinicalData.ASSEMBLER.subset[rownames(ClinicalData.subset),]
  #exclude.samples.nat <- unique(c(exclude.samples.nat,rownames(clinicalData.ASSEMBLER.subset[which(clinicalData.ASSEMBLER.subset$history_neoadjuvant_treatment == "Yes"),])))
  #exclude.samples.history <- unique(c(exclude.samples.history,rownames(clinicalData.ASSEMBLER.subset[which(clinicalData.ASSEMBLER.subset$history_other_malignancy == "Yes"),])))
}
  
# append exclusion parameters
  exclude.samples.preclust <- unique(exclude.samples.nat) 
  exclude.samples.postclust <- unique(exclude.samples.history)
  ClinicalData.subset$exclude.pre <- "No"
  ClinicalData.subset$exclude.post <- "No"
  ClinicalData.subset[exclude.samples.preclust,"exclude.pre"] <-"Yes"
  ClinicalData.subset[exclude.samples.postclust,"exclude.post"] <-"Yes"
  print ("exclusion parameter added...")
  
# export data to txt and excel
  write.csv (ClinicalData.subset, file = paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_",DL.Method,"_subset_clinicaldata.csv"),row.names = TRUE);
 # write.xlsx (ClinicalData.subset, file = paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.xlsx"), sheetName ="RNASeq subset clinical data", row.names=TRUE);
  print (paste0("Data on all Samples are saved in TCGA.",Cancerset,".RNASeq_",DL.Method,"_subset_clinicaldata.txt."))
  
  }
