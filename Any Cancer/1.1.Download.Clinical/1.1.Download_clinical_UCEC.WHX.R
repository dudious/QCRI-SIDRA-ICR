#################################################################
###
### This script downloads the Liver cancer Clinical Data and
### Biospecimen Data from the TCGA database.
### It will proccess a subset of the data into a  single table. 
### Data is saved :
### ..\2 DATA\Clinical Information\UCEC\...
### File to use :
### "selected_clinicaldata.txt"
###
#################################################################

# Setup environment
## install java for xlsx export
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
required.packages <- c("xlsx","RCurl","httr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)

source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
library (xlsx)  #xlsx needs java installed

# Download de-identified clinical information of UCEC patients in the Biotab format
DownloadClinicalData(traverseResultFile = "./2 DATA/DirectoryTraverseResult_May-06-2015.rda", 
                     saveFolderName = "./2 DATA/Clinical Information/UCEC/RawData",
                     cancerType = "UCEC",
                     clinicalDataType = c("patient",                                          
                                          "cqcf",
                                          "drug",
                                          "radiation",
                                          "nte",                                          
                                          "follow_up"));

DownloadBiospecimenData(traverseResultFile = "./2 DATA/DirectoryTraverseResult_May-06-2015.rda",
                        saveFolderName = "./2 DATA/Biospecimeninfo/UCEC/",
                        cancerType = "UCEC",
                        biospecimenDataType = c("normal_control", "tumor_sample"));

# Extract and merge relevant data
## rename files
  file.list.old  <- list.files("./2 DATA/Clinical Information/UCEC/RawData",full.names = TRUE)
  file.list.new <- gsub ("nationwidechildrens.org_clinical_","",file.list.old)
  file.list.new <- gsub ("_ucec","",file.list.new)
  file.rename (file.list.old,file.list.new)
print ("Files renamed ...")


## load data
  cqcf.table <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/cqcf.txt", header = TRUE, sep="\t", as.is=TRUE)
  cqcf.table <- (cqcf.table [-c(1,2),]) # delete first 2 rows
  row.names(cqcf.table) <- NULL

  drug.table <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/drug.txt", header = TRUE, sep="\t", as.is=TRUE)
  drug.table <- (drug.table [-c(1,2),]) # delete first 2 rows
  row.names(drug.table) <- NULL

  nte.table <- read.csv("./2 DATA/Clinical Information/UCEC/RawData/nte.txt", header = TRUE, sep="\t", as.is=TRUE)
  nte.table <- (nte.table [-c(1,2),]) # delete first 2 rows
  row.names(nte.table) <- NULL

  patient.table <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/patient.txt", header = TRUE, sep="\t", as.is=TRUE)
  patient.table <- (patient.table [-c(1,2),]) # delete first 2 rows
  row.names(patient.table) <- NULL

  radiation.table <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/radiation.txt", header = TRUE, sep="\t", as.is=TRUE)
  radiation.table <- (radiation.table [-c(1,2),]) # delete first 2 rows
  row.names(radiation.table) <- NULL

  follow_up.table.1.7 <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/follow_up_v1.7.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.1.7 <- (follow_up.table.1.7 [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.1.7) <- NULL

  follow_up.table.2.0 <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/follow_up_v2.0.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.2.0 <- (follow_up.table.2.0 [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.2.0) <- NULL

  follow_up.table.4.0 <- read.csv ("./2 DATA/Clinical Information/UCEC/RawData/follow_up_v4.0.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.4.0 <- (follow_up.table.4.0 [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.4.0) <- NULL

  follow_up.table.4.0_nte <- read.csv ("./2 DATA/Clinical Information/LIHC/RawData/follow_up_v4.0_nte.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.4.0_nte <- (follow_up.table.4.0_nte [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.4.0_nte) <- NULL


## selection of variables
#paste(colnames (patient.table),collapse=",")
  patient.vars <- c("bcr_patient_barcode","gender","race","ethnicity","history_other_malignancy","history_neoadjuvant_treatment",
                 "tumor_status","vital_status",
                 "pathologic_T","pathologic_N","pathologic_M","pathologic_stage",
                 "new_tumor_event_dx_indicator",
                 "age_at_diagnosis","birth_days_to","clinical_M","clinical_N","clinical_T","clinical_stage","death_days_to",
                 "last_contact_days_to")
  cqcf.vars <- c("bcr_patient_barcode")
  drug.vars <- c("bcr_patient_barcode") 
  radiation.vars <- c("bcr_patient_barcode")
  nte.vars <- c("bcr_patient_barcode")
  follow_up.vars <- c("bcr_patient_barcode","tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")

## build table of selected clinical data
patient.vars[-which(patient.vars %in% colnames(patient.table))] # check for unavailable variables
  # Patient + cqcf [dim : ]
  clinicaldata.table <- merge(patient.table[patient.vars],
                            cqcf.table[cqcf.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add drug [dim : ]
  clinicaldata.table <- merge(clinicaldata.table,
                            drug.table[drug.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add radiation [dim : ]
  clinicaldata.table <- merge(clinicaldata.table,
                            radiation.table[radiation.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add nte (new tumour event)  [dim : ]
  clinicaldata.table <- merge(clinicaldata.table,
                            nte.table[nte.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add Follow up
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.1.7[c("bcr_patient_barcode","tumor_status","vital_status","last_contact_days_to","death_days_to")], #no time for NTE "new_tumor_event_dx_days_to"
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  colnames(clinicaldata.table)[(ncol(clinicaldata.table)-3):ncol(clinicaldata.table)] <- 
                              paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.7")
  clinicaldata.table <- unique(clinicaldata.table)

  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.2.0[c("bcr_patient_barcode","tumor_status","vital_status","last_contact_days_to","death_days_to")], #no time for NTE "new_tumor_event_dx_days_to"
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  colnames(clinicaldata.table)[(ncol(clinicaldata.table)-3):ncol(clinicaldata.table)] <- 
                              paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_2.0")
  clinicaldata.table <- unique(clinicaldata.table)

  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.4.0[c("bcr_patient_barcode","tumor_status","vital_status","last_contact_days_to","death_days_to")], #no time for NTE "new_tumor_event_dx_days_to"
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  colnames(clinicaldata.table)[(ncol(clinicaldata.table)-3):ncol(clinicaldata.table)] <- 
                              paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_4.0")
  clinicaldata.table <- unique(clinicaldata.table)
   
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.nte[c("bcr_patient_barcode")], # this is histo data on nte tumours
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  
  # select latest follow up data
  # original data from patient file = 1.0
  colnames (clinicaldata.table)[c(which(colnames(clinicaldata.table) %in% 
                                c("tumor_status.x","vital_status.x","last_contact_days_to.x","death_days_to.x")))] <-
                                paste0(c("tumor_status","vital_status","death_days_to","last_contact_days_to"),"_1.0")
  # make all time columns numerical
  clinicaldata.table[c(paste0("last_contact_days_to",c("_1.0","_1.7","_2.0","_4.0")))] <- as.numeric(as.character(unlist(clinicaldata.table[c(paste0("last_contact_days_to",c("_1.0","_1.7","_2.0","_4.0")))])))
  clinicaldata.table[c(paste0("death_days_to",c("_1.0","_1.7","_2.0","_4.0")))] <- as.numeric(as.character(unlist(clinicaldata.table[c(paste0("death_days_to",c("_1.0","_1.7","_2.0","_4.0")))])))  

  # latest data = 1.0/1.5
  clinicaldata.table <- cbind (clinicaldata.table,clinicaldata.table[c(paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.0"))])
  colnames (clinicaldata.table)[(ncol(clinicaldata.table)-3):ncol(clinicaldata.table)] <-
                              c("tumor_status","vital_status","last_contact_days_to","death_days_to")
  # update to 1.7 ( if last_contact_days_to_1.7 > last_contact_days_to OR death_days_to_1.7 > death_days_to)
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_1.7 > clinicaldata.table$last_contact_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to")]<-
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_1.7 > clinicaldata.table$last_contact_days_to)),
                     paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.7")]
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_1.7 > clinicaldata.table$death_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to")]<-
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_1.7 > clinicaldata.table$death_days_to)),
                     paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.7")]

  # update to 2.0 ( if last_contact_days_to_2.0 > last_contact_days_to OR death_days_to_2.0 > death_days_to)
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_2.0 > clinicaldata.table$last_contact_days_to)),
                     c("tumor_status","vital_status","last_contact_days_to","death_days_to")]<-
                      clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_2.0 > clinicaldata.table$last_contact_days_to)),
                      paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_2.0")]
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_2.0 > clinicaldata.table$death_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to")]<-
                      clinicaldata.table[c(which (clinicaldata.table$death_days_to_2.0 > clinicaldata.table$death_days_to)),
                      paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_2.0")]

  # update to 4.0 ( if last_contact_days_to_4.0 > last_contact_days_to OR death_days_to_4.0 > death_days_to)
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_4.0 > clinicaldata.table$last_contact_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to")]<-
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_4.0 > clinicaldata.table$last_contact_days_to)),
                     paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_4.0")]
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_4.0 > clinicaldata.table$death_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to")]<-
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_4.0 > clinicaldata.table$death_days_to)),
                     paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_4.0")]

  # remove redundant follow up data
  clinicaldata.table <- clinicaldata.table[-c(which(colnames(clinicaldata.table) %in% 
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.0")))]
  clinicaldata.table <- clinicaldata.table[-c(which(colnames(clinicaldata.table) %in% 
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.7")))]  
  clinicaldata.table <- clinicaldata.table[-c(which(colnames(clinicaldata.table) %in% 
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_2.0")))]
  clinicaldata.table <- clinicaldata.table[-c(which(colnames(clinicaldata.table) %in% 
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_4.0")))]

  #cleanup
  # dim : 566 x 21 (26 doubles)
  clinicaldata.table <- unique(clinicaldata.table)
  # dim : 565 x 21 (25 doubles)
  clinicaldata.table$last_contact_days_to <- as.numeric(as.character(clinicaldata.table$last_contact_days_to))
  clinicaldata.table$death_days_to <- as.numeric(as.character(clinicaldata.table$death_days_to))
  clinicaldata.table <- clinicaldata.table[order(clinicaldata.table$bcr_patient_barcode,
                                                 -abs(clinicaldata.table$last_contact_days_to)), ]    # order follow up so the longest FU time for duplicated patients comes firts
  clinicaldata.table <- clinicaldata.table [ !duplicated(clinicaldata.table$bcr_patient_barcode), ]	  # remove duplicates with identical or shorter FU time	
  # dim : 540x 21 ( NO doubles)
  row.names(clinicaldata.table) <- NULL
  print ("Data Merged ...")
  print ("Selected data : ")
  print (colnames(clinicaldata.table))

# export data to txt and excell
  write.csv (clinicaldata.table, file = "./2 DATA/Clinical Information/UCEC/selected_clinicaldata.txt", row.names=FALSE);
  write.xlsx (clinicaldata.table, file = "./2 DATA/Clinical Information/UCEC/selected_clinicaldata.xlsx", sheetName ="selected clinical data", row.names=FALSE);
  print ("Results are saved in selected_clinicaldata.xlsx and selected_clinicaldata.txt.");


