#################################################################
###
### This script downloads the Breast cancer Clinical Data and
### Biospecimen Data from the TCGA database.
### It will proccess a subset of the data into a  single table. 
### Data is saved :
### ..\2 DATA\Clinical Information\BRCA\...
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

# Download de-identified clinical information of BRCA patients in the Biotab format
DownloadClinicalData(traverseResultFile = "./2 DATA/DirectoryTraverseResult_Jan-07-2015.rda", 
                     saveFolderName = "./2 DATA/Clinical Information/BRCA/RawData",
                     cancerType = "BRCA",
                     clinicalDataType = c("patient",                                          
                                          "cqcf",
                                          "drug",
                                          "radiation",
                                          "nte",                                          
                                          "follow_up"));

DownloadBiospecimenData(traverseResultFile = "./2 DATA/DirectoryTraverseResult_Jan-07-2015.rda",
                        saveFolderName = "./2 DATA/Biospecimeninfo/BRCA/",
                        cancerType = "BRCA",
                        biospecimenDataType = c("normal_control", "tumor_sample"));

# Extract and merge relevant data
## rename files
  file.list.old  <- list.files("./2 DATA/Clinical Information/BRCA/RawData",full.names = TRUE)
  file.list.new <- gsub ("nationwidechildrens.org_clinical_","",file.list.old)
  file.list.new <- gsub ("_brca","",file.list.new)
  file.rename (file.list.old,file.list.new)
print ("Files renamed ...")


## load data
  cqcf.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/cqcf.txt", header = TRUE, sep="\t", as.is=TRUE)
  cqcf.table <- (cqcf.table [-c(1,2),]) # delete first 2 rows
  row.names(cqcf.table) <- NULL

  drug.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/drug.txt", header = TRUE, sep="\t", as.is=TRUE)
  drug.table <- (drug.table [-c(1,2),]) # delete first 2 rows
  row.names(drug.table) <- NULL

  nte.table <- read.csv("./2 DATA/Clinical Information/BRCA/RawData/nte.txt", header = TRUE, sep="\t", as.is=TRUE)
  nte.table <- (nte.table [-c(1,2),]) # delete first 2 rows
  row.names(nte.table) <- NULL

  patient.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/patient.txt", header = TRUE, sep="\t", as.is=TRUE)
  patient.table <- (patient.table [-c(1,2),]) # delete first 2 rows
  row.names(patient.table) <- NULL

  radiation.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/radiation.txt", header = TRUE, sep="\t", as.is=TRUE)
  radiation.table <- (radiation.table [-c(1,2),]) # delete first 2 rows
  row.names(radiation.table) <- NULL

  follow_up.table.1.5 <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/follow_up_v1.5.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.1.5 <- (follow_up.table.1.5 [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.1.5) <- NULL

  follow_up.table.2.1 <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/follow_up_v2.1.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.2.1 <- (follow_up.table.2.1 [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.2.1) <- NULL

  follow_up.table.4.0 <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/follow_up_v4.0.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.4.0 <- (follow_up.table.4.0 [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.4.0) <- NULL

  follow_up.table.4.0_nte <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/follow_up_v4.0_nte.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table.4.0_nte <- (follow_up.table.4.0_nte [-c(1,2),]) # delete first 2 rows
  row.names(follow_up.table.4.0_nte) <- NULL



## selection of variables
  patient.vars <- c("bcr_patient_barcode","birth_days_to","gender","menopause_status","race","history_other_malignancy",
                    "history_neoadjuvant_treatment","tumor_status","vital_status","last_contact_days_to","death_days_to",
                    "initial_pathologic_dx_year","age_at_diagnosis","histological_type","lymph_nodes_examined_he_count","ajcc_tumor_pathologic_pt",
                    "ajcc_nodes_pathologic_pn","ajcc_metastasis_pathologic_pm","ajcc_pathologic_tumor_stage","metastasis_site",
                    "er_status_by_ihc","pr_status_by_ihc","her2_status_by_ihc","her2_fish_status","her2_copy_number",
                    "cent17_copy_number","her2_cent17_ratio")
  cqcf.vars <- c("bcr_patient_barcode")
  drug.vars <- c("bcr_patient_barcode") 
  radiation.vars <- c("bcr_patient_barcode")
  nte.vars <- c("bcr_patient_barcode")
  follow_up.vars <- c("bcr_patient_barcode","tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")

## build table of selected clinical data
  # Patient + cqcf [dim : 1070  x 27]
  clinicaldata.table <- merge(patient.table[patient.vars],
                            cqcf.table[cqcf.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add drug [dim : 1070  x 27]
  clinicaldata.table <- merge(clinicaldata.table,
                            drug.table[drug.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add radiation [dim : 1070  x 27]
  clinicaldata.table <- merge(clinicaldata.table,
                            radiation.table[radiation.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add nte (new tumour event)  [dim : 1070  x 27]
  clinicaldata.table <- merge(clinicaldata.table,
                            nte.table[nte.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  # add Follow up 
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.1.5[follow_up.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  colnames(clinicaldata.table)[(ncol(clinicaldata.table)-4):ncol(clinicaldata.table)] <- 
                              paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_1.5")
  clinicaldata.table <- unique(clinicaldata.table)
  # [dim : 1071 x 32 ] -> "TCGA-A2-A04P" 1 double 
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.2.1[follow_up.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  colnames(clinicaldata.table)[(ncol(clinicaldata.table)-4):ncol(clinicaldata.table)] <- 
                              paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_2.1")
  clinicaldata.table <- unique(clinicaldata.table)
  # [dim : 1076 x 37 ] -> 5 extra doubles 
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.4.0[c("bcr_patient_barcode","tumor_status","vital_status","last_contact_days_to","death_days_to")], #no time for NTE "new_tumor_event_dx_days_to"
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  colnames(clinicaldata.table)[(ncol(clinicaldata.table)-3):ncol(clinicaldata.table)] <- 
                              paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_4.0")
  clinicaldata.table <- unique(clinicaldata.table)
  # [dim : 1088 x 41 ] -> 12 extra doubles 
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table.4.0_nte[c("bcr_patient_barcode")], # this is histo data on nte tumours
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  
  # select latest follow up data
  # original data from patient file = 1.0
  colnames (clinicaldata.table)[c(which(colnames(clinicaldata.table) %in% 
                                c("tumor_status.x","vital_status.x","last_contact_days_to.x","death_days_to.x")))] <-
                                paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.0")
  # make all time columns numerical
  clinicaldata.table[c(paste0("last_contact_days_to",c("_1.0","_1.5","_2.1","_4.0")))] <- as.numeric(as.character(unlist(clinicaldata.table[c(paste0("last_contact_days_to",c("_1.0","_1.5","_2.1","_4.0")))])))
  
  # latest data = 1.0/1.5
  clinicaldata.table <- cbind (clinicaldata.table,clinicaldata.table[c(paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_1.0"),"new_tumor_event_dx_days_to_1.5")])
  colnames (clinicaldata.table)[(ncol(clinicaldata.table)-4):ncol(clinicaldata.table)] <-
                              c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")
  # update to 1.5 ( if last_contact_days_to_1.5 > last_contact_days_to OR death_days_to_1.5 > death_days_to)
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_1.5 > clinicaldata.table$last_contact_days_to)),
                     c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")]<-
                      clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_1.5 > clinicaldata.table$last_contact_days_to)),
                      paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_1.5")]
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_1.5 > clinicaldata.table$death_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")]<-
                      clinicaldata.table[c(which (clinicaldata.table$death_days_to_1.5 > clinicaldata.table$death_days_to)),
                      paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_1.5")]
  # update to 2.1 ( if last_contact_days_to_2.1 > last_contact_days_to OR death_days_to_2.1 > death_days_to)
  clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_2.1 > clinicaldata.table$last_contact_days_to)),
                     c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")]<-
                      clinicaldata.table[c(which (clinicaldata.table$last_contact_days_to_2.1 > clinicaldata.table$last_contact_days_to)),
                      paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_2.1")]
  clinicaldata.table[c(which (clinicaldata.table$death_days_to_2.1 > clinicaldata.table$death_days_to)),
                   c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to")]<-
                      clinicaldata.table[c(which (clinicaldata.table$death_days_to_2.1 > clinicaldata.table$death_days_to)),
                      paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_2.1")]
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
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_1.5")))]
  clinicaldata.table <- clinicaldata.table[-c(which(colnames(clinicaldata.table) %in% 
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to"),"_2.1")))]
  clinicaldata.table <- clinicaldata.table[-c(which(colnames(clinicaldata.table) %in% 
                        paste0(c("tumor_status","vital_status","last_contact_days_to","death_days_to"),"_4.0")))]

  #cleanup
  # dim : 1088 x 27 (18 doubles)
  clinicaldata.table <- unique(clinicaldata.table)
  # dim : 1083 x 27 (13 doubles)
  clinicaldata.table$last_contact_days_to <- as.numeric(as.character(clinicaldata.table$last_contact_days_to))
  clinicaldata.table$death_days_to <- as.numeric(as.character(clinicaldata.table$death_days_to))
  clinicaldata.table <- clinicaldata.table[order(clinicaldata.table$bcr_patient_barcode,
                                                 -abs(clinicaldata.table$last_contact_days_to)), ]    # order follow up so the longest FU time for duplicated patients comes firts
  clinicaldata.table <- clinicaldata.table [ !duplicated(clinicaldata.table$bcr_patient_barcode), ]	  # remove duplicates with identical or shorter FU time	
  # dim : 1070 x 27 ( NO doubles)
  row.names(clinicaldata.table) <- NULL
  print ("Data Merged ...")
  print ("Selected data : ")
  print (colnames(clinicaldata.table))

# export data to txt and excell
  write.csv (clinicaldata.table, file = "./2 DATA/Clinical Information/BRCA/selected_clinicaldata.txt", row.names=FALSE);
  write.xlsx (clinicaldata.table, file = "./2 DATA/Clinical Information/BRCA/selected_clinicaldata.xlsx", sheetName ="selected clinical data", row.names=FALSE);
  print ("Results are saved in selected_clinicaldata.xlsx and selected_clinicaldata.txt.");


