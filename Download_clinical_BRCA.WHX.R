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
## install.packages("xlsx")
## download TCGA assembler scripts
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
required.packages <- c("xlsx")
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
  cqcf.table <- (cqcf.table [-c(1,2),])
  row.names(cqcf.table) <- NULL

  drug.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/drug.txt", header = TRUE, sep="\t", as.is=TRUE)
  drug.table <- (drug.table [-c(1,2),])
  row.names(drug.table) <- NULL

  nte.table <- read.csv("./2 DATA/Clinical Information/BRCA/RawData/nte.txt", header = TRUE, sep="\t", as.is=TRUE)
  nte.table <- (nte.table [-c(1,2),])
  row.names(nte.table) <- NULL

  patient.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/patient.txt", header = TRUE, sep="\t", as.is=TRUE)
  patient.table <- (patient.table [-c(1,2),])
  row.names(patient.table) <- NULL

  radiation.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/radiation.txt", header = TRUE, sep="\t", as.is=TRUE)
  radiation.table <- (radiation.table [-c(1,2),])
  row.names(radiation.table) <- NULL

  follow_up.table <- read.csv ("./2 DATA/Clinical Information/BRCA/RawData/follow_up_v4.0.txt", header = TRUE, sep="\t", as.is=TRUE)
  follow_up.table <- (follow_up.table [-c(1,2),])
  row.names(follow_up.table) <- NULL

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
  follow_up.vars <- c("bcr_patient_barcode")

## build table of selected clinical data
  clinicaldata.table <- merge(patient.table[patient.vars],
                            cqcf.table[cqcf.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  clinicaldata.table <- merge(clinicaldata.table,
                            drug.table[drug.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  clinicaldata.table <- merge(clinicaldata.table,
                            radiation.table[radiation.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  clinicaldata.table <- merge(clinicaldata.table,
                            nte.table[nte.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  clinicaldata.table <- merge(clinicaldata.table,
                            follow_up.table[follow_up.vars],
                            by.x="bcr_patient_barcode", by.y="bcr_patient_barcode",all.x=TRUE,all.y=FALSE);
  clinicaldata.table <- unique(clinicaldata.table)
  row.names(clinicaldata.table) <- NULL
  print ("Data Merged ...")
  print ("Selected data : ")
  print (colnames(clinicaldata.table))

# export data to txt and excell
  write.csv (clinicaldata.table, file = "./2 DATA/Clinical Information/BRCA/selected_clinicaldata.txt", row.names=FALSE);
  write.xlsx (clinicaldata.table, file = "./2 DATA/Clinical Information/BRCA/selected_clinicaldata.xlsx", sheetName ="selected clinical data", row.names=FALSE);
  print ("Results are saved in selected_clinicaldata.xlsx and selected_clinicaldata.txt.");

