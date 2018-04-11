#################################################################
###
### This script generates updated survival files (including vars "bcr_patient_barcode", "vital_status","days_to_last_followup","days_to_death")
### from original patient.txt files.
###
### Data is saved :
### .../3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/updatedsurvivaldata.csv"
###
### (Visual check of correct updating performed!)
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                  # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                        # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages <- c("RCurl","httr", "rjson", "stringr", "HGNChelper")
ipak(required.packages)

# Set Parameters
download.method = "TCGA_Assembler"
CancerTYPES     = "ALL"
Cancer_skip     = c("")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)

if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets <- length(CancerTYPES)

#Update the patient.txt files with follow-up data. When a patient appears twice in a follow-up table, the most recent observation is always the lowest row in the table with this
#patient barcode (see variable form completion date).

k = 3
for(k in 1:N.sets){
  Cancer = CancerTYPES[k]
  if(Cancer %in% Cancer_skip){next}
  Cancer_path = paste0 ("./2_Data/TCGA_Assembler/",Cancer,"/BiospecimenClinicalData/")
  
  ##Create vectors with file paths to patient.txt files and to followup txt files
  file.list.all <- list.files(Cancer_path, full.names = TRUE)
  patient.file <- file.list.all[grep("/patient.txt$",file.list.all)]
  file.list.followup <- file.list.all[grep("follow_up",file.list.all)]
  file.list.followup <- file.list.followup[grep("nte",file.list.followup, invert=TRUE)] #remove the new tumor event followup file (=different type of file with different columns)
  file.list.followup <- file.list.followup[order(file.list.followup)]
  
  if (length(file.list.followup)==0) { 
    assign(paste0(Cancer, ".patient.table"), patient.table)
    next}
  patient.table = read.csv(patient.file,header = FALSE, sep="\t", as.is=TRUE, skip = 3,stringsAsFactors = FALSE)
  colnames(patient.table) = read.csv(patient.file, sep="\t", nrows = 1, as.is=TRUE)
  
  #survival.table <- patient.table[, c("bcr_patient_barcode", "vital_status","days_to_last_followup","days_to_death")]
  N.files <- length(file.list.followup)
  
  for(l in 1:N.files){
    followup.file = file.list.followup[l]
    followup.table = read.csv(followup.file,header = FALSE, sep="\t", as.is=TRUE, skip = 3,stringsAsFactors = FALSE)
    colnames(followup.table) = read.csv(followup.file, sep="\t", nrows = 1, as.is=TRUE)
    
    rev.followup.table <- followup.table[nrow(followup.table):1,] #reverse order of rows to get most recent observation to first appear in the rows.
    rev.followup.table2 <- rev.followup.table[!duplicated(rev.followup.table$bcr_patient_barcode),c("bcr_patient_barcode", "vital_status","days_to_last_followup","days_to_death")] #remove rows that have the same patient_barcode as a row higher in dataframe.
    
    to.update <- which(patient.table$bcr_patient_barcode %in% rev.followup.table2$bcr_patient_barcode)
    updates <- match(patient.table$bcr_patient_barcode[to.update], rev.followup.table2$bcr_patient_barcode)
    patient.table[to.update,c("bcr_patient_barcode", "vital_status","days_to_last_followup","days_to_death")] <- rev.followup.table2[updates,]
  }
  
  assign(paste0(Cancer, ".variables"), colnames(patient.table))
  assign(paste0(Cancer, ".patient.table"), patient.table)
  #dir.create("./3_DataProcessing/",showWarnings = FALSE)
  #dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
  #dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer),showWarnings = FALSE)
  #dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/"),showWarnings = FALSE)
  #write.csv(survival.table, file = paste0("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/updatedsurvivaldata.csv"))
}

variable_vectors = grep(".variables", names(.GlobalEnv), value = TRUE)
variables_list = do.call("list", mget(variable_vectors))

variables_in_all_Cancers = Reduce(intersect, variables_list)
variables_table = as.data.frame(table(unname(unlist(variables_list))))
variables_table= variables_table[order(variables_table$Freq, decreasing = TRUE),]

