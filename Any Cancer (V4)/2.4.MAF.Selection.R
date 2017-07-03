#################################################################
###
### This script generates an overview of the available MAFs 
### for somatic mutation data. Based on this information, a 
### selection of MAFs to analyze can be made.
### (Performed between script 1.6 and 1.7)
###
### Output file:
### "./1_Log_Files/2.4_MAF.selection/2.4.MAF.file.Analysis.", gsub(":",".",gsub(" ","_",date())), ".csv"
#################################################################

##Download the SomaticMutationData data using TCGA Assembler

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                  # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                        # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_A.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_B.R"))

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper", "Hmisc")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                     # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
Log_file = paste0("./1_Log_Files/1.6_SomaticMutation_Download/SomaticMutation_Download_Log_File_",                      # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
                                                                                                                        # in the Manual of Assembler v2.0.3 and was saved as csv file.
  
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

#Create overview of all downloaded Somatic Mutation files

N.sets = length(CancerTYPES)
SomaticMutationFiles = data.frame(Cancertype = CancerTYPES, Number.MAFs = 0, Number.Curated.MAFs = 0, Number.pairs.aggregated = 0, Centers = NA,
                                  Sample.count = NA, Patient.count = NA)
  
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/SomaticMutationData")
  file_list = list.files(Cancer_path)
  
  all_files = paste0(Cancer_path, "/", file_list)
  N.files = length(all_files)
  centers = "NA"
  centers = centers[-1]
  sample.count = "NA"
  sample.count = centers[-1]
  patient.count = "NA"
  patient.count = centers[-1]
  
  for (j in 1:N.files){
    file_j = read.csv(all_files[j], sep = "\t", stringsAsFactors = FALSE)
    centers = c(centers,paste(unique(unlist(strsplit(as.character(file_j$Center),";"))),collapse = ";"))
    sample.count = c(sample.count, length(unique(file_j$Tumor_Sample_Barcode)))
    patient.count = c(patient.count, length(unique(substring(file_j$Tumor_Sample_Barcode,1,12))))
    rm(file_j)
   }
  Data = c(length(file_list), 
           length(grep("curated", file_list)), 
           paste0(length(grep("pairs.aggregated.capture.tcga.uuid.automated", file_list)),
                 "(",grep("pairs.aggregated.capture.tcga.uuid.automated", file_list),")"), 
           paste(centers, collapse = " / "),
           paste(sample.count, collapse = " / "),
           paste(patient.count, collapse = " / "))
  SomaticMutationFiles[SomaticMutationFiles$Cancertype == Cancer, c(2:7)] = Data
}

dir.create("./1_Log_Files/2.4_MAF.selection", showWarnings = FALSE)
write.csv(SomaticMutationFiles, file = paste0("./1_Log_Files/2.4_MAF.selection/2.4.MAF.file.Analysis.", gsub(":",".",gsub(" ","_",date())), ".csv"))
