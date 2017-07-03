#################################################################
###
### This script downloads the Clinical Data and
### Biospecimen Data from the TCGA database.
### 
### Data is saved :
### .../2_Data/",download.method,"/",Cancer,"/BiospecimenClinicalData/
###
#################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# and indicate the location.

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                  # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                        # Set code path to the location were the R code is located.

source(paste0(code_path, "R tools/ipak.function.R")) 
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_A.R"))

required.packages <- c("xlsx","RCurl","httr", "rjson", "stringr", "HGNChelper")
ipak(required.packages)                                                                                               # Install and load required packages     


# Set Parameters
download.method = "TCGA_Assembler"                                                                                    # Specify download method (this information to be used when saving the file)
CancerTYPES     = "ALL"
Cancer_skip     = c("")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)

# Download clinical and Biospecimen information in the Biotab format
if (CancerTYPES == "ALL") { 
      CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets <- length(CancerTYPES)

for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/BiospecimenClinicalData/")
  DownloadBiospecimenClinicalData(cancerType = Cancer,
                                  saveFolderName = Cancer_path,
                                  outputFileName = "")
  print (paste0("Clinical and Biospecimen data downloaded to ",Cancer_path))
  # Extract and merge relevant data
  ## rename files
  file.list.old  <- list.files(Cancer_path,full.names = TRUE)
  file.list.new <- gsub ("nationwidechildrens.org_clinical_","",file.list.old)
  file.list.new <- gsub (paste0("_",tolower(Cancer)),"",file.list.new)
  file.rename (file.list.old,file.list.new)
  print ("Files renamed ...")
  

for (j in 1:length(file.list.new)){
  file = file.list.new[j]
  name = gsub(paste0("./2_Data/TCGA_Assembler/", Cancer, "BiospecimenClinicalData/"),"",file)
  name = paste0(gsub(".txt","",name),".table")
  table = read.csv(file,header = TRUE, sep="\t", as.is=TRUE,skip=2)
  assign(name,table)
  rm(table)
}
}

