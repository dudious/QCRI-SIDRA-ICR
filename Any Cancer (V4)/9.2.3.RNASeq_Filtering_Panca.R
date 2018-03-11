#################################################################
###
### This script extracts tissue specific samples from 
### normalized RNASeq Data from ANY-CANCER from the TCGA database.
### Filtering to obtain a single sample for each patient
### is an option.
### 
### It will process the data into an Rdata file. 
### Data is saved in:
### "./3_Dataprocessing/",download.method, "/", Cancer, "/RNASeqData/"
###
#################################################################

##Download the RNASeq data using TCGA Assembler

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                  # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                        # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_B.R"))

required.packages = c("base64enc", "HGNChelper","RCurl","httr","stringr","digest","bitops",
                      "rjson")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                 # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                     # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"
assay.platform = "gene_RNAseq"                                                                                            
Log_file = paste0("./1_Log_Files/", download.method, "/9.2.3_RNASeq_Filtering_Panca/RNASeq_Filtering_Log_File_",                                      # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
TCGASampleTypeFile = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/TCGASampleType.txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/9.2.3_RNASeq_Filtering_Panca/"), showWarnings = FALSE)
cat("This is a log file for filtering of normalized RNASeq data",                                                       # Set-up logfile
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),                                                          
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    paste0("assay.platform = ", assay.platform),
    "",
    "Scripts output :",
    "",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

start.time.all <- Sys.time()

load(paste0("./3_DataProcessing/TCGA_Assembler/Pancancer/RNASeqData/Pancancer_gene_RNAseq_normalized.Rdata"), verbose = TRUE)

cat("Number of genes in RNASeq.NORM.quantiles is ", nrow(RNASeq.NORM.Panca.quantiles), ".\n",
    "Number of samples in RNASeq.NORM.quantiles is", ncol(RNASeq.NORM.Panca.quantiles), ".\n")

tissue.type = "TP"

filtered.norm.RNAseqData.Panca = ExtractTissueSpecificSamples(inputData = RNASeq.NORM.Panca.quantiles,
                                                              tissueType = tissue.type,
                                                              singleSampleFlag = TRUE,
                                                              sampleTypeFile = TCGASampleTypeFile)

number.deleted.samples = ncol(RNASeq.NORM.Panca.quantiles) - ncol(filtered.norm.RNAseqData.Panca)
cat("Number of deleted samples are ", number.deleted.samples, ".", "\n")

load("./3_DataProcessing/TCGA_Assembler/SKCM/RNASeqData/SKCM_gene_RNAseq_normalized_TPandTM_filtered.Rdata")
SKCM_patients = colnames(filtered.norm.RNAseqData)
TP_patients_PANCA = colnames(filtered.norm.RNAseqData.Panca)
TM_SKCM_patients = SKCM_patients[-which(SKCM_patients %in% TP_patients_PANCA)]

cat("Number of primary samples in SKCM is ", length(SKCM_patients[which(SKCM_patients %in% TP_patients_PANCA)]),
    "\nNumber of metastasis samples in SKCM is ", length(SKCM_patients[-which(SKCM_patients %in% TP_patients_PANCA)]))

### Add the 366 SKCM TM samples to the matrix

filtered.norm.RNAseqData.Panca = cbind(filtered.norm.RNAseqData.Panca, RNASeq.NORM.Panca.quantiles[, which(colnames(RNASeq.NORM.Panca.quantiles) %in%
                                                                                                             TM_SKCM_patients)])

## Check if there are any duplicated patients
duplicated.patients = substring(colnames(filtered.norm.RNAseqData.Panca)[duplicated(substring(colnames(filtered.norm.RNAseqData.Panca),1,12))], 1, 12)
  

save(filtered.norm.RNAseqData.Panca,geneInfo, file = 
"./3_DataProcessing/TCGA_Assembler/Pancancer/RNASeqData/Pancancer_gene_RNAseq_normalized_TissueType_Filtered.Rdata")
