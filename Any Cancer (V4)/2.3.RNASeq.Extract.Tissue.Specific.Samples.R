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
Log_file = paste0("./1_Log_Files/2.3_RNASeq_Filtering/RNASeq_Filtering_Log_File_",                                      # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
TCGASampleTypeFile = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/TCGASampleType.txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/2.3_RNASeq_Filtering/"), showWarnings = FALSE)
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

# Filtering
for (i in 1:N.sets) {
  start.time.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  cat(paste0(Cancer, ": "), file = Log_file, append = TRUE, sep = "\n")
  if (Cancer %in% Cancer_skip) {next}
  tissue.type = "TP"
  if (Cancer == "SKCM") {tissue.type = c("TP", "TM")}
  print(paste0 ("Filtering ",Cancer,"."))
  
  sink(Log_file, append=TRUE, split=TRUE)
  
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
    {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
                "\n",
                "-----------------------------------------------------------------------------------------------------------",
                "\n"))
    sink()
    next}
  
  load(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
              assay.platform, "_", "normalized.Rdata"), verbose = TRUE)
  
  cat("Number of genes in RNASeq.NORM.quantiles is ", nrow(RNASeq.NORM.quantiles), ".\n",
      "Number of samples in RNASeq.NORM.quantiles is", ncol(RNASeq.NORM.quantiles), ".\n")
  
  filtered.norm.RNAseqData = ExtractTissueSpecificSamples(inputData = RNASeq.NORM.quantiles,
                                                         tissueType = tissue.type,
                                                         singleSampleFlag = TRUE,
                                                         sampleTypeFile = TCGASampleTypeFile)
  
  number.deleted.samples = ncol(RNASeq.NORM.quantiles) - ncol(filtered.norm.RNAseqData)
  cat("Number of deleted samples are ", number.deleted.samples, ".", "\n")
  
  
  if(length(unique(substring(colnames(filtered.norm.RNAseqData),1,12))) == length(colnames(filtered.norm.RNAseqData))){
    cat("Number of samples is equal to total number of patients.\n")
  } else
    {cat("Number of samples is ", length(colnames(filtered.norm.RNAseqData)), " and number of patients is ",
           length(unique(substring(colnames(filtered.norm.RNAseqData),1,12))), ".\n")
      
      duplicated.patients = substring(colnames(filtered.norm.RNAseqData)[duplicated(substring(colnames(filtered.norm.RNAseqData),1,12))], 1, 12)
      filtered.norm.RNAseqData = filtered.norm.RNAseqData[,-which(substring(colnames(filtered.norm.RNAseqData),1,12) %in% duplicated.patients & substring(colnames(filtered.norm.RNAseqData),14,15) == "06")]}
  
  sink()
  
  if(Cancer== "SKCM"){Rdata.file = paste0("./3_Dataprocessing/",download.method, "/", Cancer, "/RNASeqData/",Cancer, "_", 
                                                       assay.platform, "_normalized_TPandTM_filtered.Rdata")}
  else{Rdata.file = paste0("./3_Dataprocessing/",download.method, "/", Cancer, "/RNASeqData/",Cancer, "_", 
                      assay.platform, "_normalized_", tissue.type, "_filtered.Rdata")}
  
  save(filtered.norm.RNAseqData,geneInfo, file= Rdata.file)
  
  
  end.time.cancer = Sys.time()
  time = substring(as.character(capture.output(round(end.time.cancer - start.time.cancer, 2))),20,100)
  msg = paste0("Filtering time for ", Cancer, ": ",time, ".", "\n", "Outputfile is ", Rdata.file, "\n",
               " which contains RNASeq.NORM.quantiles.filtered and geneInfo.", "\n",
               "-----------------------------------------------------------------------------------------------------------")
  cat(msg)
  cat(msg, file= Log_file,sep = "\n",append=TRUE)
}

end.time.all = Sys.time ()
time = substring(as.character(capture.output(round(end.time.all - start.time.all, 2))),20,100)
msg = paste0("Filtering time for all cancertypes: ",time, "\n", 
             "---------------------------------------------------------------")
cat(msg)   
cat(msg, file= Log_file,sep = "\n",append=TRUE)
