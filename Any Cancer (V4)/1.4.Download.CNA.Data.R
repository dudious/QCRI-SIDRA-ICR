#################################################################
###
### This script downloads ANY-CANCER copy number Data 
### from the TCGA database
### It will download and process the data. 
### Downloaded raw data is saved in:
### "./2_Data/",download.method,"/",Cancer,"/CNAData/"
### Processed data is saved in:
###("./3_DataProcessing/",download.method,"/",Cancer,"/CNAData/")
###
#################################################################

##Download the CNA/CNV data using TCGA Assembler

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                     # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                           # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_A.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_B.R"))

required.packages <- c("RCurl","httr", "rjson", "stringr", "HGNChelper")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
Log_file = paste0("./1_Log_Files/1.4_CNA_Processing/CNA-Process_Log_File_",
                  gsub(":",".",gsub(" ","_",date())),".txt")
GenomeFileHg18 = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/Hg18GenePosition.txt")                 # Load RefGenomeFiles from Assembler_v2.0.3
GenomeFileHg19 = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/Hg19GenePosition.txt")                 # Load RefGenomeFiles from Assembler_v2.0.3

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
                                                                                                                         # in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/1.4_CNA_Processing/"), showWarnings = FALSE)
cat("This is a log file for Processing CNV Data",
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    file = Log_file,
    append = FALSE, sep= "\n")

N.sets = length(CancerTYPES)

#Download data
start.time <- Sys.time ()
for (i in 1:N.sets) {
  if(i %in% c(1:19)){next}
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/CNAData/")
  DownloadCNAData(cancerType = CancerTYPES[i],
                  assayPlatform = NULL,
                  tissueType = NULL,
                  saveFolderName = Cancer_path,
                  outputFileName = "",
                  inputPatientIDs = NULL)
  print (paste0("CNA Data are downloaded to ",Cancer_path))
}
end.time <- Sys.time ()
time <- end.time - start.time
print (time)   

## Process CNV data (See Manual: function does the following: 1. calculate gene-level copy number value, 2. check and correct the gene identifiers to official gene symbols, 
# 3. Draw and save a box plot of gene-level copy number data for QC, 4. Save tab-delimited .txt file, 5. Save data as R data file .rda)

start.time <- Sys.time ()
for (j in 1:N.sets) {
  start.time.loop.cancer = Sys.time()
  Cancer = CancerTYPES[j]
  if (Cancer %in% Cancer_skip) {next}
  print (paste0 ("Processing ",Cancer,"."))
  cat(paste0("Cancer : ",Cancer), file= Log_file,sep = "\n",append=TRUE)
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/CNAData/")
  file.list.all <- list.files(Cancer_path, full.names = TRUE)
  dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer),showWarnings = FALSE)
  dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer,"/CNAData/"),showWarnings = FALSE)
  folder = paste0("./3_DataProcessing/",download.method,"/",Cancer,"/CNAData/")
  
  N.files = length(file.list.all)
  if (N.files>0) {
    for(i in 1:N.files){
      start.time.loop.file = Sys.time()
      file = file.list.all[i]
      print(paste0("Processing ",file,"."))
      if(length(grep("hg19",file))){ref.file = GenomeFileHg19}
      if(length(grep("hg18",file))){ref.file = GenomeFileHg18}
      location.input = paste0("./2_Data/TCGA_Assembler/", Cancer, "/CNAData/")
      outputname = gsub(location.input, "", file)
      outputname = gsub(".txt", "_Processed_GeneLevel", outputname)
      ProcessCNAData(inputFilePath = file,
                    outputFileName = outputname,
                    outputFileFolder = folder,
                    refGenomeFile = ref.file)
      end.time <- Sys.time ()
      time <- end.time - start.time.loop.file
      msg = paste0("Processing time for ", file, ":",time,"\n", "outputfile is ", folder,outputname, ".")
      cat (msg)
      cat(msg, file = Log_file, sep = "\n",append=TRUE)
    }
  }
  end.time <- Sys.time ()
  time <- end.time - start.time.loop.cancer
  msg = paste0("Processing time for ",Cancer," :",time)
  print(msg)
  cat(msg,"", file= Log_file,sep = "\n",append=TRUE)
}
end.time <- Sys.time ()
time <- end.time - start.time
msg = paste0("Processing time for script :",time)
print(msg)
cat(msg, file= Log_file,sep = "\n",append=TRUE)


