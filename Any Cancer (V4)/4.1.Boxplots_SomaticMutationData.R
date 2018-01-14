#################################################################
###
### This script generates boxplots for ANY-CANCER Somatic Mutation Data 
### from the TCGA database.
### 
### Input data:
###("./3_DataProcessing/",download.method,"/",Cancer,"/SomaticMutationData/")
### Output dat:
###("./5_Figures/SomaticMutationBoxplots/", download.method, "/")
#################################################################

## Process SomaticMutationData data using TCGA Assembler

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                   # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                         # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))
required.packages = c("ggplot2")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
Log_file = paste0("./1_Log_Files/4.1_SomaticMutation_Boxplots/SomaticMutation_Boxplots_Log_File_",                      # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)

# Create folders
dir.create("./5_Figures",showWarnings = FALSE)                                                                          # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./5_Figures/SomaticMutationBoxplots"),showWarnings = FALSE)
dir.create(paste0("./5_Figures/SomaticMutationBoxplots/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/4.1_SomaticMutation_Boxplots/"), showWarnings = FALSE)
cat("This is a log file for the generation of boxplots with Somatic Mutation data",                                                   # Set-up logfile
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
    "",
    "Scripts output :",
    "",
    "Downloads",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

start.time.script <- Sys.time()

## Process data
start.time.process.all = Sys.time()
msg = paste0("Processing", "\n")
cat(msg, file= Log_file,sep = "\n",append=TRUE)

i=1
for (i in 1:N.sets) {
  start.time.loop.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  print (paste0 ("Making boxplots for ",Cancer,"."))
  load(paste0("./3_DataProcessing/",download.method,"/",Cancer,"/SomaticMutationData/", Cancer, "_Pairs_Aggregated_Capture_Processed_mutationLevel.Rdata"))
  
  # Generate mutation matrix
  rownames(Data) = Des$GeneSymbol
  N.rows = nrow(Data)
  N.columns = ncol(Data)
  
  j=1
  k=1
  
  for(j in 1:N.rows){
    for(k in 1:N.columns){
      if(Data[j, k] == 1){
        Data[j,k] = Des$Variant_Classification[j]
      }
      if(Data[j,k] == 0){
        Data[j,k] = ""
      }
    }
  }
  

end.time.process.all <- Sys.time ()
time <- substring(as.character(capture.output(round(end.time.process.all - start.time.process.all, 2))),20,100)
msg = paste0("\n","Processing time for all cancertypes: ",time, " min.", "\n", "---------------------------------------------------------------")
cat(msg)
cat(msg, file = Log_file,"",sep = "\n",append=TRUE)