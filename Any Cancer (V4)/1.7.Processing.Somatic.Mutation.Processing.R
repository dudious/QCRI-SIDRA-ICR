#################################################################
###
### This script processes ANY-CANCER Somatic Mutation Data 
### from the TCGA database.
### MAF-files for processing are selected based on the output of script 2.4. 
### It will process the data. 
### Processed data is saved in:
###("./3_DataProcessing/",download.method,"/",Cancer,"/SomaticMutationData/")
###
#################################################################

##Download the SomaticMutationData data using TCGA Assembler

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())
#setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")
setwd("D:/Jessica/Dropbox (TBI-Lab)/TCGA Analysis pipeline/") 
required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper")


# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)


#Path.R.Tools = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/R tools/"                               # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Path.R.Tools = "D:/Jessica/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/R tools/"
Log_file = paste0("./1_Log_Files/1.7_SomaticMutation_Processing/SomaticMutation_Processing_Log_File_",                  # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
support_folder = paste0(Path.R.Tools, "TCGA-Assembler_v2.0.3/SupportingFiles")

# Load data
TCGA.cancersets = read.csv ("./TCGA.datasets.csv",stringsAsFactors = FALSE)                                             # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
source(paste0(Path.R.Tools, "TCGA-Assembler_v2.0.3/Module_A.R"))
source(paste0(Path.R.Tools, "TCGA-Assembler_v2.0.3/Module_B.R"))
source(paste0(Path.R.Tools, "ipak.function.R"))                                                                         


#Install and load required packages
ipak(required.packages)

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/1.7_SomaticMutation_Processing/"), showWarnings = FALSE)
cat("This is a log file for the Processing of Somatic Mutation data",                                                   # Set-up logfile
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


for (i in 1:N.sets) {
  start.time.loop.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  print (paste0 ("Processing ",Cancer,"."))
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/SomaticMutationData")
  file_list_all = list.files(Cancer_path, full.names = TRUE)
  pairs_aggregated_file = file_list_all[grep("pairs.aggregated.capture.tcga.uuid.automated", file_list_all)]
  
  dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer),showWarnings = FALSE)
  dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer,"/SomaticMutationData"),showWarnings = FALSE)
  
  folder = paste0("./3_DataProcessing/",download.method,"/",Cancer,"/SomaticMutationData")
  
  N.files = length(pairs_aggregated_file)
  if (N.files>0) {
    for(k in 1:N.files){
      start.time.process.cancer = Sys.time()
      file = pairs_aggregated_file[k]
      print(paste0("Processing ",file,"."))
      outputname = paste0(Cancer, "_", "Pairs_Aggregated_Capture_Processed")
      ProcessSomaticMutationData(inputFilePath = file,
                                 outputFileName = outputname,
                                 outputFileFolder = folder)
      outputfiles <- list.files(folder, full.names = TRUE)
      if(length(grep(".txt",outputfiles))){
        file.remove(outputfiles[grep(".txt",outputfiles)])
      }
      if(length(grep(".rda",outputfiles))){
        old.name <- outputfiles[grep(".rda",outputfiles)]
        new.name <- gsub(".rda", ".Rdata", old.name)
        file.rename(old.name, new.name)
      }
      end.time.process.cancer <- Sys.time ()
      time = substring(as.character(capture.output(round(end.time.process.cancer - start.time.process.cancer, 2))),20,100)
      msg = paste0("Processing time for ", file, ": ",time, ".", "\n", "Outputfile is ", outputname, ".Rdata.")
      cat(msg)
      cat(msg, file = Log_file, sep = "\n",append=TRUE)
    }
  }
}
end.time.process.all <- Sys.time ()
time <- substring(as.character(capture.output(round(end.time.process.all - start.time.process.all, 2))),20,100)
msg = paste0("\n","Processing time for all cancertypes: ",time, " min.", "\n", "---------------------------------------------------------------")
cat(msg)
cat(msg, file = Log_file,"",sep = "\n",append=TRUE)