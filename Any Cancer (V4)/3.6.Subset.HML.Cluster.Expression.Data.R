#################################################################
###
### This script extracts expression data for samples in individual
### ICR clusters.
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c()
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/3.6_Subset.HML.Expression.Data/3.6_Subset_Expressiondata_Log_File_",                   # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")                                                            # Specify version of manual correction to perform in this script
version = "v3"

# Load data
load (paste0(code_path, "Datalists/ICR_genes.RData"))
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders and log file
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.6_Subset.HML.Expression.Data"), showWarnings = FALSE)
cat("This is a log file for subsetting expression data on HML cluster",                                       # Set-up logfile
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
    paste0("version of manual correction file for HML allocation used to generate heatmaps = ", version),
    "",
    "Scripts output :",
    "",
    "Creating heatmaps",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create heatmaps
start.time.process.all = Sys.time()
msg = paste0("Making subsets", "\n")
cat(msg)

for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Subsetting ",Cancer,"."))
  
  ## load RNASeq data
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  if(Cancer == "SKCM"){
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  samples_with_med_or_low_ICR = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR Medium" | table_cluster_assignment$HML_cluster == "ICR Low")]
  Med_Low_ICR_RNAseqData = filtered.norm.RNAseqData[, samples_with_med_or_low_ICR]
  
  samples_with_low_ICR = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR Low")]
  Low_ICR_RNAseqData = filtered.norm.RNAseqData[, samples_with_low_ICR]
  
  samples_with_med_ICR = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR Medium")]
  Med_ICR_RNAseqData = filtered.norm.RNAseqData[, samples_with_med_ICR]
  
  samples_with_high_ICR = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR High")]
  High_ICR_RNAseqData = filtered.norm.RNAseqData[, samples_with_high_ICR]
  
  save(Med_Low_ICR_RNAseqData, Low_ICR_RNAseqData, Med_ICR_RNAseqData, High_ICR_RNAseqData, 
       file = paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_TP_filtered_ICR_subsetted.Rdata"))
}
  
  
