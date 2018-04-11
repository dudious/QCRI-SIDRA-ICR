####################################################################
###
### This Script copy pastes the Rdata files in the TCGA Assembler 
### folder with table cluster assignment and ICR score to the 
### Pancancer_matrix folder.
### Unfinished script!
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/stefanofunctions.R"))                                                                 # Used for calinsky function and plot

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper")
required.bioconductor.packages = "ConsensusClusterPlus"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                    # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                    # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Pancancer_matrix"                                                                                   # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
Log_file = paste0("./1_Log_Files/3.1.2_Consensus_Clustering/3.1.2_Consensus_Clustering_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.1.2_Consensus_Clustering/"), showWarnings = FALSE)
cat("This is a log file for Consensus Clustering of RNASeq data",                                                       # Set-up logfile
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
    "Clustering",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

## Perform clusting
start.time.process.all = Sys.time()
msg = paste0("Clustering", "\n")
cat(msg)

i=1
for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Clustering ",Cancer,"."))
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  dir.create(paste0("./4_Analysis"), showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/", download.method), showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/", download.method, "/", Cancer), showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering"), showWarnings = FALSE)
  load(paste0("./4_Analysis/TCGA_Assembler/", Cancer, "/Clustering/", Cancer, ".TCGA_Assembler.EDASeq.ICR.reps5000/", Cancer, "_ICR_cluster_assignment_k2-6.Rdata"))
  load(paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  ICR_subset_RNAseq_log2 = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% ICR_genes, ])
  table_cluster_assignment$ICRscore_Pancancer = rowMeans(ICR_subset_RNAseq_log2)
