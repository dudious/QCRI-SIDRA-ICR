#################################################################
###
### This script adds a scaled ICR gene expression variable to the 
### table cluster assignment table in 
### "./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", 
### download.method, ".EDASeq.ICR.reps5000/", Cancer, "_ICR_cluster_assignment_k2-6.Rdata".
### The initial running of this script includes the re-ordering of
### columns (table_cluster_assignment = table_cluster_assignment[,c(1, 12, 2:11)]).
### 
#################################################################

rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("plyr")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/3.1_Consensus_Clustering/3.1.1_Adding_Scaled_Gene_Expression_Log_File_",               # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Log file
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.1_Consensus_Clustering/"), showWarnings = FALSE)
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

for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}

cat(paste0("Adding scaled gene expression to ", Cancer, ".", file = Log_file, append = TRUE, sep = "\n"))
Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                      Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
load(Cluster_file)

table_cluster_assignment$Scaled_ICRscore = round((100-1)/(max(table_cluster_assignment$ICRscore)-min(table_cluster_assignment$ICRscore))*(table_cluster_assignment$ICRscore-max(table_cluster_assignment$ICRscore))+100,1)
#table_cluster_assignment = table_cluster_assignment[,c(1, 12, 2:11)] ## Only re-order after initial run of the 3.1.1 script directly after clustering
                                                                      ## Formula for scaling can be adapted to revise scaling, but no extra columns should be added and columns should maintain their order to keep full pipeline functional.
save(table_cluster_assignment,optimal.calinsky, file = Cluster_file)
}
