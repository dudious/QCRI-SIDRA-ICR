#################################################################
###
### This script creates an overview of the ICR cluster assignment
### to decide on how to generate a High-Medium-Low ICR classification
### for all cancer types.
###
### Outputfile:
### "./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/Cluster_assignment_analysis
### gsub(":",".",gsub(" ","_",date()))," .csv"
### This file can be manually corrected before HML-classification in
### script 3.3.
#################################################################

# Setup environment
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
Log_file = paste0("./1_Log_Files/3.2_HMLclusters/3.2_HMLclusters_",                                                     # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
                                                                                                                        # in the Manual of Assembler v2.0.3 and was saved as csv file.

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create folders and log file
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.2_HMLclusters"), showWarnings = FALSE)
cat("This is a log file to decide how to allocate samples to HML clusters",                                             # Set-up logfile
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
    "Creating HML clusters",
    file = Log_file,
    append = FALSE, sep= "\n")

start.time.process.all = Sys.time()
msg = paste0("Processing", "\n")
cat(msg)

cluster_assignment_analysis = data.frame(Cancertype = CancerTYPES, Optimal.K = NA, mean.ICR.scores = NA, 
                                         mean.scaled.ICR.scores = NA, patient.numbers = NA, patient.percentages = NA, auto.action = NA, chosen.K = NA, 
                                         manual.correction = "no", reason.manual.correction = NA, manual.action = NA)

for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    cluster_assignment_analysis[cluster_assignment_analysis$Cancertype == Cancer, -1] = NA
    next}
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  table_cluster_assignment$Scaled_ICRscore = round((100-1)/(max(table_cluster_assignment$ICRscore)-min(table_cluster_assignment$ICRscore))*(table_cluster_assignment$ICRscore-max(table_cluster_assignment$ICRscore))+100,1)
  cluster_assignment_analysis$Optimal.K[cluster_assignment_analysis$Cancertype == Cancer] = optimal.calinsky
  cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer] = optimal.calinsky
  
  cluster = paste0("ICR_cluster_k", cluster_assignment_analysis$Optimal.K[cluster_assignment_analysis$Cancertype == Cancer])
  ICR_scores = round(aggregate(ICRscore~get(cluster),data = table_cluster_assignment, FUN=mean)$ICRscore, 2)
  cluster_assignment_analysis$mean.ICR.scores[cluster_assignment_analysis$Cancertype == Cancer] = gsub(","," / ",toString(ICR_scores))
  scaled_ICR_scores = round(aggregate(Scaled_ICRscore~get(cluster),data = table_cluster_assignment, FUN=mean)$Scaled_ICRscore, 2)
  cluster_assignment_analysis$mean.scaled.ICR.scores[cluster_assignment_analysis$Cancertype == Cancer] = gsub(","," / ",toString(scaled_ICR_scores))
  
  patient_table = count(table_cluster_assignment, cluster)
  patient_table$pct = round(patient_table$freq/nrow(table_cluster_assignment)*100, 2)
  
  cluster_assignment_analysis$patient.numbers[cluster_assignment_analysis$Cancertype == Cancer] = gsub(",", " / ", toString(patient_table$freq))
  cluster_assignment_analysis$patient.percentages[cluster_assignment_analysis$Cancertype == Cancer] = gsub(",", " / ", toString(patient_table$pct))
  
  if(cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer] == "3"){
    cluster_assignment_analysis$auto.action[cluster_assignment_analysis$Cancertype == Cancer] = "none"
  } else{
    combine_location = which.min(diff(ICR_scores))
    all_ICR_clusters = as.character(patient_table[,cluster])
    ICR_clusters_to_combine = c(all_ICR_clusters[combine_location], all_ICR_clusters[combine_location+1])
    cluster_assignment_analysis$auto.action[cluster_assignment_analysis$Cancertype == Cancer] = paste0("Combine clusters ", paste(ICR_clusters_to_combine,collapse = " & "), ".")
  }
}

dir.create(paste0("./4_Analysis/"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/TCGA_Assembler/"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering"), showWarnings = FALSE)
write.csv(cluster_assignment_analysis ,file = paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/Cluster_assignment_analysis",
                                                     gsub(":",".",gsub(" ","_",date())), ".csv"))
