#################################################################
###
### This script re-classifies patients to High-Medium-Low ICR
### clusters as specified in 
### "./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/Cluster_assignment_analysis.csv"
### Input- & Outputfiles:
### "./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
### Cancer, "_ICR_cluster_assignment_k2-6.Rdata"
### (Rdata files get updated using this script)
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("plyr")
required.bioconductor.packages = c()
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/3.3_HMLclustering/3.3_HMLclustering_",                                                 # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data and R scripts
cluster_assignment_analysis = read.csv("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/Cluster_assignment_analysis.csv",
                                       stringsAsFactors = FALSE)
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
                                                                                                                        # in the Manual of Assembler v2.0.3 and was saved as csv file.

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create folders and log file
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.3_HMLclustering"), showWarnings = FALSE)
cat("This is a log file for creating HML clusters",                                                                     # Set-up logfile
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

cluster_assignment_analysis$mean.ICR.scores.chosen.k = NA
cluster_assignment_analysis$mean.scaled.ICR.scores.chosen.k = NA
cluster_assignment_analysis$patient.numbers.chosen.k = NA
cluster_assignment_analysis$patient.percentages.chosen.k = NA

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
    cluster_assignment_analysis[cluster_assignment_analysis$Cancertype == Cancer, -2] = NA
    next}
  
  cat(paste0("Assigning HML clusters to ", Cancer), file = Log_file, sep = "\n", append = TRUE)
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  if(cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer] == "3"){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$ICR_cluster_k3)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR1"] = "ICR Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Medium"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR High"
  }
  if(cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer] == "4"){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$ICR_cluster_k4)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR1"] = "ICR Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR4"] = "ICR High"
    if(grepl("ICR1 & ICR2", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer])){
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Low"
    }else{
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Medium"
    }
    if(grepl("ICR3 & ICR4", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer])){
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR High"
    }else{
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR Medium"
    }
  }
  if(cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer] == "5"){
    table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$ICR_cluster_k5)
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR1"] = "ICR Low"
    table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR5"] = "ICR High"
    if(grepl("ICR1 & ICR2", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer])){
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Low"
    }else{
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Medium"
    }
    if(grepl("ICR4 & ICR5", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer])){
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR4"] = "ICR High"
    }else{
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR4"] = "ICR Medium"
    }
    if(grepl("ICR1 & ICR2", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer]) & 
       grepl("ICR2 & ICR3", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer])){
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR Low"
      }
    if(grepl("ICR3 & ICR4", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer]) & 
       grepl("ICR4 & ICR5", cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer])){
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR High"
      }else{
      table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR Medium"
      }
  }
  cluster = paste0("ICR_cluster_k", cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer])
  ICR_scores = round(aggregate(ICRscore~get(cluster),data = table_cluster_assignment, FUN=mean)$ICRscore, 2)

  cluster_assignment_analysis$mean.ICR.scores.chosen.k[cluster_assignment_analysis$Cancertype == Cancer] = gsub(","," / ",toString(ICR_scores))
  scaled_ICR_scores = round(aggregate(Scaled_ICRscore~get(cluster),data = table_cluster_assignment, FUN=mean)$Scaled_ICRscore, 2)
  cluster_assignment_analysis$mean.scaled.ICR.scores.chosen.k[cluster_assignment_analysis$Cancertype == Cancer] = gsub(","," / ",toString(scaled_ICR_scores))
  
  patient_table = count(table_cluster_assignment, cluster)
  patient_table$pct = round(patient_table$freq/nrow(table_cluster_assignment)*100, 2)
  
  cluster_assignment_analysis$patient.numbers.chosen.k[cluster_assignment_analysis$Cancertype == Cancer] = gsub(",", " / ", toString(patient_table$freq))
  cluster_assignment_analysis$patient.percentages.chosen.k[cluster_assignment_analysis$Cancertype == Cancer] = gsub(",", " / ", toString(patient_table$pct))
  
  save(table_cluster_assignment,optimal.calinsky, file = Cluster_file)
}

dir.create(paste0("./4_Analysis/"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/TCGA_Assembler/"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering"), showWarnings = FALSE)
cluster.analysis.file = paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/Cluster_assignment_analysis_", 
                               gsub(":",".",gsub(" ","_",date())), ".csv")
write.csv(cluster_assignment_analysis ,file = cluster.analysis.file, row.names = FALSE)
