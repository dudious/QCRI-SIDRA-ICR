#################################################################
###
### This script creates .......
###
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
Log_file = paste0("./1_Log_Files/3.2_HMLclusters/3.2_HMLclusters_",                                                     # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data and R scripts
load (paste0(code_path, "Datalists/ICR_genes.RData"))
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
msg = paste0("Create heatmaps", "\n")
cat(msg)

TCGA.cancersets$HML_k     = 0
TCGA.cancersets$optimal_k = 0

i=5
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
  next}
  #cat (paste0 ("Creating heatmap ",Cancer,"."))
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  j=3
  for(j in 3:5){
    group = paste0("ICR_cluster_k",j)
    HML_table = count(table_cluster_assignment, group)
    HML_table$pct = HML_table$freq/nrow(table_cluster_assignment)*100
    HML_table$mean_ICRs = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean)[,2]
    k=j
    if(HML_table$pct[1] > 20 & HML_table$pct[nrow(HML_table)] > 20){break}
    k=NA
  }
  TCGA.cancersets$HML_k[i] = k
  TCGA.cancersets$optimal_k[i] = optimal.calinsky
}
CancerTYPES_filtered = TCGA.cancersets$cancerType[is.na(TCGA.cancersets$HML_k)]
i=1
for (i in 1:length(CancerTYPES_filtered)) {
  Cancer = CancerTYPES_filtered[i]
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  j=4
  for(j in 4:5){
    group = paste0("ICR_cluster_k",j)
    HML_table = count(table_cluster_assignment, group)
    HML_table$ICR_cluster_k4 = as.character(HML_table$ICR_cluster_k4)
    HML_table$pct = HML_table$freq/nrow(table_cluster_assignment)*100
    HML_table$mean_ICRs = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean)[,2]
    HML_table[1,1] = "ICR_low"
    # ICR1 = low
    HML_table[nrow(HML_table),1] = "ICR_high"
    # ICRmax = high
    if(HML_table$pct[1] < 20){
      HML_table$pct[1]=HML_table$pct[1]+HML_table$pct[2]
      HML_table <- HML_table[-2,]
      # table cluter assignmnet ICR2 = low
    }
    if(HML_table$pct[nrow(HML_table)] < 20){
      HML_table$pct[nrow(HML_table)]=HML_table$pct[nrow(HML_table)]+HML_table$pct[nrow(HML_table)-1]
      HML_table <- HML_table[-c(nrow(HML_table)-1),]
      # table cluter assignmnet ICR max - 1) = high
    }
    k=j
    if(HML_table$pct[1] >= 20 & HML_table$pct[nrow(HML_table)] >= 20){break}
    k=NA
  }
TCGA.cancersets$HML_k[TCGA.cancersets$cancerType==Cancer] = k
TCGA.cancersets$optimal_k[TCGA.cancersets$cancerType==Cancer] = optimal.calinsky
}
CancerTYPES_filtered = TCGA.cancersets$cancerType[is.na(TCGA.cancersets$HML_k)]
