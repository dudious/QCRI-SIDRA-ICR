####################################################################
###
### This Script clusters the RNASeq date from TCGA that have been
### normalized using EDASeq and filtered for specific tissue type
### and 1 sample per patient. This script includes a 
### log transformation of the data. Additionally, optimal kalinsky
### is calculated.
### 
### Input data:
### ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/")
### Output data are saved as Rdata file:
### ./, Cancer, ".", download.method, ".EDASeq.ICR.reps5000/" , Cancer, "_ICR_cluster_assignment_k2-6.Rdata"
### which includes: (1) table_cluster_assignment and (2) optimal.calinsky.
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
Cancer_skip = c("")                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
Log_file = paste0("./1_Log_Files/3.1_Consensus_Clustering/3.1_Consensus_Clustering_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/3.1_Consensus_Clustering"), showWarnings = FALSE)
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
  if(Cancer == "SKCM"){
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  
  ICR_subset_RNAseq = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% ICR_genes, ])
  ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)
  
  dir.create(paste0("./4_Analysis/",download.method,"/",Cancer),showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/",download.method,"/",Cancer,"/Clustering"),showWarnings = FALSE)
  
  setwd(paste0(getwd(),"/4_Analysis/",download.method,"/",Cancer,"/Clustering"))
  
  
  ddist = dist(ICR_subset_RNAseq_log2)
  
  ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                                 maxK = 6,                                                                              # set K
                                                 pItem = 0.8,
                                                 reps=5000,                                                                             # set repeats
                                                 title=paste0(Cancer,".", download.method,".EDASeq.ICR.reps5000"),              # Output filename (no path)
                                                 clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                                 innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                                 finalLinkage = "complete",
                                                 plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                                 writeTable = TRUE,
                                                 verbose = TRUE)
  outputfiles = list.files(paste0(Cancer,".", download.method,".EDASeq.ICR.reps5000"), full.names = TRUE)
  class_files = outputfiles[grep("consensusClass", outputfiles)]
  
  N.files = length(class_files)
  table_cluster_assignment = data.frame(ICRscore = rowMeans(ICR_subset_RNAseq_log2))
  
  for (j in 1:N.files){
    file = paste0("./", class_files[j])
    consensus_class = read.csv(file = file,header=FALSE)
    group = paste0("Group_k",j+1)
    colnames(consensus_class) = c("PatientID", group)
    rownames(consensus_class) = consensus_class$PatientID
    consensus_class$PatientID = NULL
    table_cluster_assignment[,group] = consensus_class[,group][match(rownames(table_cluster_assignment), rownames(consensus_class))]
    
    transl_table_ICR_cluster = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean)
    colnames(transl_table_ICR_cluster) = c(group,"mean_ICRscore")
    transl_table_ICR_cluster = cbind(transl_table_ICR_cluster[order(transl_table_ICR_cluster$mean_ICRscore),],ICR_name=paste0("ICR",c(1:(j+1))))
    
    ICR_cluster = paste0("ICR_cluster_k",j+1)
    table_cluster_assignment[,ICR_cluster] = transl_table_ICR_cluster$ICR_name[match(table_cluster_assignment[,group],
                                                                                          transl_table_ICR_cluster[,group])]
  }
  
  #calinsky
  sHc <- hclust(ddist, method = "ward.D2")
  aCalinsky <- calinsky(sHc,gMax=10)
  pdf(file = paste0("./", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/" , Cancer, "_ICR_cluster_assignment_k2-6.Calinsky.pdf"), width = 16, height = 6)
  plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
  text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
  dev.off()
  optimal.calinsky = which(aCalinsky == max(aCalinsky[3:5]))
  
  #save data
  outputname = paste0("./", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/" , Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  save(table_cluster_assignment,optimal.calinsky, file = outputname)
  
  #log calinsky
  setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")
  cat(paste0("Optimal calinsky for ", Cancer, " is: ", optimal.calinsky, "."), file = Log_file, append = TRUE, sep = "\n")
  
  end.time.process.cancer <- Sys.time ()
  time = substring(as.character(capture.output(round(end.time.process.cancer - start.time.process.cancer, 2))),20,100)
  msg = paste0("Processing time for ", file, ": ",time, ".", "\n", "Outputfile is ", outputname)
  cat(msg)
  cat(msg, file = Log_file, sep = "\n",append=TRUE)
}
end.time.process.all <- Sys.time ()
time <- substring(as.character(capture.output(round(end.time.process.all - start.time.process.all, 2))),20,100)
msg = paste0("\n","Processing time for all cancertypes: ",time, " min.", "\n", "---------------------------------------------------------------")
cat(msg)
cat(msg, file = Log_file,"",sep = "\n",append=TRUE)
