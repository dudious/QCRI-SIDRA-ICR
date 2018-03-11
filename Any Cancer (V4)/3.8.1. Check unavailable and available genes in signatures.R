#################################################################
###
### Create RNASEQ based Deconvolution
### using Bindea's signatures and ssGSEA. Check gene availability.
###
### Output files:
### "./1_Log_Files/3.8_Deconvolution_Bindea/3.8_Available_Unavailable_genes_",                          # Specify complete name of the logfile that will be saved during this script
### gsub(":",".",gsub(" ","_",date())),".txt"
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path,"R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

if(!("xCell" %in% installed.packages()[,"Package"])){
  required.packages = c("devtools")
  ipak(required.packages)
}
required.bioconductor.packages = c("GSVA","heatmap3")                                                                   
ibiopak(required.bioconductor.packages)
devtools::install_github('dviraran/xCell')
library ("xCell")

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
pw_selection_version = "3.2"
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("HML ICR clusters", "Bindea clusters")
Legend = c("ICR Low","ICR Med","ICR High", "Bindea Low", "Bindea High")
Legend_colors = c("blue","green","red", "pink", "purple")
Log_file = paste0("./1_Log_Files/", download.method, "/3.8_Deconvolution_Bindea/3.8_Available_Unavailable_genes_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

cat("This is a log file for Deconvolution using Bindeas gene signatures, xCell and Hallmark pathways on RNASeq data",   # Set-up logfile
    "_________________________________________________________________________________",
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
    "Calculating Deconvolution scores",
    file = Log_file,
    append = FALSE, sep= "\n")

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0(code_path, "Datalists/Selected.pathways.",pw_selection_version,".Rdata"))
load(paste0(code_path, "Datalists/marker_list_bindea2.RData"))


# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

Cancer = "ACC"
if (Cancer %in% Cancer_skip) {next}
  
  ## load RNASeq data
  
Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  
  ## load cluster data
Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
load(Cluster_file)
  
Expression.data = log(filtered.norm.RNAseqData +1, 2)
available_genes = rownames(Expression.data)
available_bindea_genes = unlist(marker_list)[which(unlist(marker_list) %in% rownames(Expression.data))]
unavailable_bindea_genes_RNAseq = unlist(marker_list)[-which(unlist(marker_list) %in% rownames(Expression.data))]
  
cat(paste0("Bindea ssGSEA ", Cancer, ". Total number of Bindea genes is ", length(unlist(marker_list)), ".",
             " Of which ", length(unlist(marker_list)[unlist(marker_list) %in% available_genes]), 
             " genes are available in expression data."), file = Log_file, append = TRUE, sep = "\n")
cat(paste0("\nAvailable bindea genes = "), file = Log_file, append = TRUE, sep = "\n")
cat(paste0(toString(available_bindea_genes)), file = Log_file, append = TRUE, sep = "\n")
cat(paste0("\nUnavailable bindea genes = \n", toString(unavailable_bindea_genes_RNAseq)), file = Log_file, append = TRUE, sep = "\n")
  
available_Hallmark_genes_RNAseq = unlist(Selected.pathways)[which(unlist(Selected.pathways) %in% rownames(Expression.data))]
unavailable_Hallmark_genes_RNAseq = unlist(Selected.pathways)[-which(unlist(Selected.pathways) %in% rownames(Expression.data))]
  
cat(paste0("Total number of hallmark genes is ", length(unlist(Selected.pathways)), ".",
             " Of which ", length(unlist(Selected.pathways)[unlist(Selected.pathways) %in% available_genes]), 
             " genes are available in expression data.\n"), file = Log_file, append = TRUE, sep = "\n")
cat(paste0("\nAvailable Hallmark genes = \n", toString(available_Hallmark_genes_RNAseq)), file = Log_file, append = TRUE, sep = "\n")
cat(paste0("\nUnavailable Hallmark genes = \n", toString(unavailable_Hallmark_genes_RNAseq)), file = Log_file, append = TRUE, sep = "\n")

  