#################################################################
###
### Create RNASEQ based Deconvolution
### using Bindea's signatures and ssGSEA
###
### Input files:
### ./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData
### "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"
### Output files:
### 
###
#################################################################

#### Xcell is not included in the log file !!!!

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"                                                                # Setwd to location were output files have to be saved.

source(paste0(code_path,"R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.bioconductor.packages = c("GSVA","heatmap3", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Pancancer_matrix"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
pw_selection_version = "3.2"
Log_file = paste0("./1_Log_Files/3.8_Deconvolution_Bindea/3.8_Deconvolution_Bindea_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
#ColsideLabels = c("HML ICR clusters", "Bindea clusters")
#Legend = c("ICR Low","ICR Med","ICR High", "Bindea Low", "Bindea High")
#Legend_colors = c("blue","green","red", "pink", "purple")

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0(code_path, "Datalists/Selected.pathways.",pw_selection_version,".Rdata"))
load(paste0(code_path, "Datalists/marker_list_bindea3.RData"))

# Create folders and log file
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folders to save Rdata.files
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.8_Deconvolution_Bindea"), showWarnings = FALSE)
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
    paste0("pathway selection version used =", pw_selection_version),
    "",
    "Scripts output :",
    "",
    "Calculating Deconvolution scores",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

start.time.process.all = Sys.time()
msg = paste0("Calculating deconvolution scores and generating heatmaps", "\n")
cat(msg)

i=1
for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Calculating Deconvolution scores ",Cancer,"."))
  
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
  
  Expression.data = t(log(RNAseq.matrix +1, 2))
  available_genes = rownames(Expression.data)
  unavailable_genes_RNAseq = unlist(Bindea_ORIG)[-which(unlist(Bindea_ORIG) %in% rownames(Expression.data))]
  
  cat(paste0("Bindea ssGSEA ", Cancer, ". Total number of Bindea genes is ", length(unlist(Bindea_ORIG)), ".",
             " Of which ", length(unlist(Bindea_ORIG)[unlist(Bindea_ORIG) %in% available_genes]), 
             " genes are available in expression data."), file = Log_file, append = TRUE, sep = "\n")
  
  ## Bindea ssGSEA
  Bindea.enrichment.score = gsva(Expression.data,Bindea_ORIG,method="ssgsea")
  Bindea.enrichment.z.score = Bindea.enrichment.score 
  for(j in 1: nrow(Bindea.enrichment.z.score))  {
    Bindea.enrichment.z.score[j,] = (Bindea.enrichment.score[j,]-mean(Bindea.enrichment.score[j,]))/sd(Bindea.enrichment.score[j,]) # z-score the enrichment matrix
  }
  
  ## Patrick bindea_patrick ssGSEA
  bindea_patrick.enrichment.score = gsva(Expression.data,Bindea_PATRICK,method="ssgsea")
  bindea_patrick.enrichment.z.score = bindea_patrick.enrichment.score 
  for(j in 1: nrow(bindea_patrick.enrichment.z.score))  {
    bindea_patrick.enrichment.z.score[j,] = (bindea_patrick.enrichment.score[j,]-mean(bindea_patrick.enrichment.score[j,]))/sd(bindea_patrick.enrichment.score[j,]) # z-score the enrichment matrix
  }
  
  ## Save Scores
  dir.create(paste0("./4_Analysis/",download.method, "/", Cancer),showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment"))
  
  save(Bindea.enrichment.score, Bindea.enrichment.z.score, bindea_patrick.enrichment.score, bindea_patrick.enrichment.z.score, file = paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
                                                     "_Bindea.Rdata"))
}

