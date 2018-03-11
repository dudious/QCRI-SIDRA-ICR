####################################################################
###
### This Script creates a complete PanCancer expression matrix
### from the expression matrixes that have been split from an
### original PanCancer expression matrix.
### 
### Input data:
### "~/Dropbox (TBI-Lab)/BREAST_QATAR/2 DATA/TCGA RNAseq/RNASeq_PANCANCER/CLEAN/TCGA.", Cancer , ".RNASeq.PANCANCER.Split.RData" files
### Output data are saved as Rdata file:
### "~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/3_DataProcessing/TCGA_Assembler/Pancancer/RNASeqData/Pancancer_gene_RNAseq_normalized_CLEAN.Rdata"
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())
"~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/3_DataProcessing/TCGA_Assembler/Pancancer/RNASeqData/Pancancer_gene_RNAseq_normalized_CLEAN.Rdata"

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

#required.packages = c("")
#ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./3_DataProcessing/",download.method, "/Pancancer"), showWarnings = FALSE)
dir.create(paste0("./3_DataProcessing/",download.method, "/Pancancer/RNASeqData"), showWarnings = FALSE)

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

## Combine all split files to a single expression matrix

i= 3
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("~/Dropbox (TBI-Lab)/BREAST_QATAR/2 DATA/TCGA RNAseq/RNASeq_PANCANCER/", Cancer,
                         ".RNASeq.TCGA.PANCANCER.SPLIT.DATA.RData"))) {next}
  load(paste0("~/Dropbox (TBI-Lab)/BREAST_QATAR/2 DATA/TCGA RNAseq/RNASeq_PANCANCER/", Cancer,
              ".RNASeq.TCGA.PANCANCER.SPLIT.DATA.RData"))
}

