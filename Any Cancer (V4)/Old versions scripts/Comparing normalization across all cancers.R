#################################################################
###
### This script makes boxplots of ICR scores in all cancers.
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("ggplot2")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Pancancer_matrix"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"

#ICR_classification_k = "HML_classification"
#basis_ordering = "Mean ICR"
#subset = "all"                                                                                                          # Subset can be "ICR High", "ICR Low", "All"
# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create folders and log file
load("./3_DataProcessing/Pancancer_matrix/ACC/RNASeqData/ACC_gene_RNAseq_normalized_TP_filtered.Rdata")
test = data.frame(matrix(nrow = nrow(filtered.norm.RNAseqData), ncol = length(CancerTYPES)), row.names = rownames(filtered.norm.RNAseqData))
colnames(test) = CancerTYPES

i=2
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {next}
  if(Cancer == "SKCM"){
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  #if(download.method == "Pancancer"){
    #load("./4_Analysis/Pancancer_matrix/Pancancer_matrix_QC/Genes_with_NA.Rdata")
    #filtered.norm.RNAseqData = filtered.norm.RNAseqData[-which(rownames(filtered.norm.RNAseqData) %in%
                                                                  #all_genes_that_have_NA_for_at_least_1_patient),]
  #}
  test[,i] = filtered.norm.RNAseqData[,1]
}

test$LAML = NULL

t_test = t(test)
if(download.method == "Pancancer_matrix"){
  t_test_log2 = t_test
}else{
  t_test_log2 = log(t_test +1, 2)
}
final = t(t_test_log2)

final = as.data.frame(final)

dev.new()
boxplot1 = ggplot(stack(final), aes(x = factor(ind, levels = names(final)), y = values)) + 
  geom_boxplot() + labs(x = "Cancers", y = "expression values") + theme(axis.title = element_text(size = 30),
                                                                        axis.text = element_text(size = 15))

print(boxplot1)
