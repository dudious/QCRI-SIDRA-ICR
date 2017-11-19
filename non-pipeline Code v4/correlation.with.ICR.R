#################################################################
###
### Create RNASEQ based Deconvolution
### and correlation with ICR
###
### Input files:
### ./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData
### "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"
### Output files:
### 
###
#################################################################

# Setup environment
rm(list=ls())

#setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    
#code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 

library(GSVA)
library(heatmap3)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = "LAML"                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
PL.genes = c("PLCB1","PLCB2","PLCB3","PLCB4",
             "PLCD1","PLCD3","PLCD4",
             "PLCE1","PLCG1","PLCG2",
             "PLCH1","PLCH2","PLCL1","PLCL2",
             "PLCZ1","PLD1","PLD2")

if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)
byCancerResults = cbind(data.frame(cancer = CancerTYPES),matrix(nrow = length(CancerTYPES),ncol = length(PL.genes)+1))
rownames(byCancerResults) = byCancerResults$cancer
byCancerResults$cancer = NULL
colnames(byCancerResults) = c(PL.genes,"enrichment")

for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "SKCM"){
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    load (paste0("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/",
               Cancer,"_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  load (paste0("./4_Analysis/",download.method,"/",Cancer,"/Clustering/",
               Cancer,".TCGA_Assembler.EDASeq.ICR.reps5000/",
               Cancer,"_ICR_cluster_assignment_k2-6.Rdata"))
  #select data
  PL.matrix = as.data.frame(t(filtered.norm.RNAseqData[PL.genes,rownames(table_cluster_assignment)]))
  enrichment = as.numeric(gsva(filtered.norm.RNAseqData, list(PL.genes), method="ssgsea"))
  PL.matrix$enrichment = enrichment
  #corelations
  corelations = 0
  for (j in 1:ncol(PL.matrix)){
    corelations = c(corelations,cor(table_cluster_assignment$ICRscore,PL.matrix[,j],method = "pearson"))
  }
  corelations = corelations[-1]
  byCancerResults[Cancer,]=corelations
}
byCancerResults= byCancerResults[-14,]

dev.new()
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)  
heatmap3 (t(byCancerResults),
           col = my.palette,
           #main = paste0("Pan-Cancer correlation \n between ICR and Phopholipase"),
           Rowv = NA,
           scale = "none",
           cex.main = 0.7,
           cexCol = 1.2,
           cexRow = 1.2,
           margins=c(12, 15))
