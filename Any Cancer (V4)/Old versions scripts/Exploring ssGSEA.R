rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path,"R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.bioconductor.packages = c("GSVA","heatmap3")
ibiopak(required.bioconductor.packages)

load(paste0(code_path, "Datalists/Gene_Lists.Rdata"))
load(paste0(code_path, "Datalists/marker_list_bindea2.RData"))
download.method = "TCGA_Assembler"

Cancer = "ACC"
Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                      Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
load(Cluster_file)

Expression.data = log(filtered.norm.RNAseqData +1, 2)

marker_list$ICR_genes = Gene_Lists$ICR_genes

enrichment.score.ICR = gsva(Expression.data, marker_list[25], method= "ssgsea")
enrichment.score.ICR.B = gsva(Expression.data, marker_list[c(1,25)],method="ssgsea")
enrichment.score.ICR.Bindea = gsva(Expression.data, marker_list, method = "ssgsea")

marker_list$ICR_genes2 = Gene_Lists$ICR_genes

enrichment.score.ICR.ICR.Bindea = gsva(Expression.data, marker_list, method = "ssgsea")

enrichment.z.score.ICR = enrichment.score.ICR 
for(j in 1: nrow(enrichment.z.score.ICR))  {
  enrichment.z.score.ICR[j,] = (enrichment.score.ICR[j,]-mean(enrichment.score.ICR[j,]))/sd(enrichment.score.ICR[j,]) # z-score the enrichment matrix
}

enrichment.z.score.ICR.B = enrichment.score.ICR.B 
for(j in 1: nrow(enrichment.z.score.ICR.B))  {
  enrichment.z.score.ICR.B[j,] = (enrichment.score.ICR.B[j,]-mean(enrichment.score.ICR.B[j,]))/sd(enrichment.score.ICR.B[j,]) # z-score the enrichment matrix
}

enrichment.z.score.ICR.Bindea = enrichment.score.ICR.Bindea 
for(j in 1: nrow(enrichment.z.score.ICR.Bindea))  {
  enrichment.z.score.ICR.Bindea[j,] = (enrichment.score.ICR.Bindea[j,]-mean(enrichment.score.ICR.Bindea[j,]))/sd(enrichment.score.ICR.Bindea[j,]) # z-score the enrichment matrix
}

enrichment.z.score.ICR.ICR.Bindea = enrichment.score.ICR.ICR.Bindea 
for(j in 1: nrow(enrichment.z.score.ICR.ICR.Bindea))  {
  enrichment.z.score.ICR.ICR.Bindea[j,] = (enrichment.score.ICR.ICR.Bindea[j,]-mean(enrichment.score.ICR.ICR.Bindea[j,]))/sd(enrichment.score.ICR.ICR.Bindea[j,]) # z-score the enrichment matrix
}

Sample_high_ICR = which.max(enrichment.score.ICR[1,])
Sample_low_ICR = which.min(enrichment.score.ICR[1,])
enrichment.score.ICR.B[1,Sample_high_ICR]
enrichment.score.ICR[1,Sample_high_ICR]
enrichment.score.ICR.Bindea[25, Sample_high_ICR]
enrichment.score.ICR.ICR.Bindea[25, Sample_high_ICR]

enrichment.z.score.ICR[1,Sample_high_ICR]
enrichment.z.score.ICR[1,Sample_low_ICR]

enrichment.z.score.ICR.B[2,Sample_high_ICR]
enrichment.z.score.ICR.B[2,Sample_low_ICR]

enrichment.z.score.ICR.Bindea[25,Sample_high_ICR]
enrichment.z.score.ICR.Bindea[25,Sample_low_ICR]

enrichment.z.score.ICR.ICR.Bindea[25,Sample_high_ICR]
enrichment.z.score.ICR.ICR.Bindea[25,Sample_low_ICR]

enrichment.z.score.ICR[1,Sample_high_ICR] / enrichment.z.score.ICR[1,Sample_low_ICR]
enrichment.z.score.ICR.B[2,Sample_high_ICR] / enrichment.z.score.ICR.B[2,Sample_low_ICR]
enrichment.z.score.ICR.Bindea[25,Sample_high_ICR] / enrichment.z.score.ICR.Bindea[25,Sample_low_ICR]


Hallmark.enrichment.score = gsva(Expression.data, Gene_Lists, method="ssgsea")
Hallmark.enrichment.z.score = Hallmark.enrichment.score 
for(j in 1: nrow(Hallmark.enrichment.z.score))  {
  Hallmark.enrichment.z.score[j,] = (Hallmark.enrichment.score[j,]-mean(Hallmark.enrichment.score[j,]))/sd(Hallmark.enrichment.score[j,]) # z-score the enrichment matrix
}

Hallmark.enrichment.z.score[58, Sample_high_ICR]
Hallmark.enrichment.z.score[58, Sample_low_ICR]


