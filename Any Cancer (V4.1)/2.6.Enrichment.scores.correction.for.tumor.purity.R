################################################################
###
### Correction of ES scores for tumor purity
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize")
ibiopak("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,"/2.6.Correction.for.tumor.purity/2.6.Correction.for.tumor.purity_", # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark.Pancancer.Rdata"))
Leuk_estimate = read.csv(paste0("./3_DataProcessing/External/Leuk.infil.data.clean.csv"), stringsAsFactors = FALSE)
load(paste0("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

Hallmark_enrichment_score_for_correction = as.data.frame(t(Hallmark_enrichment_score_all))
Hallmark_enrichment_score_for_correction$ICR_score = ICR_cluster_assignment_allcancers$ICRscore[match(rownames(Hallmark_enrichment_score_for_correction),
                                                                                                      rownames(ICR_cluster_assignment_allcancers))]
rownames(Hallmark_enrichment_score_for_correction) = substring(rownames(Hallmark_enrichment_score_for_correction), 1, 12)


Leuk_estimate$SampleID = substring(Leuk_estimate$SampleID, 1, 12)

Hallmark_enrichment_score_for_correction$Leuk_estimate = Leuk_estimate$Leuk.Estimate[match(rownames(Hallmark_enrichment_score_for_correction), Leuk_estimate$SampleID)]
Hallmark_enrichment_score_for_correction = Hallmark_enrichment_score_for_correction[-which(is.na(Hallmark_enrichment_score_for_correction$Leuk_estimate)),]

dim(Hallmark_enrichment_score_for_correction)

Hallmark_ES_all_Leuko_corrected = Hallmark_enrichment_score_for_correction

i=2
for(i in 1:nrow(Hallmark_ES_all_Leuko_corrected)){
  Hallmark_ES_all_Leuko_corrected[i,] = Hallmark_ES_all_Leuko_corrected[i,] / Hallmark_ES_all_Leuko_corrected[i,ncol(Hallmark_ES_all_Leuko_corrected)] 
}

Hallmark_ES_all_Leuko_corrected$Leuk_estimate = NULL
Hallmark_ES_all_Leuko_corrected$ICR_score_not_corrected = ICR_cluster_assignment_allcancers$ICRscore[match(rownames(Hallmark_enrichment_score_for_correction),
                                                                                                           substring(rownames(ICR_cluster_assignment_allcancers), 1, 12))]
Hallmark_ES_all_Leuko_corrected = t(Hallmark_ES_all_Leuko_corrected)

save(Hallmark_ES_all_Leuko_corrected, 
     file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark.Pancancer.Leukocyte.Estimate.Corrected_v2.Rdata"))

