#################################################################
###
###
### Data input:
###
### Output :
### 
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

required.packages <- c("gtools", "circlize")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"
Pathway_skip = ""                                                                                                        
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,
                  "/5.7.ICR.benefit.score/5.7.Pancancer.benefit.score_",          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/Hallmark_and_ICR_cluster_assignment_allcancers.Rdata"))

load(paste0("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Survival_Analysis/ICR_All_Pathway_High_Survival_analysis_High_vs_Low_Groups_Oncogenic_pathways.Rdata"))
All_survival_analysis_data_ONCO_High = All_survival_analysis_data_ONCO
rm(All_survival_analysis_data_ONCO)
load(paste0("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Survival_Analysis/ICR_All_Pathway_Low_Survival_analysis_High_vs_Low_Groups_Oncogenic_pathways.Rdata"))
All_survival_analysis_data_ONCO_Low = All_survival_analysis_data_ONCO
rm(All_survival_analysis_data_ONCO)

merged = merge(All_survival_analysis_data_ONCO_High, All_survival_analysis_data_ONCO_Low, by = "Oncogenic_Pathway",
               suffixes = c("_High", "_Low"))

merged$score_contribution = NA
merged$score_contribution[which(merged$HR_High >1 & merged$HR_Low <1 & merged$p_value_High < 0.05 & merged$p_value_Low < 0.05)] = "enabling"
merged$score_contribution[which(merged$HR_High <1 & merged$HR_Low >1 & merged$p_value_High < 0.05 & merged$p_value_Low < 0.05)] = "disabling"

## Delete redundant scores here by subsetting

enabling_pathways = as.character(merged$Oncogenic_Pathway[which(merged$score_contribution == "enabling")])
disabling_pathways = as.character(merged$Oncogenic_Pathway[which(merged$score_contribution == "disabling")])

score_calculation_df = Hallmark_and_ICR_cluster_assignment_allcancers[, which(colnames(Hallmark_and_ICR_cluster_assignment_allcancers) %in% enabling_pathways |
                                                                                colnames(Hallmark_and_ICR_cluster_assignment_allcancers) %in% disabling_pathways)]
for (i in 1:length(enabling_pathways)){
  Pathway_cluster = enabling_pathways[i]
  Pathway_name = gsub("_cluster_Pancancer", "", Pathway_cluster)
  score_calculation_df[,Pathway_cluster][which(score_calculation_df[,Pathway_cluster] == paste0(Pathway_name, " High"))] = 1
  score_calculation_df[,Pathway_cluster][which(score_calculation_df[,Pathway_cluster] == paste0(Pathway_name, " Low"))] = -1
}
for (j in 1:length(disabling_pathways)){
  Pathway_cluster = disabling_pathways[j]
  Pathway_name = gsub("_cluster_Pancancer", "", Pathway_cluster)
  score_calculation_df[,Pathway_cluster][which(score_calculation_df[,Pathway_cluster] == paste0(Pathway_name, " High"))] = -1
  score_calculation_df[,Pathway_cluster][which(score_calculation_df[,Pathway_cluster] == paste0(Pathway_name, " Low"))] = 1
}
score_calculation_df = as.matrix(score_calculation_df)
mode(score_calculation_df) = "numeric"
score_calculation_df = as.data.frame(score_calculation_df) 
score_calculation_df$score = rowSums(score_calculation_df)

Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_score = score_calculation_df$score[match(rownames(Hallmark_and_ICR_cluster_assignment_allcancers),
                                                                                                        rownames(score_calculation_df))]
save(Hallmark_and_ICR_cluster_assignment_allcancers, file = paste0("./4_Analysis/", download.method ,"/Pan_Cancer/Clustering/Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefitv1.Rdata"))                                                                                                       
