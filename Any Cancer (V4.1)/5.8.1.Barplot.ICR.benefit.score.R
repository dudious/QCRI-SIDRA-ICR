####################################################################
###
### This Script calculates 
### 
### Input data:
### ("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/")
### Output data are saved as Rdata file:
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
ipak(required.packages)

# Set Parameters
Pathways = "ALL"
CancerTYPES = "ALL"
Pathway_skip = ""                                                                                                      # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
order = "Enabled_disabled"
Log_file = paste0("./1_Log_Files/", download.method ,"/5.8.1.Pancancer_Survival_Analysis/5.8.1.Pancancer_Survival_Analysis_Log_File_",                              # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
Cutoff_HR = 1
ICR_benefit_score_cutoff = 0

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefitv1.Rdata"))
load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", "Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))

ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])

# Create directories
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Bar_plots_Oncogenic_clusters/"), showWarnings = FALSE)

if(order == "Enabled_disabled"){
  Cancer_order = c(ICR_enabled_cancers, ICR_disabled_cancers)
  Hallmark_and_ICR_cluster_assignment_allcancers$Cancer = factor(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer, levels = Cancer_order)
}

Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_cluster[which(Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_score < ICR_benefit_score_cutoff)] = "ICR_benefit_score_Low"
Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_cluster[which(Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_score > ICR_benefit_score_cutoff)] = "ICR_benefit_score_High"

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Bar_plots_Oncogenic_clusters/ICR_benefit_v1"), showWarnings = FALSE)
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Bar_plots_Oncogenic_clusters/ICR_benefit_v1/", "ICR_benefit_v1_",
           "_", order, "_Pancancer.png"),res=600,height=6,width=20,unit="in")

plot = ggplot(Hallmark_and_ICR_cluster_assignment_allcancers, aes(Cancer)) + geom_bar(aes(fill=ICR_benefit_cluster)) + labs(fill = "ICR_benefit_cluster") +
  theme(text = element_text(size = 15, colour = "black"))

plot(plot)
dev.off()
