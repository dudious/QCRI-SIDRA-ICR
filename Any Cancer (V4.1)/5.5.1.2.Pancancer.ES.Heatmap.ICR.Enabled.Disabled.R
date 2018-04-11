#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
### Manual adjustment of min and max
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
Log_file = paste0("./1_Log_Files/", download.method ,"/5.5.1.Pancancer_Bindea_Heatmap/5.5.1.Pancancer.Bindea.Heatmap_", # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("Cancers", "ICR_Enabled/Disabled")
Legend = c("ICR Low","ICR Med","ICR High")
Deconvolution_matrix = "Hallmark.enrichment.score"                                                                 # "Bindea.enrichment.score", "bindea_patrick.enrichment.score"
Cutoff_HR = 1
Tumor_purity_correction = "No_correction"
pathway_filter = "significant_and_inverse_pathways"
include_proliferation = "include_proliferation"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/5.5.1.Pancancer_Bindea_Heatmap"), showWarnings = FALSE)

cat("This is a log file for creating heatmap bindea enrichment plots",
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    file = Log_file,
    append = FALSE, sep= "\n")

N.sets = length(CancerTYPES)

load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", "Survival_analysis_High_vs_Low_Groups", 
            "HML_classification", ".Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]
load("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata")

if(Tumor_purity_correction == "No_correction" & include_proliferation == "include_proliferation"){
  load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark&Proliferation.Pancancer.Rdata"))
  Enrichment_score_all = Hallmark_enrichment_score_all
}

if(Tumor_purity_correction == "Leukocyte_estimate_correction"){
  load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark.Pancancer.Leukocyte.Estimate.Corrected_v2.Rdata"))
  Enrichment_score_all = Hallmark_ES_all_Leuko_corrected
  rownames(ICR_cluster_assignment_allcancers) = substring(rownames(ICR_cluster_assignment_allcancers), 1, 12)
}

###
#Leuk_estimate = read.csv("./3_DataProcessing/External/Leuk.infil.data.clean.csv", stringsAsFactors = FALSE)
#ICR_cluster_assignment_allcancers$Leuk_estimate = Leuk_estimate$Leuk.Estimate[match(rownames(ICR_cluster_assignment_allcancers), substring(Leuk_estimate$SampleID, 1, 12))]
#ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[-which(is.na(ICR_cluster_assignment_allcancers$Leuk_estimate)),]
#ICR_cluster_assignment_allcancers$Leuk_estimate[which(ICR_cluster_assignment_allcancers$Leuk_estimate >= 0.2)] = "High"
#ICR_cluster_assignment_allcancers$Leuk_estimate[which(ICR_cluster_assignment_allcancers$Leuk_estimate < 0.2)] = "Low"
#annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer", "Leuk_estimate")]
#annotation$LE.col[which(annotation$Leuk_estimate == "Low")] = "lightblue"
#annotation$LE.col[which(annotation$Leuk_estimate == "High")] = "pink"
###


annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer")]
annotation$ICR_ED = NA
annotation$ICR_ED[which(annotation$Cancer %in% ICR_enabled_cancers)] = "ICR_enabled"
annotation$ICR_ED[which(annotation$Cancer %in% ICR_disabled_cancers)] = "ICR_disabled"

# Bindea classification
annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"

annotation$Cancer.col = Cancer_color_table$color[match(annotation$Cancer, Cancer_color_table$Group.1)]

annotation$ICR_ED.col[annotation$ICR_ED == "ICR_enabled"] = "orange"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_disabled"] = "purple"

All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Cancer_order = as.character(All_survival_analysis_data$Cancertype[-which(All_survival_analysis_data$Cancertype == "LAML")])
ICR_order = c("ICR Low","ICR Medium","ICR High")
ICR_ED_order = c("ICR_enabled", "ICR_disabled")
##
#Leuk_est_order = c("Low", "High")
#annotation = annotation[order(match(annotation$ICR_ED, ICR_ED_order), match(annotation$Leuk_estimate, Leuk_est_order)),]
#annotation = annotation[order(match(annotation$HML_cluster, ICR_order), match(annotation$Cancer, Cancer_order)),]

#annotation = annotation[order(match(annotation$HML_cluster, ICR_order)),]
annotation = annotation[order(match(annotation$ICR_ED, ICR_ED_order), match(annotation$Cancer, Cancer_order)),]

annotation.blot = as.matrix(annotation[,c("Cancer.col", "ICR_ED.col"), drop = FALSE])
annotation.blot = annotation.blot[which(rownames(annotation.blot) %in% colnames(Enrichment_score_all)),]
#annotation.blot = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
#Expression.data = Expression.data[colnames(annotation.blot),]

Enrichment_score_all.z.score = Enrichment_score_all 
j=1
for(j in 1: nrow(Enrichment_score_all.z.score))  {
  Enrichment_score_all.z.score[j,] = (Enrichment_score_all[j,]-mean(Enrichment_score_all[j,]))/sd(Enrichment_score_all[j,]) # z-score the enrichment matrix
}
Enrichment_score_all.z.score = Enrichment_score_all.z.score[,rownames(annotation.blot)]

ICR_cluster_assignment_allcancers$ICR_ED = NA
ICR_cluster_assignment_allcancers$ICR_ED[ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers] = "ICR_enabled"
ICR_cluster_assignment_allcancers$ICR_ED[ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers] = "ICR_disabled"

## Determine order rows/signatures
to_aggregate = t(Enrichment_score_all.z.score)
to_aggregate = as.data.frame(to_aggregate)
to_aggregate$ICR_ED = ICR_cluster_assignment_allcancers$ICR_ED[match(rownames(to_aggregate), rownames(ICR_cluster_assignment_allcancers))]
mean_ES_all = aggregate(.~ ICR_ED, data = to_aggregate, FUN = mean)
mean_ES_all = t(mean_ES_all)
colnames(mean_ES_all) = mean_ES_all[1,]
mean_ES_all = mean_ES_all[-1,]
mode(mean_ES_all) = "numeric"
mean_ES_all = as.data.frame(mean_ES_all)
mean_ES_all$DeltaED = c(mean_ES_all$ICR_enabled - mean_ES_all$ICR_disabled)
mean_ES_all = mean_ES_all[order(mean_ES_all$DeltaED, decreasing = TRUE),]
cell_order = rownames(mean_ES_all)
Enrichment_score_all.z.score = Enrichment_score_all.z.score[c(cell_order),]

if(pathway_filter == "significant_and_inverse_pathways"){
  load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/ICR_All_Pathway_High_Survival_analysis_High_vs_Low_Groups_Oncogenic_pathways.Rdata"))
  All_survival_analysis_data_ONCO_High = All_survival_analysis_data_ONCO
  rm(All_survival_analysis_data_ONCO)
  load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/ICR_All_Pathway_Low_Survival_analysis_High_vs_Low_Groups_Oncogenic_pathways.Rdata"))
  All_survival_analysis_data_ONCO_Low = All_survival_analysis_data_ONCO
  rm(All_survival_analysis_data_ONCO)
  
  merged = merge(All_survival_analysis_data_ONCO_High, All_survival_analysis_data_ONCO_Low, by = "Oncogenic_Pathway",
                 suffixes = c("_High", "_Low"))
  
  merged$score_contribution = NA
  merged$score_contribution[which(merged$HR_High >1 & merged$HR_Low <1 & merged$p_value_High < 0.05 & merged$p_value_Low < 0.05)] = "enabling"
  merged$score_contribution[which(merged$HR_High <1 & merged$HR_Low >1 & merged$p_value_High < 0.05 & merged$p_value_Low < 0.05)] = "disabling"
  
  enabling_pathways = gsub("_cluster_Pancancer", "", as.character(merged$Oncogenic_Pathway[which(merged$score_contribution == "enabling")]))
  enabling_pathways[which(enabling_pathways == "Proliferation")] = c("[LM] Proliferation")
  disabling_pathways = gsub("_cluster_Pancancer", "", as.character(merged$Oncogenic_Pathway[which(merged$score_contribution == "disabling")]))
  
  Enrichment_score_all.z.score = Enrichment_score_all.z.score[c(enabling_pathways, disabling_pathways),]
}

## Plot prep
Legend = Cancer_order
rownames(Cancer_color_table) = Cancer_color_table$Group.1
Cancer_color_table = Cancer_color_table[Cancer_order,]
Legend_colors = c(Cancer_color_table$color)

#Legend2 = c("ICR Low","ICR Med","ICR High")
#Legend_colors2 = c("blue","green","red")

Legend3 = c("ICR enabled", "ICR disabled")
Legend_colors3 = c("orange", "purple")

### Plotting
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/", Deconvolution_matrix, "_", Tumor_purity_correction,
           "_", pathway_filter,"_ICR_E_D_seperated_v2", "_Heatmap_RNASeq_Pancancer_HML.png"), 
    res = 600, height = 15, width = 15, unit = "in")
heatmap.3((as.matrix(Enrichment_score_all.z.score)),
          main= paste0("Pancancer enrichment scores \nssGSEA/Oncogenic pathway signatures \n", Tumor_purity_correction),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          #Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.25,
          cexRow = 1.3,
          margins = c(13, 30))

title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                        "Oncogenic pathway enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
legend("bottomleft",legend = Legend,
       col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
#legend("top", legend = Legend2,
       #col = Legend_colors2, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("topright", legend = Legend3,
       col = Legend_colors3, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)

dev.off()

## Color Key

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Color_Key_", Deconvolution_matrix, "_", Tumor_purity_correction,
           "_", pathway_filter,"_ICR_E_D_seperated_v2", "_Heatmap_RNASeq_Pancancer_HML.png"), 
    res = 600, height = 7, width = 7, unit = "in")
heatmap.3((as.matrix(Enrichment_score_all.z.score)),
          main= paste0("Pancancer enrichment scores \nssGSEA/Oncogenic pathway signatures \n", Tumor_purity_correction),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.25,
          cexRow = 1.3,
          margins = c(13, 30))
dev.off()
