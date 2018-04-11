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
ColsideLabels = c("ICR clusters", "Cancers", "ICR Enabled/Disabled")
Legend = c("ICR Low","ICR Med","ICR High")
Deconvolution_matrix = "Bindea.enrichment.score"                                                                 # "Bindea.enrichment.score", "bindea_patrick.enrichment.score"
Cutoff_HR = 1

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

load(paste0("./4_Analysis/", download.method, "/ACC/Signature_Enrichment/GSEA_ACC_Bindea_xCell_HallmarkPathways.Rdata"))
Enrichment_score_all = get(Deconvolution_matrix)

for (i in 2:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer == "LAML"){next}
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, "_Bindea_xCell_HallmarkPathways.Rdata"))
  Enrichment_score_all = cbind(Enrichment_score_all, get(Deconvolution_matrix))
}
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer")]

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]
annotation$ICR_ED = NA
annotation$ICR_ED[which(annotation$Cancer %in% ICR_enabled_cancers)] = "ICR_enabled"
annotation$ICR_ED[which(annotation$Cancer %in% ICR_disabled_cancers)] = "ICR_disabled"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_enabled"] = "orange"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_disabled"] = "purple"

# Bindea classification
annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"

annotation$Cancer.col = Cancer_color_table$color[match(annotation$Cancer, Cancer_color_table$Group.1)]

All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Cancer_order = as.character(All_survival_analysis_data$Cancertype[-which(All_survival_analysis_data$Cancertype == "LAML")])
ICR_order = c("ICR Low","ICR Medium","ICR High")

annotation = annotation[order(match(annotation$HML_cluster,ICR_order), match(annotation$Cancer, Cancer_order)),]

annotation.blot = as.matrix(annotation[,c("HML_cluster.col","Cancer.col", "ICR_ED.col"), drop = FALSE])
#annotation.blot = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
#Expression.data = Expression.data[colnames(annotation.blot),]

Enrichment_score_all.z.score = Enrichment_score_all 
for(j in 1: nrow(Enrichment_score_all.z.score))  {
  Enrichment_score_all.z.score[j,] = (Enrichment_score_all[j,]-mean(Enrichment_score_all[j,]))/sd(Enrichment_score_all[j,]) # z-score the enrichment matrix
}
Enrichment_score_all.z.score = Enrichment_score_all.z.score[,rownames(annotation.blot)]

## Determine order rows/signatures
to_aggregate = t(Enrichment_score_all.z.score)
to_aggregate = as.data.frame(to_aggregate)
to_aggregate$ICR_cluster = ICR_cluster_assignment_allcancers$HML_cluster[match(rownames(to_aggregate), rownames(ICR_cluster_assignment_allcancers))]
mean_ES_all = aggregate(.~ ICR_cluster, data = to_aggregate, FUN = mean)
mean_ES_all = t(mean_ES_all)
colnames(mean_ES_all) = mean_ES_all[1,]
mean_ES_all = mean_ES_all[-1,]
mode(mean_ES_all) = "numeric"
mean_ES_all = as.data.frame(mean_ES_all)
mean_ES_all$DeltaHL = c(mean_ES_all$`ICR High` - mean_ES_all$`ICR Low`)
mean_ES_all = mean_ES_all[order(mean_ES_all$DeltaHL, decreasing = TRUE),]
cell_order = rownames(mean_ES_all)

## Plot prep
Enrichment_score_all.z.score = Enrichment_score_all.z.score[c(cell_order),]
Legend = Cancer_order
Legend_colors = c(Cancer_color_table$color[order(match(Cancer_color_table$Group.1, Cancer_order))])

Legend2 = c("ICR Low","ICR Med","ICR High")
Legend_colors2 = c("blue","green","red")

Legend3 = c("ICR enabled", "ICR disabled")
Legend_colors3 = c("orange", "purple")

### Plotting
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/", Deconvolution_matrix,"Test_Cancers_reordered_Heatmap_RNASeq_Pancancer_HML.png"), res = 600, height = 15, width = 15, unit = "in")
heatmap.3((as.matrix(Enrichment_score_all.z.score)),
          main= paste0("Pancancer enrichment scores \nssGSEA/ Bindea"),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.3,
          cexRow = 1.3,
          margins = c(13, 30))

title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                        "Oncogenic pathway enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
legend("bottomleft",legend = Legend,
       col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("topright", legend = Legend2,
       col = Legend_colors2, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("top", legend = Legend3,
       col = Legend_colors3, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
dev.off()

if(Deconvolution_matrix == "Hallmark.enrichment.score"){
  dir.create("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment")
  Hallmark_enrichment_score_all = Enrichment_score_all
  Hallmark_enrichment_z_score_all = Enrichment_score_all.z.score
  save(Hallmark_enrichment_score_all, Hallmark_enrichment_z_score_all, file = 
         "./4_Analysis/", download.method,"/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark.Pancancer.Rdata")
}

