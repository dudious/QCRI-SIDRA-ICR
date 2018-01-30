#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
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
ibiopak("ComplexHeatmap")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = c("ALL")                                                                                                   # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/5.5_Pancancer_Bindea_xCell_Dot_Heatmap/5.5_Pancancer_Bindea_xCell_Dot_Heatmap", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
Deconvolution_method = "Bindea"                                                                                           # c("Bindea", "bindea_patrick", "xCell")
Plot_type = "Mean_all"                                                                                                  # "Mean_all" or "ICR_High_vs_Low"     
enrichment.score.type = ".enrichment.score"
scaling = "z_score"

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
dir.create(paste0("./5_Figures/Pancancer_plots/Bindea_enrichment_pancancer_Dot_Heatmap"), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/5.5_Pancancer_Bindea_xCell_Dot_Heatmap"), showWarnings = FALSE)

cat("This is a log file for creating dot heatmap bindea enrichment plots",
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

start.time = Sys.time ()

load("./4_Analysis/TCGA_Assembler/ACC/Signature_Enrichment/GSEA_ACC_Bindea_xCell_HallmarkPathways.Rdata")
if(Deconvolution_method == "xCell"){
  Selected_Deconvolution_method_score = xCell.matrix
}else{Selected_Deconvolution_method_score = get(paste0(Deconvolution_method,enrichment.score.type))}

Cancers_included = CancerTYPES[-which(CancerTYPES =="LAML")]
Deconvolution_score_per_ICR_df = data.frame(matrix(nrow = length(rownames(Selected_Deconvolution_method_score)), ncol = length(Cancers_included) * 4), row.names = rownames(Selected_Deconvolution_method_score))
ICR_Highs = paste(Cancers_included, "_ICR_High", sep = "")
ICR_Lows = paste(Cancers_included, "_ICR_Low", sep = "")
ICR_delta = paste(Cancers_included, "_ICR_High_vs_Low", sep = "")
Mean_all = paste(Cancers_included, "_Mean_all", sep = "")
colnames(Deconvolution_score_per_ICR_df) = c(ICR_Highs, ICR_Lows, ICR_delta, Mean_all)

i = 2
for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))){next}
  load(paste0("./4_Analysis/", download.method, "/", Cancer,"/Signature_Enrichment/GSEA_", Cancer, "_Bindea_xCell_HallmarkPathways.Rdata"))
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/", 
              Cancer, "_ICR_cluster_assignment_k2-6.Rdata"))
  
  if(Deconvolution_method == "xCell"){
    Selected_Deconvolution_method_score = xCell.matrix
  }else{Selected_Deconvolution_method_score = get(paste0(Deconvolution_method,enrichment.score.type))}
  ICR_high_patients = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR High")]
  ICR_low_patients = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR Low")]
  Selected_Deconvolution_method_score = as.data.frame(Selected_Deconvolution_method_score)
  Selected_Deconvolution_method_score$ICR_High_mean = rowMeans(Selected_Deconvolution_method_score[, which(colnames(Selected_Deconvolution_method_score) %in% ICR_high_patients)])
  Selected_Deconvolution_method_score$ICR_Low_mean = rowMeans(Selected_Deconvolution_method_score[, which(colnames(Selected_Deconvolution_method_score) %in% ICR_low_patients)])
  Selected_Deconvolution_method_score$Mean_all = rowMeans(Selected_Deconvolution_method_score)
  
  ICR_High_Cancer = paste0(Cancer, "_ICR_High")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_High_Cancer)] = Selected_Deconvolution_method_score$ICR_High_mean[match(rownames(Selected_Deconvolution_method_score),
                                                                                                                                                           rownames(Deconvolution_score_per_ICR_df))]
  ICR_Low_Cancer = paste0(Cancer, "_ICR_Low")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_Low_Cancer)] = Selected_Deconvolution_method_score$ICR_Low_mean[match(rownames(Selected_Deconvolution_method_score),
                                                                                                                                                          rownames(Deconvolution_score_per_ICR_df))]
  ICR_Delta_Cancer = paste0(Cancer, "_ICR_High_vs_Low")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_Delta_Cancer)] = Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_High_Cancer)] -
    Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_Low_Cancer)]
  
  Mean_all_Cancer = paste0(Cancer, "_Mean_all")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == Mean_all_Cancer)] = Selected_Deconvolution_method_score$Mean_all[match(rownames(Selected_Deconvolution_method_score),
                                                                                                                                                              rownames(Deconvolution_score_per_ICR_df))]
}

if(Plot_type == "ICR_High_vs_Low"){
  Deconvolution_score_for_plot = Deconvolution_score_per_ICR_df[, grep(pattern = "_ICR_High_vs_Low", colnames(Deconvolution_score_per_ICR_df))]
  colnames(Deconvolution_score_for_plot) = gsub("_ICR_High_vs_Low","", colnames(Deconvolution_score_for_plot))
}

if(Plot_type == "Mean_all"){
  Deconvolution_score_for_plot = Deconvolution_score_per_ICR_df[, grep(pattern = "_Mean_all", colnames(Deconvolution_score_per_ICR_df))]
  colnames(Deconvolution_score_for_plot) = gsub("_Mean_all","", colnames(Deconvolution_score_for_plot))
}

Deconvolution_score_for_plot = Deconvolution_score_for_plot[order(rowMeans(Deconvolution_score_for_plot), decreasing = TRUE),]
Deconvolution_score_for_plot = data.matrix(Deconvolution_score_for_plot)

Deconvolution_score_for_plot_z_scored = Deconvolution_score_for_plot
for(j in 1: nrow(Deconvolution_score_for_plot))  {
  Deconvolution_score_for_plot_z_scored[j,] = (Deconvolution_score_for_plot[j,]-mean(Deconvolution_score_for_plot[j,]))/sd(Deconvolution_score_for_plot[j,]) # z-score the enrichment matrix
}

if(scaling == "z_score"){
  Deconvolution_score_for_plot = Deconvolution_score_for_plot_z_scored
}

min = -5
max = 5

col_fun = circlize::colorRamp2(c(min, 0, max), c("blue", "white", "red"))

dir.create("./5_Figures", showWarnings = FALSE)
dir.create("./5_Figures/Pancancer_plots/", showWarnings = FALSE)
dir.create("./5_Figures/Pancancer_plots/DotHeatmaps", showWarnings = FALSE)
dir.create(paste("./5_Figures/Pancancer_plots/DotHeatmaps/", download.method), showWarnings = FALSE)

#DOT HEATMAP
#png(paste0("./5_Figures/Pancancer_plots/DotHeatmaps/", download.method, "/", Bindea, "_DotHeatmap_ICR_High_vs_Low.png"), res=600,height=10,width=10,unit="in")
dev.new()
Heatmap(Deconvolution_score_for_plot,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(type = "none"),
        column_names_side = "top",
        column_names_gp = gpar(fontsize =20),
        row_names_gp = gpar(fontsize = 20),
        col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
        heatmap_legend_param = list(at= c(min, 0, max), labels = c(min, "0", max),                             ## Make sure this is the same as col_fun!!
          legend_direction = "horizontal", legend_height = unit(5, "in"), title_position = "lefttop"),
        #top_annotation = ha_column,
        name = "Mean enrichment score (z-scored by row)",
        #gap = unit(5, "mm"),
        row_title_gp = gpar(fontsize = 20),
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        grid.circle(x = x, y = y, r = 0.02,gp = gpar(fill = col_fun(Deconvolution_score_for_plot[i, j]), col = NA))
        }
)

#title(main = paste0(Bindea, " signature expression in ICR High tumors compared to ICR Low tumors"))
#dev.off()
