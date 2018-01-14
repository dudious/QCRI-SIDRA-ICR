#################################################################
###
### Make oncoPrint plots of hallmark pathways
### Input files:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata"
### "./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
### Cancer, "_ICR_cluster_assignment_k2-6.Rdata"
### Output files:
### "./5_Figures/OncoPrints/Hallmark_OncoPrints/", download.method, "/Hallmark_OncoPrint_", Cancer, ".png"
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

required.bioconductor.packages = c("ComplexHeatmap", "circlize", "gridExtra", "rlist")
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/oncoPrint.R"))                                                                        # source adapted version of oncoPrint version (Other functions of ComplexHeatmap packages are still required still required)
source(paste0(code_path, "R tools/heatmap.3.R"))

# Set Parameters
CancerTYPES = c("BRCA", "COAD")                                                                                 # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
pathway_selection = "intersection"                                                                                             # Specify whether you want to combine pathways that overlap between all cancers that are used as input (intersection)
                                                                                                                        # or if you want include all pathways that appear in either one of the cancers used as input (union). 
download.method = "TCGA_Assembler"
Log_file = paste0("./1_Log_Files/5.4_OncoPrint_Combining_Cancers_RNASeq/3.12_OncoPrint_RNASeq_Log_File_",               # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
subset = "COR_COEF"                                                                                                      # Options: "ALL_SIG" or "INV_COR_SIG" or "POS_COR_SIG" or "COR_COEF"
cor_cutoff = 0
ICR_medium_excluded = "ICR_medium_excluded"                                                                              # Options: "ICR_medium_excluded" or "all_included"


# Load data
load("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/TCGA_Assembler/z_score_1.5_COR_COEF_0/All_Hallmark_pathways_Oncoprint_matrixesz_score_1.5_COR_COEF_0.Rdata")

if(subset == "ALL_SIG" | subset == "INV_COR_SIG" | subset == "POS_COR_SIG"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_only_significant.Rdata")
}
if(subset == "COR_COEF"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_irrespective_of_significance.Rdata")
}


# Combine datasets of interest
N.sets = length(CancerTYPES)

for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  assign(paste0("Complete_Oncoprint_", Cancer), get(grep(pattern = Cancer, ls(), value = TRUE)))
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  assign(paste0("Table_Cluster_Assignment_", Cancer), table_cluster_assignment)
}

# Combine oncoprint matrixes of interest
oncoprints_to_combine = grep(pattern = "Complete_Oncoprint_", ls(), value = TRUE)
oncoprints_to_combine_list = mget(oncoprints_to_combine)
oncoprint_combined = do.call("cbind", oncoprints_to_combine_list)

# Combine cluster assignment tables of interest
cluster_assignments_to_combine = grep(pattern = "Table_Cluster_Assignment_", ls(), value = TRUE)
cluster_assignments_to_combine_list = mget(cluster_assignments_to_combine)
cluster_assignments_to_combine_list = unname(cluster_assignments_to_combine_list)                                                           # Unname else the next step will generate new rownames
cluster_assignments_combined = do.call("rbind", cluster_assignments_to_combine_list)

# Combine oncoprints of interest with ICR cluster assignment
oncoprint_matrix = rbind(oncoprint_combined, cluster_assignments_combined$HML_cluster[match(colnames(oncoprint_combined), row.names(cluster_assignments_combined))])
rownames(oncoprint_matrix)[nrow(oncoprint_matrix)] = "ICR cluster"
oncoprint_matrix[oncoprint_matrix == "ICR Low" | oncoprint_matrix == "ICR Medium"] = NA
oncoprint_matrix[oncoprint_matrix == "ICR High"] = "UP"

oncoprint_matrix[is.na(oncoprint_matrix)] = ""

complete_matrix = oncoprint_matrix

# Define pathways to display in the plot
if(subset == "COR_COEF"){
  for (i in 1:N.sets){
    Cancer = CancerTYPES[i]
    assign(paste0("cor_pathways_", Cancer), row.names(pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < -cor_cutoff),]))
  }
  cor_pathways_to_combine = grep(pattern = "cor_pathways_", ls(), value = TRUE)
  cor_pathways_to_combine_list = mget(cor_pathways_to_combine)
  
  if(pathway_selection == "intersection"){
    cor_pathways = Reduce(intersect, cor_pathways_to_combine_list)
  }
  if(pathway_selection == "union"){
    cor_pathways = unique(Reduce(c, cor_pathways_to_combine_list))
  }
}

cor_pathways = c(cor_pathways, "ICR cluster")

# Filter to only plot the selected cor_pathways
complete_matrix = complete_matrix[which(rownames(complete_matrix) %in% cor_pathways),]

pathway_order = c("ICR cluster", cor_pathways[-which(cor_pathways == "ICR cluster")])

dir.create(paste0("./5_Figures/Pancancer_plots/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/OncoPrints"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_")), showWarnings = FALSE)

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = "grey", col = NA)),
  UP = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UP"], col = NA)))

col = c(UP = "red", DOWN = "blue")

if(ICR_medium_excluded == "all_included"){
  number_of_patients = nrow(cluster_assignments_combined)
}
if(ICR_medium_excluded == "ICR_medium_excluded"){
  number_of_patients = nrow(cluster_assignments_combined[cluster_assignments_combined$HML_cluster != "ICR Medium",])
}

# Filter out Medium Tumors if needed (specify in parameters in top of this script)
if(ICR_medium_excluded == "ICR_medium_excluded"){
  excluded_patients = rownames(cluster_assignments_combined)[cluster_assignments_combined$HML_cluster == "ICR Medium"]
  matrix_oncoprint_input = complete_matrix[,-which(colnames(complete_matrix) %in% excluded_patients)]
}

png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_", subset, "_", cor_cutoff, "_", ICR_medium_excluded ,".png"),res=600,height= 9,width= 18,unit="in")
oncoPrint(matrix_oncoprint_input,
          alter_fun = alter_fun,
          col = col,
          heatmap_legend_param = list(title= "", at = c("UP"),
                                      labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
          column_title = paste0("OncoPrint: ", paste(CancerTYPES, collapse = " "), " RNASeq expression", "\n N patients = ", number_of_patients),
          row_order = pathway_order,
          show_row_barplot = FALSE,
          row_barplot_width = unit(7, "in"),
          top_annotation = NULL,
          width = unit(8, "in"),
          row_title = "Hallmark pathways",
          column_title_gp = gpar(fontface = 2, fontsize = 16))
dev.off()
