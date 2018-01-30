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
#Group_yellow = c("UCS", "OV", "SARC", "THYM", "GBM", "LGG")
Group_red = c("UCEC", "MESO", "LUAD", "LUSC", "LIHC", "KIRC", "STAD", "HNSC", "ESCA", "SKCM",
              "CHOL", "CESC")
CancerTYPES = Group_red                                                                                                 # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
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
cor_cutoff = 0.1
ICR_medium_excluded = "ICR_medium_excluded"                                                                              # Options: "ICR_medium_excluded" or "all_included"
IPA_excluded = "IPA_excluded"


# Load data
load(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/TCGA_Assembler/z_score_1.5_COR_COEF_", cor_cutoff,
            "/All_Hallmark_pathways_Oncoprint_matrixesz_score_1.5_COR_COEF_", cor_cutoff, ".Rdata"))

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

if(subset == "ALL_SIG" | subset == "INV_COR_SIG" | subset == "POS_COR_SIG"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_only_significant.Rdata")
}
if(subset == "COR_COEF"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_irrespective_of_significance.Rdata")
}

N.all.cancers = length(TCGA.cancersets)

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

N.cancers = length(cluster_assignments_to_combine)

i=2
for (i in 1:N.cancers){
  Cancer = CancerTYPES[i]
  tablecluster_temp = get(grep(pattern = paste0("Table_Cluster_Assignment_", Cancer), ls(), value = TRUE))
  tablecluster_temp$HML_cluster[tablecluster_temp$HML_cluster == "ICR High"] = Cancer
  assign(paste0("Table_Cluster_Assignment_", Cancer), tablecluster_temp)
}
cluster_assignments_to_combine = grep(pattern = "Table_Cluster_Assignment_", ls(), value = TRUE)
cluster_assignments_to_combine_list = mget(cluster_assignments_to_combine)
cluster_assignments_to_combine_list = unname(cluster_assignments_to_combine_list)                                                           # Unname else the next step will generate new rownames
cluster_assignments_combined = do.call("rbind", cluster_assignments_to_combine_list)

# Combine oncoprints of interest with ICR cluster assignment
oncoprint_matrix = rbind(oncoprint_combined, cluster_assignments_combined$HML_cluster[match(colnames(oncoprint_combined), row.names(cluster_assignments_combined))])
rownames(oncoprint_matrix)[nrow(oncoprint_matrix)] = "ICR cluster"
oncoprint_matrix[oncoprint_matrix == "ICR Low" | oncoprint_matrix == "ICR Medium"] = NA

oncoprint_matrix[is.na(oncoprint_matrix)] = ""

complete_matrix = oncoprint_matrix

# Define pathways to display in the plot
if(subset == "COR_COEF"){
  for (i in 1:N.sets){
    Cancer = CancerTYPES[i]
    assign(paste0("cor_pathways_", Cancer), row.names(pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < cor_cutoff),]))
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

if(IPA_excluded == "IPA_excluded"){
  cor_pathways = cor_pathways[grep(pattern = "IPA", cor_pathways, invert = TRUE)]
}
if(length(cor_pathways) == 0){
  next
}

cor_pathways = c(cor_pathways, "ICR cluster")

# Filter to only plot the selected cor_pathways
complete_matrix = complete_matrix[which(rownames(complete_matrix) %in% cor_pathways),]

pathway_order = c("ICR cluster", cor_pathways[-which(cor_pathways == "ICR cluster")])

dir.create(paste0("./5_Figures/Pancancer_plots/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/OncoPrints"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/OncoPrints/Multicolor"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/OncoPrints/Multicolor/", paste(CancerTYPES, collapse = "_")), showWarnings = FALSE)

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = "grey", col = NA)),
  UP = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UP"], col = NA)),
  ACC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["ACC"], col = NA)),
  BLCA = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["BLCA"], col = NA)),
  BRCA = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["BRCA"], col = NA)),
  CESC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["CESC"], col = NA)),
  CHOL = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["CHOL"], col = NA)),
  COAD = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["COAD"], col = NA)),
  DLBC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["DLBC"], col = NA)),
  ESCA = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["ESCA"], col = NA)),
  GBM = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["GBM"], col = NA)),
  HNSC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["HNSC"], col = NA)),
  KICH = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["KICH"], col = NA)),
  KIRC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["KIRC"], col = NA)),
  KIRP = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["KIRP"], col = NA)),
  LGG = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["LGG"], col = NA)),
  LIHC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["LIHC"], col = NA)),
  LUAD = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["LUAD"], col = NA)),
  LUSC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["LUSC"], col = NA)),
  MESO = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["MESO"], col = NA)),
  OV = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["OV"], col = NA)),
  PAAD = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["PAAD"], col = NA)),
  PCPG = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["PCPG"], col = NA)),
  PRAD = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["PRAD"], col = NA)),
  READ = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["READ"], col = NA)),
  SARC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["SARC"], col = NA)),
  SKCM = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["SKCM"], col = NA)),
  STAD = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["STAD"], col = NA)),
  TGCT = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["TGCT"], col = NA)),
  THCA = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["THCA"], col = NA)),
  THYM = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["THYM"], col = NA)),
  UCEC = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UCEC"], col = NA)),
  UCS = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UCS"], col = NA)),
  UVM = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UVM"], col = NA)))

col = c(UP = "red", ACC = "#666666", BLCA = "#F12B7E", BRCA = "#D9D9D9", CESC = "#E5C494", CHOL = "#FC9998",
        COAD = "#D96002", DLBC = "#8DD3C7", ESCA = "#FC8D62", GBM = "#BC7FBC", HNSC = "#F1E1CC", KICH = "#387EB7",
        KIRC = "#34A02C", KIRP = "#B2DE68", LGG = "#BEB9DA", LIHC = "#F781BF", LUAD = "#7FB1D3", LUSC = "#DECAE4",
        MESO = "#FFFECC", OV = "#FDCDAC", PAAD = "#396BAF", PCPG = "#FFD82F", PRAD = "#E72A89", READ = "#FEB462",
        SARC = "#E4211E", SKCM = "#E6AB03", STAD = "#A5761D", TGCT = "#FB7F72", THCA = "#A6CEE3", THYM = "#FFED6F",
        UCEC = "#FFFD99", UCS = "#4DAE4B", UVM = "#B3B3B3")

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

png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", "Multicolor/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_", subset, "_", cor_cutoff, "_", ICR_medium_excluded, "_",
           IPA_excluded, ".png"),res=600,height= 9,width= 25,unit="in")
oncoPrint(matrix_oncoprint_input,
          alter_fun = alter_fun,
          col = col,
          #heatmap_legend_param = list(title= "", at = c("UP"),
          #labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
          heatmap_legend_param = list(title = "Upregulation in cancers: "),
          column_title = paste0("OncoPrint: ", paste(CancerTYPES, collapse = " "), " RNASeq expression", "\n N patients = ", number_of_patients,
                                "\n cor cutoff = ", cor_cutoff, " pathways included: ", pathway_selection, " of all cancers"),
          row_order = pathway_order,
          show_row_barplot = FALSE,
          row_barplot_width = unit(7, "in"),
          top_annotation = NULL,
          width = unit(8, "in"),
          row_title = "Hallmark pathways",
          column_title_gp = gpar(fontface = 2, fontsize = 16))
dev.off()
