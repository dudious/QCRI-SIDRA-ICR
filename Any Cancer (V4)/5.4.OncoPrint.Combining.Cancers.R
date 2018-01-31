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

required.bioconductor.packages = c("ComplexHeatmap", "circlize", "gridExtra", "rlist", "reshape2", "ggplot2")
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/oncoPrint.R"))                                                                        # source adapted version of oncoPrint version (Other functions of ComplexHeatmap packages are still required still required)
source(paste0(code_path, "R tools/heatmap.3.R"))

# Previously used parameters
Group_yellow = c("GBM", "LGG","OV", "SARC", "THYM","UCS")
Group_red = c("LUAD", "LUSC","MESO","STAD", "UCEC")
Group_blue = c("CESC", "ESCA", "HNSC","SKCM")
Group_green = c("BLCA", "COAD", "KICH", "KIRC", "KIRP",
                "LIHC", "PAAD", "PCPG", "READ", "THCA",
                "UVM")

Pathway_set_in_plot = c("[HM] KRAS SIGNALING DN", "[HM] HEDGEHOG SIGNALING", "[HM] NOTCH SIGNALING", "[HM] WNT BETA CATENIN SIGNALING",
                        "[TBI] BARRIER GENES", "[HM] GLYCOLYSIS", "[HM] ESTROGEN RESPONSE", "[TPW] Mesenchymal Transition",
                        "[HM] TGF BETA SIGNALING", "[HM] UV RESPONSE DN", "[TBI] MAPK UP GENES", "[HM] MITOTIC SPINDLE",
                        "[TBI] PHOSPHOLIPASE", "[TPW] NOS1 Signature", "[HM] E2F TARGETS", "[HM] G2M CHECKPOINT",
                        "[HM] OXIDATIVE PHOSPHORYLATION", "[HM] MYC TARGETS", "[HM] DNA REPAIR")

Pathway_set1 = c("BARRIER GENES", "[HM] GLYCOLYSIS", "[HM] ESTROGEN RESPONSE LATE", "[HM] ESTROGEN RESPONSE EARLY",
                 "[HM] KRAS SIGNALING DN", "[HM] WNT BETA CATENIN SIGNALING", "[HM] NOTCH SIGNALING", "[HM] HEDGEHOG SIGNALING")
Pathway_set2 = c("[HM] ESTROGEN RESPONSE EARLY",
                 "[HM] KRAS SIGNALING DN", "[HM] WNT BETA CATENIN SIGNALING", "[HM] NOTCH SIGNALING", "[HM] HEDGEHOG SIGNALING")



# Set Parameters
CancerTYPES = Group_yellow                                                                                                # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
pathway_selection = "intersection"                                                                                      # Specify whether you want to combine pathways that overlap between all cancers that are used as input (intersection)
                                                                                                                        # or if you want include all pathways that appear in either one of the cancers used as input (union). 
download.method = "TCGA_Assembler"
Log_file = paste0("./1_Log_Files/5.4_OncoPrint_Combining_Cancers_RNASeq/3.12_OncoPrint_RNASeq_Log_File_",               # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
subset = "<0.1_for_at_least_one_third_of_cancers"                                                                                   # Options: "COR_COEF" or c("[HM] WNT BETA CATENIN SIGNALING", "[HM] NOTCH SIGNALING"..etc)
minimum_number_of_cancers = ceiling(length(CancerTYPES)/3)
cor_cutoff = -0.1                                                                                                        # Cor-cutoff to define which pathways are selected for each cancer
ICR_medium_excluded = "ICR_medium_excluded"                                                                              # Options: "ICR_medium_excluded" or "all_included"
IPA_excluded = "IPA_excluded"

# Load data
load(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v5/TCGA_Assembler/z_score_1.5_COR_COEF_", cor_cutoff,
            "/All_Hallmark_pathways_Oncoprint_matrixesz_score_1.5_COR_COEF_", cor_cutoff, ".Rdata"))

if(subset == "ALL_SIG" | subset == "INV_COR_SIG" | subset == "POS_COR_SIG"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_only_significant_IPA_excluded.Rdata")
}
if(subset == "COR_COEF"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_irrespective_of_significance_IPA_excluded.Rdata")
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
}else{
  cor_pathways = subset
}

if(subset == "<0.1_for_at_least_one_third_of_cancers"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_irrespective_of_significance_IPA_excluded.Rdata")
  other_pathways = setdiff(rownames(pancancer_Hallmark_GSEA_correlation_table), Pathway_set_in_plot)
  pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[match(c(Pathway_set_in_plot, other_pathways), rownames(pancancer_Hallmark_GSEA_correlation_table)),]
  pancancer_Hallmark_GSEA_correlation_table_subset = pancancer_Hallmark_GSEA_correlation_table[,which(colnames(pancancer_Hallmark_GSEA_correlation_table) %in% CancerTYPES)]
  test_pcor = pancancer_Hallmark_GSEA_correlation_table_subset < cor_cutoff
  test_pcor = as.data.frame(test_pcor)
  test_pcor$rowsums = rowSums(test_pcor)
  pathways_to_exclude = rownames(test_pcor)[which(test_pcor$rowsums <= (minimum_number_of_cancers-1))]                                                                                                              
  cor_pathways = setdiff(rownames(pancancer_Hallmark_GSEA_correlation_table_subset), pathways_to_exclude)
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

oncoPrint(matrix_oncoprint_input,
          alter_fun = alter_fun,
          col = col,
          heatmap_legend_param = list(title= "", at = c("UP"),
                                      labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
          column_title = paste0("OncoPrint: ", paste(CancerTYPES, collapse = " "), " RNASeq expression", "\n N patients = ", number_of_patients,
                                "\n cor cutoff = ", cor_cutoff, " pathways included: ", pathway_selection, " of all cancers"),
          row_order = pathway_order,
          show_row_barplot = FALSE,
          row_barplot_width = unit(7, "in"),
          top_annotation = NULL,
          width = unit(8, "in"),
          row_title = "Hallmark pathways",
          column_title_gp = gpar(fontface = 2, fontsize = 16))

load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata")
matrix_oncoprint_input = matrix_oncoprint_input[,column_order_oncoprint]
matrix_oncoprint_input = matrix_oncoprint_input[match(pathway_order, rownames(matrix_oncoprint_input)),]
cancer_order = ICR_cluster_assignment_allcancers$Cancer[match(colnames(matrix_oncoprint_input), rownames(ICR_cluster_assignment_allcancers))]

cancer_annotation_bar = as.character(cancer_order)
df_Annotation = data.frame(cancer_annotation_bar)
colors = Cancer_color_table$color
colors = setNames(colors, Cancer_color_table$Group.1)
colors = colors[order(names(colors))]

#subset, "_cor_cutoff_", cor_cutoff, 
#"_Pathway_set_2"
#dev.new()
png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_Pathway_set_from_plot", "_", ICR_medium_excluded, "_",
           IPA_excluded, "_", subset, ".png"),res=600,height= 9,width= 18,unit="in")
oncoPrint(matrix_oncoprint_input,
          alter_fun = alter_fun,
          col = col,
          heatmap_legend_param = list(title= "", at = c("UP"),
                                      labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
          column_title = paste0("OncoPrint: ", paste(CancerTYPES, collapse = " "), " RNASeq expression", "\n N patients = ", number_of_patients,
                                ", Pathways <-0.1 for at least one third of the cancers"),
                                #"\n cor cutoff for all cancers: < ", cor_cutoff, " pathways included: ", pathway_selection, " of all cancers"),
          row_order = pathway_order,
          show_row_barplot = FALSE,
          row_barplot_width = unit(7, "in"),
          top_annotation = HeatmapAnnotation(df = df_Annotation,
                                             col = list(cancer_annotation_bar = colors),
                                             annotation_height = unit(2, "cm")),
          width = unit(8, "in"),
          row_title = "Hallmark pathways",
          column_title_gp = gpar(fontface = 2, fontsize = 16))
dev.off()

#### Extra part of the script: create a bar graph corresponding with the oncoprint
load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata")
matrix_oncoprint_input_t = t(matrix_oncoprint_input) 
oncoprint_df = as.data.frame(matrix_oncoprint_input_t)
oncoprint_df$Oncogenic_pathways = do.call(paste, c(oncoprint_df[2:ncol(oncoprint_df)], sep = ""))
oncoprint_df$Cancer = ICR_cluster_assignment_allcancers$Cancer[match(rownames(oncoprint_df), rownames(ICR_cluster_assignment_allcancers))]

ICR_High_Onco_High = oncoprint_df[oncoprint_df$`ICR cluster`=="UP" & oncoprint_df$Oncogenic_pathways != "",]
ICR_High_Onco_Low = oncoprint_df[oncoprint_df$`ICR cluster`=="UP" & oncoprint_df$Oncogenic_pathways == "",]
ICR_Low_Onco_High = oncoprint_df[oncoprint_df$`ICR cluster`=="" & oncoprint_df$Oncogenic_pathways != "",]
ICR_Low_Onco_Low = oncoprint_df[oncoprint_df$`ICR cluster`=="" & oncoprint_df$Oncogenic_pathways == "",]

png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_Pathway_set_from_plot", "_", ICR_medium_excluded, "_",
           IPA_excluded, "_", subset, "_boxplot_IHOH", ".png"),res=600,height= 5,width= 9,unit="in")
IHOH = table(ICR_High_Onco_High$Cancer)
IHOH = IHOH[which(names(IHOH) %in% CancerTYPES)]
IHOH = IHOH[order(names(IHOH))]
colors_barplot = colors[names(IHOH)]
barplot(IHOH, ylab = "Number of patients", col = colors_barplot, cex.axis = 1, cex.names = 1, cex.lab = 1, ylim = c(0, 180))
dev.off()

png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_Pathway_set_from_plot", "_", ICR_medium_excluded, "_",
           IPA_excluded, "_", subset, "_boxplot_IHOL", ".png"),res=600,height= 5,width= 9,unit="in")
IHOL = table(ICR_High_Onco_Low$Cancer)
IHOL = IHOL[which(names(IHOL) %in% CancerTYPES)]
IHOL = IHOL[order(names(IHOL))]
colors_barplot = colors[names(IHOL)]
barplot(IHOL, ylab = "Number of patients", col = colors_barplot, cex.axis = 1, cex.names = 1, cex.lab = 1, ylim = c(0, 180))
dev.off()

png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_Pathway_set_from_plot", "_", ICR_medium_excluded, "_",
           IPA_excluded, "_", subset, "_boxplot_ILOH", ".png"),res=600,height= 5,width= 9,unit="in")
ILOH = table(ICR_Low_Onco_High$Cancer)
ILOH = ILOH[which(names(ILOH) %in% CancerTYPES)]
ILOH = ILOH[order(names(ILOH))]
colors_barplot = colors[names(ILOH)]
barplot(ILOH, ylab = "Number of patients", col = colors_barplot, cex.axis = 1, cex.names = 1, cex.lab = 1, ylim = c(0, 180))
dev.off()

png(paste0("./5_Figures/Pancancer_plots/OncoPrints/", paste(CancerTYPES, collapse = "_"),
           "/Hallmark_OncoPrint_", paste(CancerTYPES, collapse = "_"), "_", pathway_selection, "_Pathway_set_from_plot", "_", ICR_medium_excluded, "_",
           IPA_excluded, "_", subset, "_boxplot_ILOL", ".png"),res=600,height= 5,width= 9,unit="in")
ILOL = table(ICR_Low_Onco_Low$Cancer)
ILOL = ILOL[which(names(ILOL) %in% CancerTYPES)]
ILOL = ILOL[order(names(ILOL))]
colors_barplot = colors[names(ILOL)]
barplot(ILOL, ylab = "Number of patients", col = colors_barplot, cex.axis = 1, cex.names = 1, cex.lab = 1, ylim = c(0, 180))
dev.off()

### Relevant percentages graph

ICR_High = IHOH + IHOL
ICR_Low = ILOH + ILOL

Percentages_Onco_High_in_ICR_High = IHOH/ICR_High * 100
Percentages_Onco_High_in_ICR_Low = ILOH/ICR_Low * 100

Plot_df = data.frame(Percentages_Onco_High_in_ICR_High, Percentages_Onco_High_in_ICR_Low)
Plot_df = Plot_df[,c(1,2,4)]
colnames(Plot_df) = c("Cancer", "Percentage Onco High in ICR High", "Percentage Onco High in ICR Low")
Plot_df_long = melt(Plot_df)

dev.new()
ggplot(Plot_df_long, aes(Cancer, value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("Percentage") +
  xlab("Cancer")
