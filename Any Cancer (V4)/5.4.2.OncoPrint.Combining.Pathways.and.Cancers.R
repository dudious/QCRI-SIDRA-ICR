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
#rm(list=ls())
rm(list=setdiff(ls(), "Cancers_for_oncoprint"))

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
Group_1 = c("[HM] KRAS SIGNALING DN", "[HM] HEDGEHOG SIGNALING", "[HM] NOTCH SIGNALING",
            "[HM] WNT BETA CATENIN SIGNALING")

Group_2 = c("[TBI] BARRIER GENES", "[HM] GLYCOLYSIS", "[HM] ESTROGEN RESPONSE",
            "[TPW] Mesenchymal Transition", "[HM] TGF BETA SIGNALING",
            "[HM] UV RESPONSE DN")

Group_3 = c("[TBI] MAPK UP GENES", "[HM] MITOTIC SPINDLE", "[TBI] PHOSPHOLIPASE",
            "[TPW] NOS1 Signature", "[HM] E2F TARGETS", "[HM] G2M CHECKPOINT",
            "[HM] OXIDATIVE PHOSPHORYLATION", "[HM] MYC TARGETS", "[HM] DNA REPAIR")

#Group_yellow = c("GBM", "LGG","OV", "SARC", "THYM","UCS")
#Group_red = c("LUAD", "LUSC","MESO","STAD", "UCEC")
#Group_blue = c("CESC", "ESCA", "HNSC","SKCM")
#Group_green = c("BLCA", "COAD", "KICH", "KIRC", "KIRP",
                #"LIHC", "PAAD", "PCPG", "READ", "THCA",
                #"UVM")

# Set Parameters
CancerTYPES = Cancers_for_oncoprint                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("DLBC", "PRAD", "CHOL")
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                                                                                                       # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
#pathway_selection = "intersection"                                                                                      # Specify whether you want to combine pathways that overlap between all cancers that are used as input (intersection)
                                                                                                                        # or if you want include all pathways that appear in either one of the cancers used as input (union). 
download.method = "TCGA_Assembler"
Log_file = paste0("./1_Log_Files/5.4_OncoPrint_Combining_Cancers_RNASeq/3.12_OncoPrint_RNASeq_Log_File_",               # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
#subset = "<0.1_for_at_least_one_third_of_cancers"                                                                                   # Options: "COR_COEF" or c("[HM] WNT BETA CATENIN SIGNALING", "[HM] NOTCH SIGNALING"..etc)
#minimum_number_of_cancers = length(CancerTYPES)/3
cor_cutoff = -0.1                                                                                                        # Cor-cutoff to define which pathways are selected for each cancer
ICR_medium_excluded = "ICR_medium_excluded"                                                                              # Options: "ICR_medium_excluded" or "all_included"
IPA_excluded = "IPA_excluded"

# Load data
load(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v5/TCGA_Assembler/z_score_1.5_COR_COEF_", cor_cutoff,
            "/All_Hallmark_pathways_Oncoprint_matrixesz_score_1.5_COR_COEF_", cor_cutoff, ".Rdata"))

load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_irrespective_of_significance_IPA_excluded.Rdata")

# Combine datasets of interest
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))){next}
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

# Filter to only plot the selected cor_pathways
cor_pathways = c(Group_3, "ICR cluster")
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
                                "\n cor cutoff = ", cor_cutoff, " pathways included: ", paste(Group_3, collapse = " ") , " of selected cancers"),
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
#png(paste0("./5_Figures/Pancancer_plots/OncoPrints_selecting_pathways_first/All_cancers/All_cancers_", "Combined_Oncoprint_", 
           #paste(Group_3, collapse = "_"), ".png"),res=600,height= 9,width= 18,unit="in")
png(paste0("./5_Figures/Pancancer_plots/OncoPrints_selecting_pathways_first/", paste(Group_3, collapse = "_"), "Combined_Oncoprint", ".png"),res=600,height= 9,width= 18,unit="in")
oncoPrint(matrix_oncoprint_input,
          alter_fun = alter_fun,
          col = col,
          heatmap_legend_param = list(title= "", at = c("UP"),
                                      labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
          column_title = paste0("OncoPrint: ", "All cancers RNASeq expression", "\n N patients = ", number_of_patients,
                                "\n cor cutoff = ", cor_cutoff, " pathways included: ", paste(Group_3, collapse = " ") , " of selected cancers"),
                                #"\n cor cutoff for all cancers: < ", cor_cutoff, " pathways included: ", pathway_selection, " of selected cancers"),
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

included_cancers = setdiff(CancerTYPES, Cancer_skip)

IHOH = table(ICR_High_Onco_High$Cancer)
IHOH = IHOH[which(names(IHOH) %in% included_cancers)]
IHOH = IHOH[order(names(IHOH))]

IHOL = table(ICR_High_Onco_Low$Cancer)
IHOL = IHOL[which(names(IHOL) %in% included_cancers)]
IHOL = IHOL[order(names(IHOL))]

ILOH = table(ICR_Low_Onco_High$Cancer)
ILOH = ILOH[which(names(ILOH) %in% included_cancers)]
ILOH = ILOH[order(names(ILOH))]

ILOL = table(ICR_Low_Onco_Low$Cancer)
ILOL = ILOL[which(names(ILOL) %in% included_cancers)]
ILOL = ILOL[order(names(ILOL))]

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

Plot_df$Select = Plot_df$`Percentage Onco High in ICR High` < Plot_df$`Percentage Onco High in ICR Low`

Cancers_for_oncoprint = as.character(Plot_df$Cancer[Plot_df$Select == TRUE])
### Use the selected Cancers_for_oncoprint to run this script again CancerTYPES = Cancers_for_oncoprint

