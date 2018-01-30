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


## In pipeline: Run script 5.1 first!

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))

required.bioconductor.packages = c("ComplexHeatmap", "circlize", "gridExtra")
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/oncoPrint.R"))                                                                        # source adapted version of oncoPrint version (Other functions of ComplexHeatmap packages are still required still required)
source(paste0(code_path, "R tools/heatmap.3.R"))

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"
#ColsideLabels = c("HML ICR clusters", "Bindea clusters")
#Legend = c("ICR Low","ICR Med","ICR High", "Bindea Low", "Bindea High")
#Legend_colors = c("blue","green","red", "pink", "purple")
Log_file = paste0("./1_Log_Files/3.12_OncoPrint_RNASeq/3.12_OncoPrint_RNASeq_Log_File_",                                # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"

expression_units = "z_score"                                                                                             # ("z_score" or "quantiles". Define the matrix to use for generation of OncoPrint (z-score matrix or quantile expression values)
z_score_upregulation = 1.5
z_score_downregulation = -1.5
subset = "COR_COEF"                                                                                                      # Options: "ALL_SIG" or "INV_COR_SIG" or "POS_COR_SIG" or "COR_COEF"
cor_cutoff = 0.1
ICR_medium_excluded = "ICR_medium_excluded"                                                                              # Options: "ICR_medium_excluded" or "all_included"
IPA_excluded = "IPA_excluded"                                                                                            # If all IPA pathways need to be excluded, set IPA_excluded

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)
if(subset == "ALL_SIG" | subset == "INV_COR_SIG" | subset == "POS_COR_SIG"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_only_significant.Rdata")
}
if(subset == "COR_COEF"){
  load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_irrespective_of_significance.Rdata")
}

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

i=2
for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))){next}
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  load(paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
              "_Bindea_xCell_HallmarkPathways.Rdata"))
  
  # Selection of pathways to display in plot
  if(subset == "ALL_SIG"){
    cor_pathways = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < 0 | pancancer_Hallmark_GSEA_correlation_table[,Cancer] > 0), ]
  }
  if(subset == "INV_COR_SIG"){
    cor_pathways = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < 0), ]
  }
  if(subset == "POS_COR_SIG"){
    cor_pathways = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] > 0), ]
  }
  if(subset == "ALL"){
    cor_pathways = pancancer_Hallmark_GSEA_correlation_table
  }
  if(subset == "COR_COEF"){
    cor_pathways_neg = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < cor_cutoff),] 
    if(length(cor_pathways_neg) == 0){
      next} 
    if(length(cor_pathways_neg) > 0){
      cor_pathways = rownames(pancancer_Hallmark_GSEA_correlation_table)[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < cor_cutoff)] 
    }
  }
  
  if(IPA_excluded == "IPA_excluded"){
    cor_pathways = cor_pathways[grep(pattern = "IPA", cor_pathways, invert = TRUE)]
  }
  if(length(cor_pathways) == 0){
    next
  }
  
  if("ICR_genes" %in% cor_pathways & "ICR_score" %in% cor_pathways){
    pathway_order = c("ICR cluster", "ICR_genes", cor_pathways[-which(cor_pathways == "ICR_score" | cor_pathways == "ICR_genes")])
  }else{
    pathway_order = c("ICR cluster", cor_pathways)
  }
  
  if(expression_units == "z_score"){
    ## Translate z_score to upregulation or downregulation
    # First set all to new numeric value (100, -100 or 0)
    oncoprint_matrix = Hallmark.enrichment.z.score
    oncoprint_matrix[oncoprint_matrix > z_score_upregulation] = 100
    oncoprint_matrix[oncoprint_matrix < z_score_downregulation] = -100
    oncoprint_matrix[oncoprint_matrix >= z_score_downregulation & oncoprint_matrix <= z_score_upregulation] = 0
    
    # Convert individual values to new character value
    oncoprint_matrix[oncoprint_matrix == 100] = Cancer
    oncoprint_matrix[oncoprint_matrix == -100] = NA
    oncoprint_matrix[oncoprint_matrix == 0] = NA
    
    # Save oncoprint_matrix with all pathways to the working environment
    assign(paste0(Cancer, "_all_pathways_onco_matrix"), oncoprint_matrix)
    
    # Filter to only plot the selected cor_pathways
    oncoprint_matrix = oncoprint_matrix[which(rownames(oncoprint_matrix) %in% cor_pathways), , drop = FALSE]
    oncoprint_matrix = rbind(oncoprint_matrix, table_cluster_assignment$HML_cluster[match(colnames(oncoprint_matrix), row.names(table_cluster_assignment))])
    rownames(oncoprint_matrix)[nrow(oncoprint_matrix)] = "ICR cluster"
    oncoprint_matrix[oncoprint_matrix == "ICR Low" | oncoprint_matrix == "ICR Medium"] = NA
    oncoprint_matrix[oncoprint_matrix == "ICR High"] = "UP"
    
    oncoprint_matrix[is.na(oncoprint_matrix)] = ""
    
    ## save OncoPrint
    assign(paste0(Cancer, "_oncoprint_matrix"), oncoprint_matrix)
  }
}

dir.create(paste0("./5_Figures/OncoPrints/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/", download.method, "/", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff), showWarnings = FALSE)

StringsToSelect = "_oncoprint_matrix|_matrix_heatmap_oncoprint"
AlloncoPrintMatrixes = grep(pattern = StringsToSelect, ls(), value = TRUE)

save(list = c(AlloncoPrintMatrixes), file = paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/", download.method, "/", expression_units,"_", z_score_upregulation, "_", subset, "_", cor_cutoff,
                                                   "/Hallmark_OncoPrint_", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff, "_", ICR_medium_excluded, ".Rdata"))

StringsToSelect2 = "_all_pathways_onco_matrix"
AllPathwaysAlloncoPrintMatrixes = grep(pattern = StringsToSelect2, ls(), value = TRUE)
save(list = c(AllPathwaysAlloncoPrintMatrixes), file = paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_for_multicolor/", download.method, "/", expression_units,"_", z_score_upregulation, "_", subset, "_", cor_cutoff,
                                                              "/All_Hallmark_pathways_Oncoprint_matrixes", expression_units, "_", z_score_upregulation, "_", subset, "_",
                                                              cor_cutoff, ".Rdata"))
