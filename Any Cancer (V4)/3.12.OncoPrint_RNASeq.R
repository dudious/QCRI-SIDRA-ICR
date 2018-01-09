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

#setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

#code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

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
cor_cutoff = 0
ICR_medium_excluded = "ICR_medium_excluded"                                                                              # Options: "ICR_medium_excluded" or "all_included"

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

i=3
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
    cor_pathways_neg = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < -cor_cutoff),] 
    if(is.matrix(cor_pathways_neg)){
      cor_pathways = cor_pathways_neg
    }else{next}
  }
  cor_pathways = row.names(cor_pathways)
  
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
    oncoprint_matrix[oncoprint_matrix == 100] = "UP"
    oncoprint_matrix[oncoprint_matrix == -100] = NA
    oncoprint_matrix[oncoprint_matrix == 0] = NA
    
    # Filter to only plot the selected cor_pathways
    oncoprint_matrix = oncoprint_matrix[which(rownames(oncoprint_matrix) %in% cor_pathways),]
    
    oncoprint_matrix = rbind(oncoprint_matrix, table_cluster_assignment$HML_cluster[match(colnames(oncoprint_matrix), row.names(table_cluster_assignment))])
    rownames(oncoprint_matrix)[nrow(oncoprint_matrix)] = "ICR cluster"
    oncoprint_matrix[oncoprint_matrix == "ICR Low" | oncoprint_matrix == "ICR Medium"] = NA
    oncoprint_matrix[oncoprint_matrix == "ICR High"] = "UP"
    
    oncoprint_matrix[is.na(oncoprint_matrix)] = ""
    
    complete_matrix = oncoprint_matrix
  }
  
  #if(expression_units == "quantiles"){
  #quantfun = function(x) as.integer(cut(x, breaks = quantile(x, probs = 0:4/4), include.lowest=TRUE))
  #Hallmark.expression.quantiles =  apply(Hallmark.enrichment.score, MARGIN = 1, FUN = quantfun)
  #Hallmark.expression.quantiles = t(Hallmark.expression.quantiles)
  #colnames(Hallmark.expression.quantiles) = colnames(Hallmark.enrichment.score)
  
  #Hallmark.expression.quantiles[Hallmark.expression.quantiles != 4] = NA
  #Hallmark.expression.quantiles[Hallmark.expression.quantiles == 4] = "UP"
  
  # Filter to only plot the selected cor_pathways
  #Hallmark.expression.quantiles = Hallmark.expression.quantiles[which(rownames(Hallmark.expression.quantiles) %in% cor_pathways),]
  
  #Hallmark.expression.quantiles = rbind(Hallmark.expression.quantiles, table_cluster_assignment$HML_cluster[match(colnames(Hallmark.expression.quantiles), row.names(table_cluster_assignment))])
  #rownames(Hallmark.expression.quantiles)[nrow(Hallmark.expression.quantiles)] = "ICR cluster"
  #Hallmark.expression.quantiles[Hallmark.expression.quantiles == "ICR Low" | Hallmark.expression.quantiles == "ICR Medium"] = NA
  #Hallmark.expression.quantiles[Hallmark.expression.quantiles == "ICR High"] = "UP"
  
  #Hallmark.expression.quantiles[is.na(Hallmark.expression.quantiles)] = ""
  
  #matrix = Hallmark.expression.quantiles
  #}
  
  dir.create(paste0("./5_Figures/OncoPrints/"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/", download.method), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/", download.method), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/", download.method, "/", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff), showWarnings = FALSE)
  
  alter_fun = list(
    background = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = "grey", col = NA)),
    UP = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UP"], col = NA)))
  
  col = c(UP = "red", DOWN = "blue")
  
  if(ICR_medium_excluded == "all_included"){
    number_of_patients = nrow(table_cluster_assignment)
  }
  if(ICR_medium_excluded == "ICR_medium_excluded"){
    number_of_patients = nrow(table_cluster_assignment[table_cluster_assignment$HML_cluster != "ICR Medium",])
  }
  
  # Run oncoPrint function first (without saving in png) to get the column_order_oncoprint (=order of patients/tumors), which is needed for the heatmap next to the oncoprint
  # Very important to use the complete cohort (without filtering out ICR Medium Tumors), because column_order_oncoprint is an integer vector generated in oncoPrint function.
  # First use this to order the heatmap matrix, to subsequently filter out medium tumors (if needed).
  First_oncoprint = oncoPrint(complete_matrix,
                              alter_fun = alter_fun,
                              col = col,
                              heatmap_legend_param = list(title= "", at = c("UP"),
                                                          labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
                              column_title = paste0("OncoPrint: ", Cancer, " RNASeq expression", "\n N patients = ", number_of_patients),
                              row_order = pathway_order,
                              show_row_barplot = FALSE,
                              row_barplot_width = unit(7, "in"),
                              top_annotation = NULL,
                              width = unit(8, "in"),
                              row_title = "Hallmark pathways",
                              column_title_gp = gpar(fontface = 2, fontsize = 16))
  
  # Generate matrix for heatmap from Hallmark enrichment z scores
  pathway_order_heatmap = pathway_order
  pathway_order_heatmap[1] = "ICR_genes"
  matrix_heatmap_oncoprint = Hallmark.enrichment.z.score[match(pathway_order_heatmap, rownames(Hallmark.enrichment.z.score)), column_order_oncoprint]
  #rownames(matrix_heatmap_oncoprint)[1] = "ICR cluster"
  
  # Filter out Medium Tumors if needed (specify in parameters in top of this script)
  if(ICR_medium_excluded == "ICR_medium_excluded"){
    excluded_patients = rownames(table_cluster_assignment)[table_cluster_assignment$HML_cluster == "ICR Medium"]
    matrix_heatmap_oncoprint = matrix_heatmap_oncoprint[,-which(colnames(matrix_heatmap_oncoprint) %in% excluded_patients)]
  }
  
  complete_matrix = complete_matrix[rownames(matrix_heatmap_oncoprint),column_order_oncoprint]
  
  if(ICR_medium_excluded == "all_included"){
    oncoprint_matrix = complete_matrix
  }
  
  if(ICR_medium_excluded == "ICR_medium_excluded"){
    oncoprint_matrix = complete_matrix[, -which(colnames(complete_matrix) %in% excluded_patients)]
  }
  
  # Run oncoPrint function for a second time to generate the figure, no algorithm used in this function: 
  png(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/", download.method, "/", expression_units,"_", z_score_upregulation, "_", subset, "_", cor_cutoff,
             "/Hallmark_OncoPrint_", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff, "_", ICR_medium_excluded, "_", Cancer, ".png"),res=600,height= 9,width= 18,unit="in")
  
  #col_fun_oncoprint = colorRamp2(c("UP", ""), c("red", "grey"))
  col_fun_heatmap = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  plot1 = oncoPrint(oncoprint_matrix,
                    alter_fun = alter_fun,
                    col = col,
                    column_order = NULL,
                    #heatmap_legend_param = list(title= "", at = c("UP"),
                    #labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
                    show_row_barplot = FALSE,
                    show_heatmap_legend = FALSE,
                    top_annotation = NULL,
                    width = unit(8, "in"),
                    show_pct = FALSE,
                    row_order = 1:nrow(oncoprint_matrix),
                    #row_title = "Hallmark pathways",
                    column_title_gp = gpar(fontface = 2, fontsize = 16))
  
  plot2 = Heatmap(matrix_heatmap_oncoprint, name = "expression", 
                  cluster_rows = FALSE, cluster_columns = FALSE, 
                  show_row_dend = FALSE,  show_column_dend = FALSE,
                  column_order = 1:ncol(matrix_heatmap_oncoprint), row_order = 1:nrow(matrix_heatmap_oncoprint),
                  col = col_fun_heatmap,
                  #row_title = "Hallmark pathways",
                  show_heatmap_legend = FALSE,
                  show_column_names = FALSE,
                  width = unit(8, "in"),
                  column_title_gp = gpar(fontface = 2, fontsize = 16))
  
  #lgd_oncoprint = Legend(c("UP"), col = "red", labels = c("Upregulation", "No Upregulation"), nrow = 1, title_position = "leftcenter")
  lgd_heatmap = Legend(at = c(-2, -1, 0, 1, 2), col_fun = col_fun_heatmap, title = "heatmap legend", title_position = "topleft")
  
  pushViewport(viewport(layout = grid.layout(ncol = 2, nrow = 3, widths = unit(c(8, 7.5), c("in", "in")),
                                             heights = unit(c(0.5,3.5, 3.5), c("in", "in", "in")))))
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(plot1, newpage = FALSE)
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  draw(plot2, newpage = FALSE)
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid.draw(lgd_heatmap)
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(paste0(Cancer, " RNASeq expression (OncoPrint)", "\n N patients = ", number_of_patients, "\n",
                   ICR_medium_excluded), x = unit(3, "in"), gp = gpar(fontface = "bold"))
  upViewport()
  
  dev.off()
  
  ## save OncoPrint
  assign(paste0(Cancer, "_oncoprint_matrix"), oncoprint_matrix)
  assign(paste0(Cancer, "_matrix_heatmap_oncoprint"), matrix_heatmap_oncoprint)
}

StringsToSelect = "_oncoprint_matrix|_matrix_heatmap_oncoprint"
AlloncoPrintMatrixes = grep(pattern = StringsToSelect, ls(), value = TRUE)

save(list = c(AlloncoPrintMatrixes), file = paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints_v3/", download.method, "/", expression_units,"_", z_score_upregulation, "_", subset, "_", cor_cutoff,
                                                 "/Hallmark_OncoPrint_", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff, "_", ICR_medium_excluded, ".Rdata"))
