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

rm(list = ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.bioconductor.packages = c("ComplexHeatmap")
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/oncoPrint.R"))                                                                        # source adapted version of oncoPrint version (Other functions of ComplexHeatmap packages are still required still required)

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
subset = "COR_COEF"                                                                                                      # Options: "ALL_SIG" or "INV_COR_SIG" or "POS_COR_SIG"
cor_cutoff = 0

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

i=1
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
    #cor_pathways_pos = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] > cor_cutoff),] #comment-out *1, when no # comment-out *2
    cor_pathways_neg = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < -cor_cutoff),] 
    if(is.matrix(cor_pathways_neg)){
      #cor_pathways = rbind(cor_pathways_pos, cor_pathways_neg) #comment-out *1
      cor_pathways = cor_pathways_neg #comment-out *2
    }else{next
      #cor_pathways = cor_pathways_pos
    }
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
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score > z_score_upregulation] = 100
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score < z_score_downregulation] = -100
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score >= z_score_downregulation & Hallmark.enrichment.z.score <= z_score_upregulation] = 0
    
    # Convert individual values to new character value
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score == 100] = "UP"
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score == -100] = NA
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score == 0] = NA
    
    # Filter to only plot the selected cor_pathways
    Hallmark.enrichment.z.score = Hallmark.enrichment.z.score[which(rownames(Hallmark.enrichment.z.score) %in% cor_pathways),]
    
    Hallmark.enrichment.z.score = rbind(Hallmark.enrichment.z.score, table_cluster_assignment$HML_cluster[match(colnames(Hallmark.enrichment.z.score), row.names(table_cluster_assignment))])
    rownames(Hallmark.enrichment.z.score)[nrow(Hallmark.enrichment.z.score)] = "ICR cluster"
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score == "ICR Low" | Hallmark.enrichment.z.score == "ICR Medium"] = NA
    Hallmark.enrichment.z.score[Hallmark.enrichment.z.score == "ICR High"] = "UP"
    
    Hallmark.enrichment.z.score[is.na(Hallmark.enrichment.z.score)] = ""
    
    matrix = Hallmark.enrichment.z.score
  }
  
  if(expression_units == "quantiles"){
    quantfun = function(x) as.integer(cut(x, breaks = quantile(x, probs = 0:4/4), include.lowest=TRUE))
    Hallmark.expression.quantiles =  apply(Hallmark.enrichment.score, MARGIN = 1, FUN = quantfun)
    Hallmark.expression.quantiles = t(Hallmark.expression.quantiles)
    colnames(Hallmark.expression.quantiles) = colnames(Hallmark.enrichment.score)
    
    Hallmark.expression.quantiles[Hallmark.expression.quantiles != 4] = NA
    Hallmark.expression.quantiles[Hallmark.expression.quantiles == 4] = "UP"
    
    # Filter to only plot the selected cor_pathways
    Hallmark.expression.quantiles = Hallmark.expression.quantiles[which(rownames(Hallmark.expression.quantiles) %in% cor_pathways),]
   
    Hallmark.expression.quantiles = rbind(Hallmark.expression.quantiles, table_cluster_assignment$HML_cluster[match(colnames(Hallmark.expression.quantiles), row.names(table_cluster_assignment))])
    rownames(Hallmark.expression.quantiles)[nrow(Hallmark.expression.quantiles)] = "ICR cluster"
    Hallmark.expression.quantiles[Hallmark.expression.quantiles == "ICR Low" | Hallmark.expression.quantiles == "ICR Medium"] = NA
    Hallmark.expression.quantiles[Hallmark.expression.quantiles == "ICR High"] = "UP"
    
    Hallmark.expression.quantiles[is.na(Hallmark.expression.quantiles)] = ""
    
    matrix = Hallmark.expression.quantiles
  }
  
  
  dir.create(paste0("./5_Figures/OncoPrints/"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/", download.method), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/", download.method), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/", download.method, "/", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff), showWarnings = FALSE)
  
  png(paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/", download.method, "/", expression_units,"_", z_score_upregulation, "_", subset, "_", cor_cutoff,
             "/Hallmark_OncoPrint_", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff, "_", Cancer, ".png"),res=600,height= 7,width= 18,unit="in")
  alter_fun = list(
    background = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = "grey", col = NA)),
    UP = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.9, gp = gpar(fill = col["UP"], col = NA)))
  
  col = c(UP = "red", DOWN = "blue")
    
  plot = oncoPrint(matrix,
                   alter_fun = alter_fun,
                   col = col,
                   heatmap_legend_param = list(title= "", at = c("UP"),
                                        labels = c("Upregulation"), nrow = 1, title_position = "leftcenter"),
                   column_title = paste0("OncoPrint: ", Cancer, " RNASeq expression", "\n N patients = ", nrow(table_cluster_assignment)),
                   row_order = pathway_order,
                   show_row_barplot = FALSE,
                   row_barplot_width = unit(7, "in"),
                   top_annotation = NULL,
                   width = unit(8, "in"),
                   row_title = "Hallmark pathways",
                   column_title_gp = gpar(fontface = 2, fontsize = 16))

  draw(plot, heatmap_legend_side = "bottom", row_sub_title_side = "left")
  dev.off()
  
  oncoPrintMatrix = matrix[match(pathway_order, rownames(matrix)),column_order_oncoprint]
  oncoPrintMatrix = rbind(oncoPrintMatrix, table_cluster_assignment$HML_cluster[match(colnames(oncoPrintMatrix), row.names(table_cluster_assignment))])
  assign(paste0(Cancer, "_OncoPrintMatrix"), oncoPrintMatrix)
}

rm(list = "oncoPrintMatrix")
AlloncoPrintMatrixes = grep("OncoPrintMatrix", ls(), value = TRUE)
save(list = AlloncoPrintMatrixes , file = paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/", download.method, "/", expression_units,"_", z_score_upregulation, "_", subset, "_", cor_cutoff,
                                                 "/Hallmark_OncoPrint_", expression_units, "_", z_score_upregulation, "_", subset, "_", cor_cutoff, ".Rdata"))

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"

TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)

CancerTYPES = "ALL"

if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

for(i in 1:N.sets){
  Cancer = CancerTYPES[i]
  AlloncoPrintMatrixes = grep("OncoPrintMatrix", ls(), value = TRUE)
  OncoPrintMatrix = get(AlloncoPrintMatrixes[i])
  write.csv(OncoPrintMatrix, file = paste0("./5_Figures/OncoPrints/Hallmark_OncoPrints/TCGA_Assembler/z_score_1.5_COR_COEF_0/", Cancer,"_OncoPrintMatrix.csv"))
}
