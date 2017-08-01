#################################################################
###
### Make heatmaps of hallmark pathways that have a significant
### inverse relation with ICR gene expression.
### Input files:
### "./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_only_significant.Rdata"
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata"
### "./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
### Cancer, "_ICR_cluster_assignment_k2-6.Rdata"
### Output files:
### "./5_Figures/Heatmaps/ICR_low_pathway_Heatmaps/", download.method, "/", Cancer, "/ICR_low_Pathway_Heatmap.3_RNASeq_",Cancer,".png"
###
#################################################################

rm(list = ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path,"R tools/heatmap.3.R"))

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"
my.palette = colorRampPalette(c("blue", "white", "red"))(n = 297)
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297) 
ColsideLabels = c("HML ICR clusters", "Bindea clusters")
Legend = c("ICR Low","ICR Med","ICR High", "Bindea Low", "Bindea High")
Legend_colors = c("blue","green","red", "pink", "purple")
Log_file = paste0("./1_Log_Files/3.11_Heatmap_ICR_low_pathways/3.11_Heatmap_ICR_low_pathways_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)
load("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark_only_significant.Rdata")

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
  
  sign_inv_pathways = pancancer_Hallmark_GSEA_correlation_table[which(pancancer_Hallmark_GSEA_correlation_table[,Cancer] < 0), ]
  sign_inv_pathways = rownames(sign_inv_pathways)
  sign_inv_pathways = c("ICR_genes", sign_inv_pathways)
  
  Hallmark.enrichment.z.score.ICR.low.pathways = Hallmark.enrichment.z.score[which(rownames(Hallmark.enrichment.z.score) %in% sign_inv_pathways),]
  
  annotation.blot = as.matrix(annotation[,c("HML_cluster.col","Bindea_cluster.col"), drop= FALSE])
  annotation.blot = annotation.blot[row.names(annotation),]
  
  dir.create(paste0("./5_Figures/Heatmaps/ICR_low_pathway_Heatmaps"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/Heatmaps/ICR_low_pathway_Heatmaps/", download.method), showWarnings = FALSE)
  
  png(paste0("./5_Figures/Heatmaps/ICR_low_pathway_Heatmaps/", download.method, "/ICR_low_Pathway_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height =7,width=10,unit="in")
  heatmap.3(as.matrix(Hallmark.enrichment.z.score.ICR.low.pathways),
            main = paste0(Cancer, "\n Hallmark ssGSEA Heatmap"),
            cex.main = 2,
            col= my.palette,
            ColSideLabs = ColsideLabels,
            font_size_col_Labs = 1.5,
            ColSideColors= annotation.blot,
            #breaks = my.colors,
            scale = "row",
            Colv=NA,
            labCol=NA,
            side.height.fraction = 0.5,
            cexRow= 1,
            margins=c(13,25))
  par(lend = 1)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, " v2.0.3.\n",
                          "Hallmark enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.50)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.2)
  dev.off()
  
  # Make correlation plot of pathways that are inversely associated with ICR
  
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/Correlation_matrixes_Bindea_xCell_Hallmark_pearson_", Cancer, ".Rdata"))
  
  Hallmark_GSEA_inv_ICR_cor = Hallmark_GSEA_cor[which(colnames(Hallmark_GSEA_cor) %in% sign_inv_pathways), which(rownames(Hallmark_GSEA_cor) %in% sign_inv_pathways)]
  indexes = which(colnames(Hallmark_GSEA_cor) %in% sign_inv_pathways)
  
  Hallmark_GSEA_inv_ICR_cor_sign = Hallmark_GSEA_cor_sign[[1]][which(colnames(Hallmark_GSEA_cor) %in% sign_inv_pathways), which(rownames(Hallmark_GSEA_cor) %in% sign_inv_pathways)]
  
  dir.create(paste0("./5_Figures/Correlation_plots/Sign_inv_pathways_ICR_Correlation_plots/"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/Correlation_plots/Sign_inv_pathways_ICR_Correlation_plots/", download.method), showWarnings = FALSE)

  png(paste0("./5_Figures/Correlation_plots/Sign_inv_pathways_ICR_Correlation_plots/", download.method, 
             "/Sign_inv_pathways_low_ICR_Correlation_plot_",Cancer,".png"), res=600,height=6.5,width=6,unit="in")
  cex.before <- par("cex")
  par(cex = 0.45)
  lims=c(-1,1)
  if (length(Hallmark_GSEA_inv_ICR_cor[Hallmark_GSEA_inv_ICR_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(Hallmark_GSEA_inv_ICR_cor),color = c(rep("black",length(sign_inv_pathways))),stringsAsFactors = FALSE)
  annotation$color[annotation$gene %in% c("ICR_genes")] = "#CC0506"
  annotation = annotation[corrMatOrder(Hallmark_GSEA_inv_ICR_cor,order="FPC"),]
  
  mean_correlation = round(mean(Hallmark_GSEA_inv_ICR_cor),2)
  corrplot.mixed (Hallmark_GSEA_inv_ICR_cor,
                  #type="lower",
                  #p.mat = Hallmark_GSEA_inv_ICR_cor_sign,                                                                            # add significance to correlations
                  col = colpattern,
                  lower = "square",
                  upper ="number",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.5,
                  tl.cex = 1.2,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(6,4.1,7,5))
  title(main = list(paste0(Cancer, "\n Correlation between pathways upregulated in ICR low\n Number of patients = ", nrow(annotation.blot)), cex = 1.5), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.65)
  par(cex = cex.before)
  dev.off()
  
  dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation"), showWarnings = FALSE)
  save(Hallmark_GSEA_inv_ICR_cor, Hallmark_GSEA_inv_ICR_cor_sign, file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/",
                                                        "Correlation_matrix_Sign_Inv_Correlating_Pathways_ICR_", test, "_", Cancer, ".Rdata"))
}