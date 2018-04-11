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
Pathways = "ALL"
CancerTYPES = "ALL"
Pathway_skip = ""                                                                                                        
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,
                  "/5.3.Pancancer.Clustering.Oncogenic.Pathways/5.3.Pancancer.Clustering.Oncogenic.Pathways_",          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("ICR clusters", "Cancers", "Proliferation cluster")
Legend = c("ICR Low","ICR Med","ICR High")

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark.Pancancer.Rdata"))
load(paste0("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

if (Pathways == "ALL") { 
  Pathways = rownames(Hallmark_enrichment_z_score_all)
}
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathway_Clustering"), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/5.3.Pancancer.Clustering.Oncogenic.Pathways"), showWarnings = FALSE)

cat("This is a log file for clustering Pancancer samples by expression of oncogenic pathways",
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("Pathways = ", Pathways),
    paste0("Pathway_skip = ", Pathway_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    file = Log_file,
    append = FALSE, sep= "\n")

N.pathways = length(Pathways)

Hallmark_and_ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers

for (i in 1:N.pathways){
  Pathway = Pathways[i]
  Enrichment_score_all.z.score = Hallmark_enrichment_z_score_all[which(rownames(Hallmark_enrichment_z_score_all) == Pathway), , drop = FALSE]
  sHc = hclust(ddist <- dist(t(Enrichment_score_all.z.score)), method = "ward.D2")
  plot(sHc,labels=FALSE)
  # Pancancer Oncogenic pathway classification
  annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer") ,drop=FALSE]
  annotation$pathway_cluster = cutree(sHc,k = 2)[match(rownames(annotation),names(cutree(sHc,k = 2)))]
  annotation = annotation[colnames(Enrichment_score_all.z.score),]
  annotation$pathway_score = Enrichment_score_all.z.score[1,]
  pathway_cluster_means = aggregate(pathway_score~pathway_cluster,data=annotation,FUN=mean)
  pathway_cluster_means = pathway_cluster_means[order(pathway_cluster_means$pathway_score),]
  pathway_cluster_means$pathway_cluster_name = c(paste0(Pathway, " Low"),paste0(Pathway, " High"))
  annotation$pathway_cluster_name = pathway_cluster_means$pathway_cluster_name[match(annotation$pathway_cluster,
                                                                                     pathway_cluster_means$pathway_cluster)]
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"
  
  annotation$Cancer.col = Cancer_color_table$color[match(annotation$Cancer, Cancer_color_table$Group.1)]
  
  #Cancer_order = CancerTYPES[-which(CancerTYPES =="LAML")]
  #ICR_order = c("ICR Low","ICR Medium","ICR High")
  #Pathway_order = c(paste0(Pathway," Low"), paste0(Pathway, " High"))
  #annotation = annotation[order(match(annotation$Cancer, Cancer_order)),]
  #annotation = annotation[order(match(annotation$HML_cluster,ICR_order)),]
  
  annotation$Pathway.col[annotation$pathway_cluster_name == paste0(Pathway, " High")] = "#B97474"
  annotation$Pathway.col[annotation$pathway_cluster_name == paste0(Pathway, " Low")] = "#ADDFAD"
  annotation.blot = as.matrix(annotation[,c("HML_cluster.col","Cancer.col", "Pathway.col"), drop = FALSE])
  Enrichment_score_all.z.score = rbind(Enrichment_score_all.z.score, rep(0, 9475))
  
  ## Plot prep
  Legend = CancerTYPES[-which(CancerTYPES == "LAML")]
  rownames(Cancer_color_table) = Cancer_color_table$Group.1
  Cancer_color_table = Cancer_color_table[CancerTYPES[-which(CancerTYPES == "LAML")],]
  Legend_colors = c(Cancer_color_table$color)
  
  Legend2 = c("ICR Low","ICR Med","ICR High")
  Legend_colors2 = c("blue","green","red")
  
  Legend3 = c(paste0(Pathway, " High"), paste0(Pathway, " Low"))
  Legend_colors3 = c("#B97474", "#ADDFAD")
  
  ### Plotting
  
  dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathway_Clustering"), showWarnings = FALSE)
  png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathway_Clustering/", gsub("/", "_", Pathway),
             "_Clusters_determined_PanCancer", "_Heatmap_RNASeq_Pancancer_HML.png"), res = 600, height = 8, width = 15, unit = "in")
  heatmap.3((as.matrix(Enrichment_score_all.z.score)),
            main= paste0("Pancancer enrichment scores \nssGSEA/", Pathway),
            col=my.palette,
            ColSideColors=annotation.blot,
            font_size_col_Labs = 1.5,
            cex.main = 10,
            ColSideLabs = ColsideLabels,
            Colv= as.dendrogram(sHc),
            dendrogram = "column",
            #Colv = NULL,
            Rowv = NULL,
            labCol=NA,
            side.height.fraction = 0.3,
            cexRow = 1.3,
            margins = c(13, 30))
  
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                          "Oncogenic pathway enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
  legend("bottomleft",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
  legend("bottomright", legend = Legend2,
         col = Legend_colors2, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
  legend("topright", legend = Legend3,
         col = Legend_colors3, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
  
  dev.off()
  
  Hallmark_and_ICR_cluster_assignment_allcancers[, paste0(Pathway, "_cluster_Pancancer")] = annotation$pathway_cluster_name[match(rownames(Hallmark_and_ICR_cluster_assignment_allcancers),
                                                                                                                                  rownames(annotation.blot))]
}

save(Cancer_color_table, Hallmark_and_ICR_cluster_assignment_allcancers, file = paste0("./4_Analysis/", download.method,
                                                                                       "/Pan_Cancer/Clustering/Hallmark_and_ICR_cluster_assignment_allcancers.Rdata"))
