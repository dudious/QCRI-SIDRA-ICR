#################################################################
###
### Create RNASEQ based Deconvolution
### using Bindea's signatures and ssGSEA
###
### Input files:
### ./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData
### "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"
### Output files:
### 
###
#################################################################

#### Xcell is not included in the log file !!!!

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"                                                                # Setwd to location were output files have to be saved.

source(paste0(code_path,"R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

if(!("xCell" %in% installed.packages()[,"Package"])){
  required.packages = c("devtools")
  ipak(required.packages)
  devtools::install_github('dviraran/xCell')
}
required.bioconductor.packages = c("GSVA","heatmap3", "gclus")                                                                   
ibiopak(required.bioconductor.packages)
library ("xCell")

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Pancancer_matrix"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
pw_selection_version = "3.2"
Log_file = paste0("./1_Log_Files/3.8_Deconvolution_Bindea/3.8_Deconvolution_Bindea_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("HML ICR clusters", "Bindea clusters")
Legend = c("ICR Low","ICR Med","ICR High", "Bindea Low", "Bindea High")
Legend_colors = c("blue","green","red", "pink", "purple")

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0(code_path, "Datalists/Selected.pathways.",pw_selection_version,".Rdata"))
load(paste0(code_path, "Datalists/marker_list_bindea3.RData"))

# Create folders and log file
dir.create(paste0("./5_Figures/"),showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/"),showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/Bindea_Heatmaps"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/Bindea_Heatmaps/", download.method),showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/bindea_patrick_Heatmaps"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/bindea_patrick_Heatmaps/", download.method),showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/xCell_Heatmaps/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/xCell_Heatmaps/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/Hallmark_Heatmaps"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/Hallmark_Heatmaps/", download.method),showWarnings = FALSE)
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folders to save Rdata.files
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.8_Deconvolution_Bindea"), showWarnings = FALSE)
cat("This is a log file for Deconvolution using Bindeas gene signatures, xCell and Hallmark pathways on RNASeq data",   # Set-up logfile
    "_________________________________________________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),                                                          
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    paste0("pathway selection version used =", pw_selection_version),
    "",
    "Scripts output :",
    "",
    "Calculating Deconvolution scores",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

start.time.process.all = Sys.time()
msg = paste0("Calculating deconvolution scores and generating heatmaps", "\n")
cat(msg)


for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Calculating Deconvolution scores ",Cancer,"."))
  
  ## load RNASeq data
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  if(Cancer == "SKCM"){
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  Expression.data = log(filtered.norm.RNAseqData +1, 2)
  available_genes = rownames(Expression.data)
  unavailable_genes_RNAseq = unlist(Bindea_ORIG)[-which(unlist(Bindea_ORIG) %in% rownames(Expression.data))]
  
  cat(paste0("Bindea ssGSEA ", Cancer, ". Total number of Bindea genes is ", length(unlist(Bindea_ORIG)), ".",
             " Of which ", length(unlist(Bindea_ORIG)[unlist(Bindea_ORIG) %in% available_genes]), 
             " genes are available in expression data."), file = Log_file, append = TRUE, sep = "\n")
  
  ## Bindea ssGSEA
  Bindea.enrichment.score = gsva(Expression.data,Bindea_ORIG,method="ssgsea")
  Bindea.enrichment.z.score = Bindea.enrichment.score 
  for(j in 1: nrow(Bindea.enrichment.z.score))  {
    Bindea.enrichment.z.score[j,] = (Bindea.enrichment.score[j,]-mean(Bindea.enrichment.score[j,]))/sd(Bindea.enrichment.score[j,]) # z-score the enrichment matrix
  }
  #Bindea.enrichment.ordered = Bindea.enrichment.z.score
  #Bindea.enrichment.ordered = rbind(Bindea.enrichment.ordered, colMeans(Bindea.enrichment.ordered))
  #rownames(Bindea.enrichment.ordered)[25] = "Bindea_enrichment_score"
  #Bindea.enrichment.ordered = Bindea.enrichment.ordered[, order(Bindea.enrichment.ordered[which(rownames(Bindea.enrichment.ordered) == "Bindea_enrichment_score"),])]
  #Bindea.enrichment.z.score = Bindea.enrichment.z.score[,colnames(Bindea.enrichment.ordered)]
  
  ## Bindea cluster (hierarchical)
  sHc = hclust(ddist <- dist(t(Bindea.enrichment.z.score)), method = "ward.D2")
  
  plot(sHc,labels=FALSE)
  
  ## Annotation for plotting
  annotation <- table_cluster_assignment[,"HML_cluster",drop=FALSE]
  # Bindea classification
  annotation$bindea_cluster = cutree(sHc,k = 2)[match(rownames(annotation),names(cutree(sHc,k = 2)))]
  annotation = annotation[colnames(Bindea.enrichment.score),]
  annotation$bindea_score = colMeans(Bindea.enrichment.score)
  bindea_cluster_means = aggregate(bindea_score~bindea_cluster,data=annotation,FUN=mean)
  bindea_cluster_means = bindea_cluster_means[order(bindea_cluster_means$bindea_score),]
  bindea_cluster_means$bindea_cluster_name = c("Bindea Low","Bindea High")
  annotation$bindea_cluster_name = bindea_cluster_means$bindea_cluster_name[match(annotation$bindea_cluster,bindea_cluster_means$bindea_cluster)]
  
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"
  annotation$Bindea_cluster.col[annotation$bindea_cluster_name=="Bindea Low"] = "pink"
  annotation$Bindea_cluster.col[annotation$bindea_cluster_name=="Bindea High"] = "purple"
  Bindea_order = c("Bindea Low", "Bindea High")
  ICR_order = c("ICR Low","ICR Medium","ICR High")
  annotation = annotation[, -which(colnames(annotation) %in% c("bindea_cluster", "bindea_score"))]
  annotation = annotation[order(match(annotation$HML_cluster,ICR_order), match(annotation$bindea_cluster_name, Bindea_order)),]
  annotation.blot = as.matrix(annotation[,c("HML_cluster.col","Bindea_cluster.col"), drop= FALSE])
  annotation.blot = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
  
  ### Bindea plotting
  png(paste0("./5_Figures/Heatmaps/Bindea_Heatmaps/", download.method, "/Bindea_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height=9,width=9,unit="in")
  Bindea.enrichment.z.score <- Bindea.enrichment.z.score[,rownames(annotation.blot)]
  heatmap.3((as.matrix(Bindea.enrichment.z.score)),
            main= paste0(Cancer, "\nssGSEA/Bindea signatures"),
            col=my.palette,
            ColSideColors=annotation.blot,
            font_size_col_Labs = 1.5,
            cex.main = 10,
            ColSideLabs = ColsideLabels,
            Colv= as.dendrogram(sHc),
            labCol=NA,
            side.height.fraction = 0.3,
            cexRow = 1.3,
            margins = c(13, 11))
  
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, " v2.0.3.\n",
                          "Bindea enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
  
  dev.off()
  
  ## Patrick bindea_patrick ssGSEA
  bindea_patrick.enrichment.score = gsva(Expression.data,Bindea_PATRICK,method="ssgsea")
  bindea_patrick.enrichment.z.score = bindea_patrick.enrichment.score 
  for(j in 1: nrow(bindea_patrick.enrichment.z.score))  {
    bindea_patrick.enrichment.z.score[j,] = (bindea_patrick.enrichment.score[j,]-mean(bindea_patrick.enrichment.score[j,]))/sd(bindea_patrick.enrichment.score[j,]) # z-score the enrichment matrix
  }
  #bindea_patrick.enrichment.ordered = bindea_patrick.enrichment.z.score
  #bindea_patrick.enrichment.ordered = rbind(bindea_patrick.enrichment.ordered, colMeans(bindea_patrick.enrichment.ordered))
  #rownames(bindea_patrick.enrichment.ordered)[15] = "bindea_patrick_enrichment_score"
  #bindea_patrick.enrichment.ordered = bindea_patrick.enrichment.ordered[, order(bindea_patrick.enrichment.ordered[which(rownames(bindea_patrick.enrichment.ordered) == "bindea_patrick_enrichment_score"),])]
  #bindea_patrick.enrichment.z.score = bindea_patrick.enrichment.z.score[,colnames(bindea_patrick.enrichment.ordered)]
  
  ## bindea_patrick cluster (hierarchical)
  sHc = hclust(ddist <- dist(t(bindea_patrick.enrichment.z.score)), method = "ward.D2")
  
  plot(sHc,labels=FALSE)
  
  ## Annotation for plotting
  annotation <- table_cluster_assignment[,"HML_cluster",drop=FALSE]
  # bindea_patrick classification
  annotation$bindea_patrick_cluster = cutree(sHc,k = 2)[match(rownames(annotation),names(cutree(sHc,k = 2)))]
  annotation = annotation[colnames(bindea_patrick.enrichment.score),]
  annotation$bindea_patrick_score = colMeans(bindea_patrick.enrichment.score)
  bindea_patrick_cluster_means = aggregate(bindea_patrick_score~bindea_patrick_cluster,data=annotation,FUN=mean)
  bindea_patrick_cluster_means = bindea_patrick_cluster_means[order(bindea_patrick_cluster_means$bindea_patrick_score),]
  bindea_patrick_cluster_means$bindea_patrick_cluster_name = c("Bindea Patrick Low","Bindea Patrick High")
  annotation$bindea_patrick_cluster_name = bindea_patrick_cluster_means$bindea_patrick_cluster_name[match(annotation$bindea_patrick_cluster,bindea_patrick_cluster_means$bindea_patrick_cluster)]
  
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"
  annotation$bindea_patrick_cluster.col[annotation$bindea_patrick_cluster_name=="Bindea Patrick Low"] = "pink"
  annotation$bindea_patrick_cluster.col[annotation$bindea_patrick_cluster_name=="Bindea Patrick High"] = "purple"
  bindea_patrick_order = c("Bindea Patrick Low", "Bindea Patrick High")
  ICR_order = c("ICR Low","ICR Medium","ICR High")
  annotation = annotation[, -which(colnames(annotation) %in% c("bindea_patrick_cluster", "bindea_patrick_score"))]
  annotation = annotation[order(match(annotation$HML_cluster,ICR_order), match(annotation$bindea_patrick_cluster_name, bindea_patrick_order)),]
  annotation.blot2 = as.matrix(annotation[,c("HML_cluster.col","bindea_patrick_cluster.col"), drop= FALSE])
  annotation.blot2 = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
  
  ### bindea_patrick plotting
  png(paste0("./5_Figures/Heatmaps/bindea_patrick_Heatmaps/", download.method, "/bindea_patrick_Patrick_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height=9,width=9,unit="in")
  bindea_patrick.enrichment.z.score <- bindea_patrick.enrichment.z.score[,rownames(annotation.blot2)]
  heatmap.3((as.matrix(bindea_patrick.enrichment.z.score)),
            main= paste0(Cancer, "\nssGSEA/ Bindea signatures by Patrick Danaher et al (2017)"),
            col=my.palette,
            ColSideColors=annotation.blot2,
            font_size_col_Labs = 1.5,
            cex.main = 10,
            ColSideLabs = ColsideLabels,
            Colv= as.dendrogram(sHc),
            labCol=NA,
            side.height.fraction = 0.3,
            cexRow = 1.3,
            margins = c(13, 11))
  
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, " v2.0.3.\n",
                          "bindea_patrick enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
  
  dev.off()
  
  ## xCell
  cat (paste0 ("Calculating xCell scores ",Cancer,"."))
  xCell.matrix <- xCellAnalysis(Expression.data)
  
  ### xCell heatmap plot
  annotation.blot = annotation.blot[rownames(annotation),]
  xCell.matrix = xCell.matrix[,rownames(annotation.blot)]
  
  png(paste0("./5_Figures/Heatmaps/xCell_Heatmaps/", download.method, "/xCell_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height =12.5,width=9,unit="in")
  heatmap.3(xCell.matrix,
            main = paste0(Cancer, "\nxCell Heatmap"),
            cex.main = 2,
            col= my.palette,
            ColSideLabs = ColsideLabels,
            font_size_col_Labs = 1.5,
            ColSideColors= annotation.blot,
            #breaks = my.colors,
            scale = "row",
            Colv=NA,
            labCol=NA,
            side.height.fraction = 0.2,
            cexRow= 1,
            margins=c(13,13))
  par(lend = 1)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, " v2.0.3.\n",
                          "xCell R package was used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.2)
  dev.off()

  ## Hallmark ssGSEA
  # Check availability of Hallmark genes in the dataset
  available_genes_RNAseq = unlist(Selected.pathways)[which(unlist(Selected.pathways) %in% rownames(Expression.data))]
  unavailable_genes_RNAseq = unlist(Selected.pathways)[-which(unlist(Selected.pathways) %in% rownames(Expression.data))]
  cat(paste0("Hallmark ssGSEA ", Cancer, ". Total number of hallmark genes is ", length(unlist(Selected.pathways)), ".",
             " Of which ", length(unlist(Selected.pathways)[unlist(Selected.pathways) %in% available_genes]), 
             " genes are available in expression data.\n"),
      file = Log_file, append = TRUE, sep = "\n")
  
  Hallmark.enrichment.score = gsva(Expression.data, Selected.pathways, method="ssgsea")
  Hallmark.enrichment.z.score = Hallmark.enrichment.score 
  for(j in 1: nrow(Hallmark.enrichment.z.score))  {
    Hallmark.enrichment.z.score[j,] = (Hallmark.enrichment.score[j,]-mean(Hallmark.enrichment.score[j,]))/sd(Hallmark.enrichment.score[j,]) # z-score the enrichment matrix
  }
  
  annotation.blot = annotation.blot[rownames(annotation),]
  Hallmark.enrichment.z.score = Hallmark.enrichment.z.score[,rownames(annotation.blot)]
  
  png(paste0("./5_Figures/Heatmaps/Hallmark_Heatmaps/", download.method, "/HallmarkGSEA_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height =12.5,width=12.5,unit="in")
  heatmap.3(as.matrix(Hallmark.enrichment.z.score),
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
            side.height.fraction = 0.2,
            cexRow= 1,
            margins=c(13,25))
  par(lend = 1)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, " v2.0.3.\n",
                          "Hallmark enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.50)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.2)
  dev.off()
  
  ## Hallmark col means
  
  #z-score matrix
  Hallmark.col.means <- data.frame(matrix(nrow=length(Selected.pathways), ncol=ncol(Expression.data)))
  colnames(Hallmark.col.means) = colnames(Expression.data)
  rownames(Hallmark.col.means) = names(Selected.pathways)
  
 
  for(i in 1:length(Selected.pathways)){
    pathway_genes = unlist(Selected.pathways[i], use.names = FALSE)
    available_pathway_genes = pathway_genes[which(pathway_genes %in% rownames(Expression.data))]
    Hallmark.col.means[i,] = colMeans(Expression.data[available_pathway_genes,])
  }
  
  annotation.blot = annotation.blot[rownames(annotation),]
  Hallmark.col.means = as.matrix(Hallmark.col.means[,rownames(annotation.blot)])
  
  png(paste0("./5_Figures/Heatmaps/Hallmark_Heatmaps/", download.method, "/Hallmark_colmeans_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height =12.5,width=12.5,unit="in")
  heatmap.3(as.matrix(Hallmark.col.means),
            main = paste0(Cancer, "\n Hallmark col means Heatmap"),
            cex.main = 2,
            col= my.palette,
            ColSideLabs = ColsideLabels,
            font_size_col_Labs = 1.5,
            ColSideColors= annotation.blot,
            #breaks = my.colors,
            scale = "row",
            Colv=NA,
            labCol=NA,
            side.height.fraction = 0.2,
            cexRow= 1,
            margins=c(13,25))
  par(lend = 1)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, " v2.0.3.\n",
                          "Hallmark col means were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.50)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.2)
  dev.off()
  
  ## Save Scores
  dir.create(paste0("./4_Analysis/",download.method, "/", Cancer),showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment"))
  
  save(Bindea.enrichment.score, Bindea.enrichment.z.score, bindea_patrick.enrichment.score, bindea_patrick.enrichment.z.score,
       xCell.matrix, Hallmark.enrichment.score, Hallmark.enrichment.z.score, 
       Hallmark.col.means, annotation, file = paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
                                         "_Bindea_xCell_HallmarkPathways.Rdata"))
}

