#################################################################
###
### This script creates heatmaps for ICR genes using EDASeq 
### normalized RNASeq data from the TCGA. This script includes a 
### log transformation of the data. Samples are ordered by ICR 
### score and a heatmap is created for ICR genes including 
### ICR cluster allocation for k = 3, 4 and 6 and the proposed
### HML classfication.
### Heatmaps are saved as png files at location: 
### ./5_Figures/Heatmaps/ICR_Heatmaps/", download.method, "/ICR_Heatmap_RNASeq_",Cancer/
###
#################################################################

## Create Heatmap of TCGA RNASeq ICR genes
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/heatmap.3.R"))

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper", "heatmap3", "plyr", "spatstat",
                      "heatmap3")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/3.4_Heatmap_ICR/3.4_Heatmap_ICR_Log_File_",                                            # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
version = "v3"                                                                                                          # Specify version of manual correction to perform in this script

# Load data
load (paste0(code_path, "Datalists/ICR_genes.RData"))
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
                                                                                                                        # in the Manual of Assembler v2.0.3 and was saved as csv file.
                                                                                                                
cluster_assignment_analysis = read.csv(paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Clustering/Cluster_assignment_analysis",
                                              "_manual_corrected_", version, ".csv"),stringsAsFactors = FALSE)

# Create folders and log file
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.4_Heatmap_ICR"), showWarnings = FALSE)
cat("This is a log file for creating heatmaps of RNASeq data based on ICR genes",                                       # Set-up logfile
    "__________________________________________",
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
    paste0("version of manual correction file for HML allocation used to generate heatmaps = ", version),
    "",
    "Scripts output :",
    "",
    "Creating heatmaps",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create heatmaps
start.time.process.all = Sys.time()
msg = paste0("Create heatmaps", "\n")
cat(msg)

i=26
for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Creating heatmap ",Cancer,"."))
  
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
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  ICR_subset_RNAseq = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% ICR_genes, ])
  ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)
  
  dir.create(paste0("./5_Figures/"),showWarnings = FALSE)
  dir.create(paste0("./5_Figures/Heatmaps/"),showWarnings = FALSE)
  dir.create(paste0("./5_Figures/Heatmaps/ICR_Heatmaps"),showWarnings = FALSE)
  dir.create(paste0("./5_Figures/Heatmaps/ICR_Heatmaps/",download.method),showWarnings = FALSE)
  
  ## plot annotation
  table_cluster_assignment = table_cluster_assignment[order(table_cluster_assignment$ICRscore),]
  annotation = table_cluster_assignment[,c(1,6,8,10,13)]
  
  levels(annotation$ICR_cluster_k3) = c("ICR Low", "ICR Medium1", "ICR High")
  levels(annotation$ICR_cluster_k4) = c("ICR Low", "ICR Medium1", "ICR Medium2", "ICR High")
  levels(annotation$ICR_cluster_k5) = c("ICR Low", "ICR Medium1", "ICR Medium2", "ICR Medium3", "ICR High")
  
  annotation.blot = as.matrix(annotation[,-1])
  annotation.blot[annotation.blot=="ICR Low"]="blue"
  annotation.blot[annotation.blot=="ICR Medium1"]="green"
  annotation.blot[annotation.blot=="ICR Medium2"]="yellow"
  annotation.blot[annotation.blot=="ICR Medium3"]="orange"
  annotation.blot[annotation.blot=="ICR High"]="red"
  annotation.blot[annotation.blot=="ICR Medium"]="green"
  
  my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
  expression.matrix = t(ICR_subset_RNAseq_log2)[,rownames(annotation.blot)]
  
  HML_actions = paste0("Optimal k for clustering was ", optimal.calinsky, ". \nSelected number of clusters for HML classification was: ",
                       cluster_assignment_analysis$chosen.K[cluster_assignment_analysis$Cancertype == Cancer],
                       ".\nThe combination action performed to reduce to 3 clusters is: \n",
                       cluster_assignment_analysis$manual.action[cluster_assignment_analysis$Cancertype == Cancer],
                       " \nResulting in a mean scaled ICR scores of ",
                       cluster_assignment_analysis$mean.scaled.ICR.scores.HML[cluster_assignment_analysis$Cancertype == Cancer],
                       " for Low / Medium / High resp.")
  
  Delta = cluster_assignment_analysis[which(cluster_assignment_analysis$Cancertype == Cancer), "delta.scaled.HL.in.HML.classification"]
  
  ## plot heatmap
  png(paste0("./5_Figures/Heatmaps/ICR_Heatmaps/", download.method, "/ICR_Heatmap3_RNASeq_",Cancer,".png"),res=600,height=8,width=11,unit="in")
  heatmap3 (expression.matrix,
            col=my.palette,
            Colv=NA,
            breaks=unique(c(seq(-4,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100))),
            scale = "row",
            legendfun=function() showLegend(legend=c("Low expression","High expression"),col=c("blue","red"),cex=1),
            ColSideColors= annotation.blot,
            ColSideLabs = c(paste0("k3 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k3")$freq))),
                            paste0("k4 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k4")$freq))),
                            paste0("k5 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k5")$freq))),
                            paste0("Proposed HML clusters: ", gsub(","," /",toString(count(annotation, "HML_cluster")$freq[c(2,3,1)])))),
            cexCol = 0.9,
            labCol=FALSE,
            cexRow = 1.5,
            margins=c(3,3))
  
  title(main = list(paste0(Cancer, ": RNASeq", "\n N patients = ", nrow(table_cluster_assignment), "\n Delta = ", Delta),cex = 1.5),
        outer = FALSE, line = -2, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was obtained from TCGA, \n using ", download.method, " v2.0.3. ",
                          "Samples are ordered by ICR score.\n", HML_actions), cex = 0.85), outer = FALSE, line = -0.75)
  legend(x = 0 , y= 0.9, legend = c("ICR High", "ICR Medium3", "ICR Medium2", "ICR Medium1", "ICR Low"),
         col = c("red","orange","yellow", "green","blue"), lty= 0.5,lwd = 0.5, cex = 0.9, pch= 15, pt.cex = 1.2)
  dev.off()
  
  png(paste0("./5_Figures/Heatmaps/ICR_Heatmaps/", download.method, "/ICR_Heatmap.3_RNASeq_",Cancer,".png"),res=600,height=9.5,width=11,unit="in")
  heatmap.3 (expression.matrix,
             col=my.palette,
             Colv=NA,
             scale = "row",
             side.height.fraction = 0.5,
             ColSideColors= annotation.blot,
             ColSideLabs = c(paste0("k3 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k3")$freq))),
                             paste0("k4 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k4")$freq))),
                             paste0("k5 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k5")$freq))),
                             paste0("Proposed HML clusters: ", gsub(","," /",toString(count(annotation, "HML_cluster")$freq[c(2,3,1)])))),
             font_size_col_Labs = 1.3,
             cexCol = 0.9,
             labCol=FALSE,
             cexRow = 1.5,
             margins=c(15,15))
  
  title(main = list(paste0(Cancer, ": RNASeq", "\n N patients = ", nrow(table_cluster_assignment), "\n Delta = ", Delta),cex = 1.5),
        outer = FALSE, line = -4 , adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was obtained from TCGA, \n using ", download.method, " v2.0.3. ",
                          "Samples are ordered by ICR score. ", HML_actions), cex = 0.80), outer = FALSE, line = -0.82, adj = 0.55)
  legend("topright", legend = c("ICR High", "ICR Medium3", "ICR Medium2", "ICR Medium1", "ICR Low"),
         col = c("red","orange","yellow", "green","blue"), lty= 0.5,lwd = 0.5, cex = 0.9, pch= 15, pt.cex = 1.2)
  dev.off()
}

