#################################################################
###
<<<<<<< HEAD
### Create ICR RNASEQ based expression Heatmap....
###
#################################################################

=======
### This script creates heatmaps for ICR genes using EDASeq 
### normalized RNASeq data from the TCGA. This script includes a 
### log transformation of the data. Samples are ordered by ICR 
### score and a heatmap is created for ICR genes including 
### ICR cluster allocation for k = 3, 4 and 6.
### Heatmaps are saved as png files at location: 
### ./5_Figures/Heatmaps/ICR_Heatmaps/", download.method, "/ICR_Heatmap_RNASeq_",Cancer/
###
#################################################################

## Create Heatmap of TCGA RNASeq ICR genes

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
>>>>>>> 1172c97ced4e87aeaf07a1c6fc003ad99a2350b9
# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")
#setwd("D:/Jessica/Dropbox (TBI-Lab)/TCGA Analysis pipeline/") 
required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper", "heatmap3", "plyr", "spatstat")
required.bioconductor.packages = c("heatmap3")

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)

Path.R.Tools = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/R tools/"                               # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
#Path.R.Tools = "D:/Jessica/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/R tools/"
Path.Pipeline.Scripts = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/Any Cancer (V4)/"
#Path.Pipeline.Scripts = "D:/Jessica/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/Any Cancer (V4)/"

Log_file = paste0("./1_Log_Files/3.2_Heatmap_ICR/3.2_Heatmap_ICR_Log_File_",                                            # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data and R scripts
TCGA.cancersets = read.csv ("./TCGA.datasets.csv",stringsAsFactors = FALSE)                                             # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.

source(paste0(Path.R.Tools, "ipak.function.R"))
source(paste0(Path.Pipeline.Scripts, "0.1.Specification_ICR_genes_for_pipeline.R"))


#Install and load required packages
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Create folders and log file
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.2_Heatmap_ICR"), showWarnings = FALSE)
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

for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Creating heatmap ",Cancer,"."))
  
  ## load RNASeq data
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
  annotation = table_cluster_assignment[,c(1,5,7,9)]
  
  levels(annotation$ICR_cluster_k3) = c("Low ICR", "Medium1 ICR", "High ICR")
  levels(annotation$ICR_cluster_k4) = c("Low ICR", "Medium1 ICR", "Medium2 ICR", "High ICR")
  levels(annotation$ICR_cluster_k5) = c("Low ICR", "Medium1 ICR", "Medium2 ICR", "Medium3 ICR", "High ICR")
  
  annotation.blot = as.matrix(annotation[,-1])
  annotation.blot[annotation.blot=="Low ICR"]="blue"
  annotation.blot[annotation.blot=="Medium1 ICR"]="green"
  annotation.blot[annotation.blot=="Medium2 ICR"]="yellow"
  annotation.blot[annotation.blot=="Medium3 ICR"]="orange"
  annotation.blot[annotation.blot=="High ICR"]="red"
  
  my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 300)
  expression.matrix = t(ICR_subset_RNAseq_log2)[,rownames(annotation.blot)]
  
  ## plot heatmap
  #dev.new()
  png(paste0("./5_Figures/Heatmaps/ICR_Heatmaps/", download.method, "/ICR_Heatmap_RNASeq_",Cancer,".png"),res=600,height=7,width=10,unit="in")
  heatmap3 (expression.matrix,
            col=my.palette,
            Colv=NA,
            ColSideColors= annotation.blot,
            ColSideLabs = c(paste0("k3 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k3")$freq))),
                            paste0("k4 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k4")$freq))),
                            paste0("k5 ICR clusters: ", gsub(","," /",toString(count(annotation, "ICR_cluster_k5")$freq)))),
            cexCol = 0.9,
            labCol=FALSE,
            cexRow = 1.5,
            margins=c(3,3))
  
  title(main = list(paste0(Cancer, ": RNASeq", "\n N patients = ", nrow(table_cluster_assignment), "\n Optimal k = ", optimal.calinsky),cex = 1.5),
        outer = FALSE, line = -2, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was obtained from TCGA, \n using ", download.method, " v2.0.3. ",
                          "Samples are ordered by ICR score."), cex = 0.8), outer = FALSE, line = 4)
  legend(x = 0 , y= 0.9, legend = c("High ICR", "Medium3 ICR", "Medium2 ICR", "Medium1 ICR", "Low ICR"),
         col = c("red","orange","yellow", "green","blue"), lty= 0.5,lwd = 0.5, cex = 0.9, pch= 15, pt.cex = 1.2)
  dev.off()
}

  
