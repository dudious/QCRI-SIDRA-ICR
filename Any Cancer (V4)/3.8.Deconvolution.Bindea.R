#################################################################
###
### Create ICR RNASEQ based Deconvolution
### using Bindea's signatures and ssGSEA
###
### Input files:
### ./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData
### "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"
### Output files:
###
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path,"R tools/ipak.function.R"))

required.packages = c("devtools")                                                                                      
required.bioconductor.packages = c("GSVA","heatmap3")                                                                   
ipak(required.packages)
ibiopak(required.bioconductor.packages)
devtools::install_github('dviraran/xCell')

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/3.8_Deconvolution_Bindea/3.8_Deconvolution_Bindea_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load (paste0(code_path, "Datalists/marker_list_bindea2.RData"))
load (paste0(code_path, "Datalists/ICR_genes.RData"))

# Create folders and log file
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.8_Deconvolution_Bindea"), showWarnings = FALSE)
cat("This is a log file for Deconvolution using Bindeas gene signatures on RNASeq data",                                # Set-up logfile
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

# Append ICR genes to bindea
marker_list$ICR_score <- ICR_genes

start.time.process.all = Sys.time()
msg = paste0("Calculating deconvolution scores and generating heatmaps", "\n")
cat(msg)

i=1
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
    load(paste0(Cancer_path, Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  
  ## Log transformation of normalized data
  ICR_subset_RNAseq = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% ICR_genes, ])
  ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2)
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  ## Annotation for plotting
  annotation <- Clustering[,"ICR_HML",drop=FALSE]
  annotation$ICR_HML.col[annotation$ICR_HML=="ICR_high"]<-"red"
  annotation$ICR_HML.col[annotation$ICR_HML=="ICR_medium"]<-"green"
  annotation$ICR_HML.col[annotation$ICR_HML=="ICR_low"]<-"blue"
  annotation.blot <- as.matrix(annotation[,c("annotation$ICR_HML.col")])
  
  ## Bindea ssGSEA
  Bindea.enrichment.score <- gsva(Expression.data,marker_list,method="ssgsea")
  Bindea.enrichment.z.core <- Bindea.enrichment.score 
  for(j in 1: nrow(Bindea.enrichment.score))
  {
    Bindea.enrichment.z.score[j,]<- (Bindea.enrichment.score[j,]-mean(Bindea.enrichment.score[j,]))/sd(Bindea.enrichment.score[j,]) # z-score the inrichment matrix
  }
  ### Bindea plotting
  dev.new()
  Bindea.enrichment.z.score <- Bindea.enrichment.z.score[,rownames(annotation.blot)]
  heatmap3((as.matrix(Bindea.enrichment.z.core)),
           main="ssGSEA/bindea signatures",
           ColSideColors=annotation.blot,
           Colv= NA, #as.dendrogram(sHc),
           #labCol=colnames(enrichment.z.scores),
           margins = c(12, 7))
  legend("topright",legend = c("ICR_low","ICR_med","ICR_high"),
         col = c("blue","green","red"),lty= 1,lwd = 5,cex = 0.7)
  
  ## xCell
  cat (paste0 ("Calculating xCell scores ",Cancer,"."))
  xCell.matrix <- xCellAnalysis(expression.data)
  
  ### xCell Forrest plot
  xCell.matrix<-xCell.matrix[,rownames(annotation.blot)]
  dev.new()
  heatmap.3(xCell.matrix,
            main = paste0("xCell Heatmap immune enrichment"),
            #col=my.palette,
            #breaks = my.colors,
            scale = "row",
            Colv=NA,
            #labCol=samples,
            ColSideColors=annotation.blot,
            cexRow=1.3,cexCol=1,margins=c(8,15)
  )
  par(lend = 1)
  legend("topright",legend = c("ICR_low","ICR_med","ICR_high"),
         col = c("blue","green","red"),lty= 1,lwd = 5,cex = 0.7)
  
  ## Save Scores
  save (Bindea.enrichment.score,Bindea.enrichment.z.core,xCell.matrix,
        file = "./4_Analysis/", download.method, "/", Cancer, "/Signature enrichment/", Cancer, ".", download.method, ".ssGSEA.xCell.Rdata")
}
  
