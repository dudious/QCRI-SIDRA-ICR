#################################################################
###
### This script creates a correlation matrix for each cancer
### to find correlations between pathways within specific ICR
### subgroups
###
### Data input:
### "./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/
### "_gene_RNAseq_normalized_TP_filtered_ICR_subsetted.Rdata"
### Output :
### "./5_Figures/Correlation_plots/Candidate_pathways_low_ICR_Correlation_plots/", download.method, 
### "/Candidate_pathways_low_ICR_Correlation_plot_",Cancer,".png"
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                     # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                           # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 

required.packages <- c("corrplot", "stringr", "GSVA")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "BRCA"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/3.10_Correlation_matrix_Candidate_Low_ICR/Correlation_matrix_Candidate_Low_ICR_Log_File_",
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"
version = "v1"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load (paste0(code_path, "Datalists/Candidate_pathways_low_ICR.Rdata"))
load (paste0(code_path, "Datalists/ICR_genes.RData"))
Candidate_List$ICR_genes = ICR_genes

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Candidate_pathways_low_ICR_Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Candidate_pathways_low_ICR_Correlation_plots/",
                  download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/3.10_Correlation_matrix_Candidate_Low_ICR"), showWarnings = FALSE)
cat("This is a log file for creating correlation plots",
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    file = Log_file,
    append = FALSE, sep= "\n")

N.sets = length(CancerTYPES)

# Make correlation plots
start.time <- Sys.time ()

i=1
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  ## load RNASeq data
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered_ICR_subsetted.Rdata"))
    
  #ssGSEA
  Candidate.pathways.enrichment.score = gsva(Med_Low_ICR_RNAseqData, Candidate_List,method="ssgsea")
  Candidate.pathways.enrichment.z.score = Candidate.pathways.enrichment.score 
  for(j in 1: nrow(Candidate.pathways.enrichment.z.score))  {
    Candidate.pathways.enrichment.z.score[j,] = (Candidate.pathways.enrichment.score[j,]-mean(Candidate.pathways.enrichment.score[j,]))/sd(Candidate.pathways.enrichment.score[j,]) # z-score the enrichment matrix
  }
  save(Candidate.pathways.enrichment.score, Candidate.pathways.enrichment.z.score, 
       file = paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, "_", version, 
                                                     "_Candidate.pathways.low.ICR.Rdata"))
  
  Candidate.pathways.enrichment.z.score.df = as.data.frame(t(Candidate.pathways.enrichment.z.score))
  
  #Correlation matrix
  Candidate_cor <- cor (Candidate.pathways.enrichment.z.score.df, method=test)
  
  # Correlation significance
  cor.mtest <- function(mat, conf.level = 0.95) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
        uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
      }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
  }
  Candidate_cor_sign <- cor.mtest(Candidate_cor, 0.95)
  
  # Correlation plot
  png(paste0("./5_Figures/Correlation_plots/Candidate_pathways_low_ICR_Correlation_plots/", download.method, 
             "/Candidate_pathways_low_ICR_Correlation_plot_",Cancer,".png"), res=600,height=6.5,width=6,unit="in")
  cex.before <- par("cex")
  par(cex = 0.45)
  lims=c(-1,1)
  if (length(Candidate_cor[Candidate_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(Candidate_cor),color = c(rep("#CC0506",26)),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(Candidate_cor,order="FPC"),]
  
  mean_correlation = round(mean(Candidate_cor),2)
  corrplot.mixed (Candidate_cor,
                  #type="lower",
                  #p.mat = Candidate_cor_sign[[1]],                                                                            # add significance to correlations
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
                  tl.cex = 1,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(6,4.1,7,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between pathways upregulated in ICR low. \n ","Mean: ", mean_correlation,
                           ".\n Number of patients: ", nrow(Candidate.pathways.enrichment.z.score.df), "."),
                    cex = 1.5), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.65)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  
  dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation"), showWarnings = FALSE)
  save(Candidate_cor, Candidate_cor_sign, file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/",
                                            "Correlation_matrix_Candidate_Low_ICR_", test, "_", Cancer, ".Rdata"))
}

end.time <- Sys.time ()
time <- end.time - start.time
print (time)
