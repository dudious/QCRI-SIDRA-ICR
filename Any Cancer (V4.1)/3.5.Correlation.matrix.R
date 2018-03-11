#################################################################
###
### This script creates a correlation matrix for the selected genes
###
### Data input:
### "./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/
### "_gene_RNAseq_normalized_TP_filtered.Rdata"
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                     # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                           # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 

required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Assembler_Panca_Normalized"                                                                                       # Specify download method (this information to be used when saving the file)
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)
selected_genes = "ICR"                                                                                                   # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/", download.method ,"/3.5_Correlation_matrix/Correlation_matrix_", selected_genes, "_Log_File_",
                  gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}
if (selected_genes == "ICR") {
  load(paste0(code_path, "Datalists/ICR_genes.RData"))
  genes_to_correlate = ICR_genes
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/", selected_genes, "_Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/", selected_genes, "_Correlation_plots/",
                  download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method ,"/3.5_Correlation_matrix"), showWarnings = FALSE)
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

mean_correlation_table = data.frame(Cancertype = CancerTYPES, Mean.correlation = 0)

i=1
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  ## load RNASeq data
  if(Cancer == "LAML") 
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
  
  # Subset RNAseq data to genes to correlate
  unavailable.genes <- genes_to_correlate[-which(genes_to_correlate %in% rownames(filtered.norm.RNAseqData))]
  subset_RNAseq = t(filtered.norm.RNAseqData[row.names(filtered.norm.RNAseqData) %in% genes_to_correlate, ])
  if(download.method == "TCGA_Assembler" | download.method == "Assembler_Panca_Normalized"){
    subset_RNAseq_log2 = log(subset_RNAseq +1, 2) 
  }
  if(download.method == "Pancancer_matrix"){
    subset_RNAseq_log2 = subset_RNAseq
  }
  
  # Correlation matrix
  ICR_cor <- cor (subset_RNAseq_log2,method=test)
  
  # Correlation significance
  cor.mtest <- function(mat, conf.level = 0.95) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = test, conf.level = conf.level)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
        uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
      }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
  }
  ICR_cor_sign <- cor.mtest(subset_RNAseq_log2, 0.95)
  
  # Correlation plot
  png(paste0("./5_Figures/Correlation_plots/", selected_genes, "_Correlation_plots/",
             download.method, "/", selected_genes, "_", test, "_Correlation_plot_", Cancer, ".png"),
      res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.45)
  lims=c(-1,1)
  if (length(ICR_cor[ICR_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(ICR_cor),color = c(rep("#CC0506",20)),stringsAsFactors = FALSE)
  annotation$color[annotation$gene %in% c("IDO1","CD274","CTLA4","FOXP3","PDCD1")] = "#41719C"
  annotation = annotation[corrMatOrder(ICR_cor,order="FPC"),]
  
  mean_correlation = round(mean(ICR_cor),2)
  corrplot.mixed (ICR_cor,
                  #type="lower",
                  #p.mat = ICR_cor_sign[[1]],                                                                      # add significance to correlations
                  lower.col = colpattern,
                  upper.col = colpattern,
                  lower = "square",
                  upper ="number",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.5,
                  tl.cex = 1/par("cex"),
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(6,4.1,7,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between ", selected_genes, " genes. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(subset_RNAseq), "."),
                    cex = 2.2), line = -2.5)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, ". \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 0)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation"), showWarnings = FALSE)
  save(ICR_cor, ICR_cor_sign, file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/",
                                            "Correlation_matrix_ICR_", test, "_", Cancer, ".Rdata"))
}

dir.create("./4_Analysis/", download.method, "/Pan_Cancer/Correlation", showWarnings = FALSE)
save(mean_correlation_table, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Correlation/Correlation_", selected_genes, ".Rdata"))
end.time <- Sys.time ()
time <- end.time - start.time
print (time)
