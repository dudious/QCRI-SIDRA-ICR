#################################################################
###
### This script creates a correlation matrix for Bindea, xCell,
### Hallmark pathways cor
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"                                                                # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 

required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Pancancer_matrix"                                                                                       # Specify download method (this information to be used when saving the file)
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 200)                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/", download.method,"/3.9_Correlation_matrix_Signatures/3.9_Correlation_matrix_signatures", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Bindea_Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/bindea_patrick_Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/xCell_Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots/ssGSEA_correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots/colMeans_correlation_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Combined_Signature_Correlation_plots/"), showWarnings = FALSE)

dir.create(paste0("./5_Figures/Correlation_plots/Bindea_Correlation_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/bindea_patrick_Correlation_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/xCell_Correlation_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots/ssGSEA_correlation_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots/colMeans_correlation_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Correlation_plots/Combined_Signature_Correlation_plots/", download.method), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/3.9_Correlation_matrix_Signatures"), showWarnings = FALSE)
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

mean_correlation_table = data.frame(Cancertype = CancerTYPES, Mean.correlation.Bindea = 0, Mean.correlation.xCell = 0, Mean.correlation.Hallmark.CM = 0,
                                    Mean.correlation.Hallmark.GSEA = 0, Mean.correlation.bindea.xCell.Hallmark = 0)

i=1
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  ## load Signature data
  if(Cancer == "LAML") 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  
  load(paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
              "_Bindea_xCell_HallmarkPathways.Rdata"))   
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  ### Bindea
  # Correlation matrix Bindea z scores
  print(paste0("Generating Bindea z score correlation matrix for ", Cancer))
  
  Bindea.enrichment.z.score.df = as.data.frame(t(Bindea.enrichment.z.score))
  Bindea.enrichment.z.score.df$ICR_score = table_cluster_assignment$ICRscore[match(row.names(Bindea.enrichment.z.score.df), row.names(table_cluster_assignment))]
  
  Bindea_cor <- cor (Bindea.enrichment.z.score.df, method=test)
  
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
  Bindea_cor_sign = cor.mtest(Bindea.enrichment.z.score.df, 0.95)
  
  print(paste0("Generating Bindea z score correlation plot for ", Cancer))
  
  # Correlation plot
  png(paste0("./5_Figures/Correlation_plots/Bindea_Correlation_plots/", download.method,
             "/", "Bindea_", test, "_Correlation_plot_", Cancer, ".png"), res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.35)
  lims=c(-1,1)
  if (length(Bindea_cor[Bindea_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(Bindea_cor),color = c(rep("#CC0506",nrow(Bindea_cor))),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(Bindea_cor,order="FPC"),]
  
  mean_correlation = round(mean(Bindea_cor),2)
  corrplot.mixed (Bindea_cor,
                  #type="lower",
                  #p.mat = Bindea_cor_sign[[1]],                                                                            # add significance to correlations
                  lower.col = "black",
                  upper.col = colpattern,
                  lower = "number",
                  upper = "square",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.8,
                  tl.cex = 1.4,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(8,4.1,6,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between Bindea signatures. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(t(Bindea.enrichment.z.score)), "."),
                    cex = 2.2), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.55)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean bindea correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation.Bindea[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  ### bindea_patrick
  # Correlation matrix bindea_patrick z scores
  print(paste0("Generating bindea_patrick z score correlation matrix for ", Cancer))
  
  bindea_patrick.enrichment.z.score.df = as.data.frame(t(bindea_patrick.enrichment.z.score))
  bindea_patrick.enrichment.z.score.df$ICR_score = table_cluster_assignment$ICRscore[match(row.names(bindea_patrick.enrichment.z.score.df), row.names(table_cluster_assignment))]
  
  bindea_patrick_cor <- cor (bindea_patrick.enrichment.z.score.df, method=test)
  
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
  bindea_patrick_cor_sign = cor.mtest(bindea_patrick.enrichment.z.score.df, 0.95)
  
  print(paste0("Generating bindea_patrick z score correlation plot for ", Cancer))
  
  # Correlation plot
  png(paste0("./5_Figures/Correlation_plots/bindea_patrick_Correlation_plots/", download.method,
             "/", "bindea_patrick_", test, "_Correlation_plot_", Cancer, ".png"), res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.35)
  lims=c(-1,1)
  if (length(bindea_patrick_cor[bindea_patrick_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(bindea_patrick_cor),color = c(rep("#CC0506",nrow(bindea_patrick_cor))),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(bindea_patrick_cor,order="FPC"),]
  
  mean_correlation = round(mean(bindea_patrick_cor),2)
  corrplot.mixed (bindea_patrick_cor,
                  #type="lower",
                  #p.mat = bindea_patrick_cor_sign[[1]],                                                                            # add significance to correlations
                  lower.col = "black",
                  upper.col = colpattern,
                  lower = "number",
                  upper = "square",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.8,
                  tl.cex = 1.4,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(8,4.1,6,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between bindea_patrick signatures. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(t(bindea_patrick.enrichment.z.score)), "."),
                    cex = 2.2), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.55)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean bindea_patrick correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation.bindea_patrick[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  ### xCell
  print(paste0("Generating xCell correlation matrix for ", Cancer))
  
  xCell.df = as.data.frame(t(xCell.matrix))
  xCell.df$ICR_score = table_cluster_assignment$ICRscore[match(row.names(xCell.df), row.names(table_cluster_assignment))]
  
  # Correlation matrix xCell
  xCell_cor <- cor (xCell.df, method=test)
  
  # Correlation significance
  
  xCell_cor_sign = cor.mtest(xCell.df, 0.95)
  
  print(paste0("Generating xCell correlation plot for ", Cancer))
  # Correlation plot
  png(paste0("./5_Figures/Correlation_plots/xCell_Correlation_plots/", download.method,
             "/", "xCell_", test, "_Correlation_plot_", Cancer, ".png"), res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.35)
  lims=c(-1,1)
  if (length(xCell_cor[xCell_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(xCell_cor),color = c(rep("#CC0506",nrow(xCell_cor))),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(xCell_cor,order="FPC"),]
  
  mean_correlation = round(mean(xCell_cor),2)
  corrplot.mixed (xCell_cor,
                  #type="lower",
                  #p.mat = xCell_cor_sign[[1]],                                                                            # add significance to correlations
                  lower.col = "black",
                  upper.col = colpattern,
                  lower = "number",
                  upper = "square",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.8,
                  tl.cex = 0.8,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(8,4.1,6,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between xCell signatures. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(t(xCell.matrix)), "."),
                    cex = 2.2), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.55)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean xCell correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation.xCell[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  ### Hallmark
  
  # Correlation matrix Hallmark_GSEA
  print(paste0("Generating Hallmark_GSEA z score correlation matrix for ", Cancer))
  
  Hallmark.enrichment.z.score.df = as.data.frame(t(Hallmark.enrichment.z.score))
  Hallmark.enrichment.z.score.df$ICR_score = table_cluster_assignment$ICRscore[match(row.names(Hallmark.enrichment.z.score.df), row.names(table_cluster_assignment))]
  
  Hallmark_GSEA_cor <- cor (Hallmark.enrichment.z.score.df,method=test)
  
  # Correlation significance
  Hallmark_GSEA_cor_sign = cor.mtest(Hallmark.enrichment.z.score.df, 0.95)
  
  # Correlation plot
  print(paste0("Generating Hallmark_GSEA z score correlation plot for ", Cancer))
  png(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots/ssGSEA_correlation_plots/", download.method,
             "/", "Hallmark_GSEA_", test, "_Correlation_plot_", Cancer, ".png"), res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.35)
  lims=c(-1,1)
  if (length(Hallmark_GSEA_cor[Hallmark_GSEA_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(Hallmark_GSEA_cor),color = c(rep("#CC0506",nrow(Hallmark_GSEA_cor))),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(Hallmark_GSEA_cor,order="FPC"),]
  
  mean_correlation = round(mean(Hallmark_GSEA_cor),2)
  corrplot.mixed (Hallmark_GSEA_cor,
                  #type="lower",
                  #p.mat = Hallmark_GSEA_cor_sign[[1]],                                                                            # add significance to correlations
                  lower.col = "black",
                  upper.col = colpattern,
                  lower = "number",
                  upper = "square",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.8,
                  tl.cex = 0.5,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(8,4.1,6,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between Hallmark GSEA signatures. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(t(Hallmark.enrichment.z.score)), "."),
                    cex = 2.2), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.55)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean Hallmark GSEA correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation.Hallmark.GSEA[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  ### Hallmark signatures
  
  ## Correlation matrix Hallmark_CM (column Means)
  print(paste0("Generating Hallmark column means correlation matrix for ", Cancer))
  
  Hallmark.col.means.df = as.data.frame(t(Hallmark.col.means))
  Hallmark.col.means.df$ICR_score = table_cluster_assignment$ICRscore[match(row.names(Hallmark.col.means.df), row.names(table_cluster_assignment))]
  
  Hallmark_CM_cor <- cor (Hallmark.col.means.df, method=test)
  
  # Correlation significance
  Hallmark_CM_cor_sign = cor.mtest(Hallmark.col.means.df, 0.95)
  
  # Correlation plot
  print(paste0("Generating Hallmark column means correlation plot for ", Cancer))
  png(paste0("./5_Figures/Correlation_plots/Hallmark_Correlation_plots/colMeans_correlation_plots/", download.method,
             "/", "Hallmark_CM_", test, "_Correlation_plot_", Cancer, ".png"), res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.35)
  lims=c(-1,1)
  if (length(Hallmark_CM_cor[Hallmark_CM_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(Hallmark_CM_cor),color = c(rep("#CC0506",nrow(Hallmark_CM_cor))),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(Hallmark_CM_cor,order="FPC"),]
  
  mean_correlation = round(mean(Hallmark_CM_cor),2)
  corrplot.mixed (Hallmark_CM_cor,
                  #type="lower",
                  #p.mat = Hallmark_CM_cor_sign[[1]],                                                                            # add significance to correlations
                  lower.col = "black",
                  upper.col = colpattern,
                  lower = "number",
                  upper = "square",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.8,
                  tl.cex = 0.5,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(8,4.1,6,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between Hallmark CM signatures. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(t(Hallmark.enrichment.z.score)), "."),
                    cex = 2.2), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.55)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean Hallmark CM correlation is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation.Hallmark.CM[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  ## Combined correlation
  
  matrix_all_GSEA = rbind(Bindea.enrichment.z.score, xCell.matrix, Hallmark.enrichment.z.score)
  
  all_GSEA.df = as.data.frame(t(matrix_all_GSEA))
  all_GSEA.df$ICR_score = table_cluster_assignment$ICRscore[match(row.names(all_GSEA.df), row.names(table_cluster_assignment))]
  
  ## Correlation matrix Hallmark_CM (column Means)
  all_GSEA_cor <- cor (all_GSEA.df, method=test)
  
  # Correlation significance
  all_GSEA_cor_sign = cor.mtest(all_GSEA.df, 0.95)
  
  # Correlation plot
  png(paste0("./5_Figures/Correlation_plots/Combined_Signature_Correlation_plots/", download.method,
             "/", "Combined_signatures_", test, "_Correlation_plot_", Cancer, ".png"), res=600,height=6,width=6,unit="in")
  
  cex.before <- par("cex")
  par(cex = 0.35)
  lims=c(-1,1)
  if (length(all_GSEA_cor[all_GSEA_cor<0]) == 0) {lims=c(0,1)}
  annotation = data.frame (gene = rownames(all_GSEA_cor),color = c(rep("#CC0506", nrow(all_GSEA_cor))),stringsAsFactors = FALSE)
  annotation = annotation[corrMatOrder(all_GSEA_cor,order="FPC"),]
  
  mean_correlation = round(mean(all_GSEA_cor),2)
  corrplot.mixed (all_GSEA_cor,
                  #type="lower",
                  #p.mat = Hallmark_CM_cor_sign[[1]],                                                                            # add significance to correlations
                  lower.col = "black",
                  upper.col = colpattern,
                  lower = "number",
                  upper = "square",
                  order="FPC",
                  cl.lim=lims,                                                                                               # only positive correlations
                  tl.pos ="lt",
                  tl.col = as.character(annotation$color),
                  insig= "pch",                                                                                              # remove insignificant correlations
                  pch = "x",
                  pch.cex= 1.8,
                  tl.cex = 0.5,
                  cl.cex = 1/par("cex"),
                  cex.main = 1/par("cex"),
                  mar=c(8,4.1,6,5))
  title(main = list(paste0(Cancer, "\n", str_to_title(test), " correlation between bindea, xCell and Hallmark signatures CM signatures. \n ","Mean: ", mean_correlation,". Number of patients: ", nrow(t(matrix_all_GSEA)), "."),
                    cex = 2.2), line = -2.5, adj = 0.55)
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data was \n obtained from TCGA, using ", download.method, " v2.0.3. \n",
                          "Significance level of correlation is represented by the size of the squares."), cex = 1.5), line = 1.5, adj = 0.55)
  par(cex = cex.before)
  dev.off()
  
  cat(paste0("For ", Cancer, " mean correlation between all signatures is ", mean_correlation), file = Log_file, append = TRUE, sep = "\n")
  mean_correlation_table$Mean.correlation.bindea.xCell.Hallmark[mean_correlation_table$Cancertype == Cancer] = mean_correlation
  
  dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation"), showWarnings = FALSE)
  save(Bindea_cor, Bindea_cor_sign, xCell_cor, xCell_cor_sign, Hallmark_CM_cor, Hallmark_CM_cor_sign,
       Hallmark_GSEA_cor, Hallmark_GSEA_cor_sign, all_GSEA_cor, all_GSEA_cor_sign, file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/",
                                                                                                 "Correlation_matrixes_Bindea_xCell_Hallmark_", test, "_", Cancer, ".Rdata"))
}

dir.create(paste0("./4_Analysis/", download.method ,"/Pan_Cancer/Correlation"), showWarnings = FALSE)
save(mean_correlation_table, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Correlation/Correlation_Bindea_xCell_Hallmark.Rdata"))
end.time <- Sys.time ()
time <- end.time - start.time
print (paste0("Between start script and completion script: ", time))

