#################################################################
###
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

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                     # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                           # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/heatmap.3.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/5.1_Pancancer_Correlation_matrix_Signatures/5.1_Pancancer_Correlation_matrix_signatures", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"
display_correlations = "irrespective_of_significance"                                                                   # Can either be "only_significant" or "irrespective_of_significance"         

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)

if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/5.1_Pancancer_Correlation_matrix_Signatures"), showWarnings = FALSE)
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

start.time <- Sys.time ()

load(file = "./4_Analysis/TCGA_Assembler/ACC/Correlation/Correlation_matrixes_Bindea_xCell_Hallmark_pearson_ACC.Rdata") 

row.names(TCGA.cancersets) = TCGA.cancersets$cancerType
pancancer_Bindea_correlation_table = t(TCGA.cancersets)[-c(1,2),]
pancancer_all_correlation_table = pancancer_xCell_correlation_table = pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_CM_correlation_table = pancancer_Bindea_correlation_table

pancancer_Bindea_correlation_table = rbind(pancancer_Bindea_correlation_table,matrix(nrow = nrow(Bindea_cor),ncol=N.sets))
rownames(pancancer_Bindea_correlation_table) = rownames(Bindea_cor)

pancancer_xCell_correlation_table = rbind(pancancer_xCell_correlation_table,matrix(nrow = nrow(xCell_cor),ncol=N.sets))
rownames(pancancer_xCell_correlation_table) = rownames(xCell_cor)

pancancer_Hallmark_GSEA_correlation_table = rbind(pancancer_Hallmark_GSEA_correlation_table,matrix(nrow = nrow(Hallmark_GSEA_cor),ncol=N.sets))
rownames(pancancer_Hallmark_GSEA_correlation_table) = rownames(Hallmark_GSEA_cor)

pancancer_Hallmark_CM_correlation_table = rbind(pancancer_Hallmark_CM_correlation_table,matrix(nrow = nrow(Hallmark_CM_cor),ncol=N.sets))
rownames(pancancer_Hallmark_CM_correlation_table) = rownames(Hallmark_CM_cor)

pancancer_all_correlation_table = rbind(pancancer_all_correlation_table, matrix(nrow = nrow(all_GSEA_cor),ncol=N.sets))
rownames(pancancer_all_correlation_table) = rownames(all_GSEA_cor)

i = 1
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  
  if(display_correlations == "irrespective_of_significance"){
    load(paste0("./4_Analysis/",download.method, "/", Cancer, "/Correlation/", "Correlation_matrixes_Bindea_xCell_Hallmark_pearson_", Cancer, ".Rdata"))
  }
  if(display_correlations == "only_significant"){
    load(paste0("./4_Analysis/",download.method, "/", Cancer, "/Correlation/", "Correlation_matrixes_Bindea_xCell_Hallmark_pearson_", Cancer, ".Rdata"))
    
    N.columns = ncol(Bindea_cor)
    N.rows = nrow(Bindea_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        Bindea_cor[i,j] = ifelse(Bindea_cor_sign[[1]][i,j] <0.05 | Bindea_cor_sign[[1]][i,j] >0.95, Bindea_cor[i,j], 0)
      }
    }
    
    N.columns = ncol(xCell_cor)
    N.rows = nrow(xCell_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        xCell_cor[i,j] = ifelse(xCell_cor_sign[[1]][i,j] <0.05 | xCell_cor_sign[[1]][i,j] >0.95, xCell_cor[i,j], 0)
      }
    }
    N.columns = ncol(Hallmark_GSEA_cor)
    N.rows = nrow(Hallmark_GSEA_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        Hallmark_GSEA_cor[i,j] = ifelse(Hallmark_GSEA_cor_sign[[1]][i,j] <0.05 | Hallmark_GSEA_cor_sign[[1]][i,j] >0.95, Hallmark_GSEA_cor[i,j], 0)
      }
    }
    N.columns = ncol(Hallmark_CM_cor)
    N.rows = nrow(Hallmark_CM_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        Hallmark_CM_cor[i,j] = ifelse(Hallmark_CM_cor_sign[[1]][i,j] <0.05 | Hallmark_CM_cor_sign[[1]][i,j] >0.95, Hallmark_CM_cor[i,j], 0)
      }
    }
    N.columns = ncol(all_GSEA_cor)
    N.rows = nrow(all_GSEA_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        all_GSEA_cor[i,j] = ifelse(all_GSEA_cor_sign[[1]][i,j] <0.05 | all_GSEA_cor_sign[[1]][i,j] >0.95, all_GSEA_cor[i,j], 0)
      }
    }
  }
  
  pancancer_Bindea_correlation_table[, Cancer] = as.numeric(Bindea_cor[,"ICR_score"])
  pancancer_xCell_correlation_table[, Cancer] = as.numeric(xCell_cor[,"ICR_score"])
  pancancer_Hallmark_GSEA_correlation_table[, Cancer] = as.numeric(Hallmark_GSEA_cor[,"ICR_score"])
  pancancer_Hallmark_CM_correlation_table[, Cancer] = as.numeric(Hallmark_CM_cor[,"ICR_score"])
  pancancer_all_correlation_table[, Cancer] = as.numeric(all_GSEA_cor[,"ICR_score"])
}

# convert to numeric matrix
mode(pancancer_Bindea_correlation_table) = "numeric"
mode(pancancer_xCell_correlation_table) = "numeric"
mode(pancancer_Hallmark_GSEA_correlation_table) = "numeric"
mode(pancancer_Hallmark_CM_correlation_table) = "numeric"
mode(pancancer_all_correlation_table) = "numeric"

# Remove LAML from all correlation tables
pancancer_Bindea_correlation_table = pancancer_Bindea_correlation_table[,-c(14)]
pancancer_xCell_correlation_table = pancancer_xCell_correlation_table[,-c(14)]
pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[,-c(14)]
pancancer_Hallmark_CM_correlation_table = pancancer_Hallmark_CM_correlation_table[,-c(14)]
pancancer_all_correlation_table = pancancer_all_correlation_table[,-c(14)]

## Bindea correlation heatmap
#png(paste0("./5_Figures/Pancancer_plots/Bindea_ICR_correlation_", display_correlations, ".png"),res=600,height= 10,width=10,unit="in")
heatmap.3 (pancancer_Bindea_correlation_table,
          col= my.palette,
          main = paste0("Pan-Cancer ", test, " correlation \n between ICR and Bindea gene signatures "),
          cex.main = 0.7,
          cexCol = 1.2,
          cexRow = 1.5,
          margins=c(14, 12))

title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -1)
#dev.off()

## xCell correlation heatmap
#png(paste0("./5_Figures/Pancancer_plots/xCell_ICR_correlation_", display_correlations, ".png"),res=600,height= 15,width=13,unit="in")
heatmap.3 (pancancer_xCell_correlation_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation \n between ICR and xCell gene signatures "),
           cex.main = 0.7,
           cexCol = 1.2,
           cexRow = 1.2,
           margins=c(12, 15))

title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -0.8)
#dev.off()

## Hallmark_GSEA correlation heatmap
#png(paste0("./5_Figures/Pancancer_plots/Hallmark_GSEA_ICR_correlation_", display_correlations, ".png"),res=600,height= 16, width=16, unit="in")
heatmap.3 (pancancer_Hallmark_GSEA_correlation_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation between \nICR and Hallmark pathways (GSEA) "),
           cex.main = 0.7,
           cexCol = 1.2,
           cexRow = 1.2,
           margins=c(12, 29))

title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -1)
#dev.off()

## Hallmark_CM correlation heatmap
#png(paste0("./5_Figures/Pancancer_plots/Hallmark_CM_ICR_correlation_", display_correlations, ".png"),res=600,height= 16, width=16, unit="in")
heatmap.3 (pancancer_Hallmark_CM_correlation_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation between \nICR and Hallmark pathways (CMs) "),
           cex.main = 0.7,
           cexCol = 1.2,
           cexRow = 1.2,
           margins=c(12, 29))

title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -1)
#dev.off()

## Hallmark_all correlation heatmap
#png(paste0("./5_Figures/Pancancer_plots/All_combined_ICR_correlation_", display_correlations, ".png"),res=600,height= 30, width=20, unit="in")
heatmap.3 (pancancer_all_correlation_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation between \nICR and other gene signatures "),
           cex.main = 0.7,
           cexCol = 1.2,
           cexRow = 1.2,
           margins=c(12, 29))

title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -1)
dev.off()

dir.create("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation", showWarnings = FALSE)
save(pancancer_Bindea_correlation_table, pancancer_xCell_correlation_table, pancancer_Hallmark_GSEA_correlation_table,
     pancancer_Hallmark_CM_correlation_table, pancancer_all_correlation_table,
     file = paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/", "Correlation_Bindea_xCell_Hallmark_", display_correlations,".Rdata"))
end.time <- Sys.time ()
time <- end.time - start.time
print (paste0("Between start script and completion script: ", time))
