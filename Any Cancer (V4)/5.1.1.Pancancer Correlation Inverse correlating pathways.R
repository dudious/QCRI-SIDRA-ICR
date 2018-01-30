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

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/heatmap.3.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = c("ALL")                                                                                                   # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("PRAD",
                "DLBC", "CHOL")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/5.1_Pancancer_Correlation_matrix_Signatures/5.1_Pancancer_Correlation_matrix_signatures", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"
display_correlations = "irrespective_of_significance"                                                                   # Can either be "only_significant" or "irrespective_of_significance"
IPA_excluded = "IPA_excluded"
MAPK_down_excluded = "MAPK_down_excluded"
pcor_pathways_excluded = "pcor_pathways_excluded"
cor_cutoff = -0.1
positive_to_zero = "positive_to_zero"
minimum_number_of_cancers = 2

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
dir.create(paste0("./5_Figures/Pancancer_plots/Cancer_clustering_ICR_vs_inverse_HallmarkPathways"), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/5.1.1_Pancancer_Correlation_matrix_Signatures"), showWarnings = FALSE)
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

load(file = "./4_Analysis/TCGA_Assembler/READ/Correlation/Correlation_matrixes_Bindea_xCell_Hallmark_pearson_READ.Rdata") 

TCGA.cancersets = TCGA.cancersets[which(TCGA.cancersets$cancerType %in% CancerTYPES),]
row.names(TCGA.cancersets) = TCGA.cancersets$cancerType

pancancer_Hallmark_GSEA_correlation_table = t(TCGA.cancersets)[-c(1,2),]

pancancer_Hallmark_GSEA_correlation_table = rbind(pancancer_Hallmark_GSEA_correlation_table,matrix(nrow = nrow(Hallmark_GSEA_cor),ncol=N.sets))
rownames(pancancer_Hallmark_GSEA_correlation_table) = rownames(Hallmark_GSEA_cor)

if(length(Cancer_skip) >= 1) {
  pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[, -which(colnames(pancancer_Hallmark_GSEA_correlation_table)
                                                                                                 %in% Cancer_skip)]
}

i=2
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
    
    N.columns = ncol(Hallmark_GSEA_cor)
    N.rows = nrow(Hallmark_GSEA_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        Hallmark_GSEA_cor[i,j] = ifelse(Hallmark_GSEA_cor_sign[[1]][i,j] <0.05 | Hallmark_GSEA_cor_sign[[1]][i,j] >0.95, Hallmark_GSEA_cor[i,j], 0)
      }
    }}
  
  pancancer_Hallmark_GSEA_correlation_table[, Cancer] = as.numeric(Hallmark_GSEA_cor[,"ICR_score"])
}

# convert to numeric matrix
mode(pancancer_Hallmark_GSEA_correlation_table) = "numeric"

# Remove LAML from correlation table
pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[,-which(colnames(pancancer_Hallmark_GSEA_correlation_table) == "LAML")]

if(IPA_excluded == "IPA_excluded"){
  pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[grep(pattern = "IPA]", rownames(pancancer_Hallmark_GSEA_correlation_table), invert = TRUE),]
}

if(MAPK_down_excluded == "MAPK_down_excluded"){
  pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[grep(pattern = "MAPK DOWN GENES", rownames(pancancer_Hallmark_GSEA_correlation_table), invert = TRUE),]
}

if(pcor_pathways_excluded == "pcor_pathways_excluded"){
  test_pcor = pancancer_Hallmark_GSEA_correlation_table < cor_cutoff
  test_pcor = as.data.frame(test_pcor)
  test_pcor$rowsums = rowSums(test_pcor)
  pathways_to_exclude = rownames(test_pcor)[which(test_pcor$rowsums <= (minimum_number_of_cancers-1))]                                                                                                              
  pancancer_Hallmark_GSEA_correlation_table = pancancer_Hallmark_GSEA_correlation_table[-which(rownames(pancancer_Hallmark_GSEA_correlation_table) %in% pathways_to_exclude),]
}

if(positive_to_zero == "positive_to_zero"){
  pancancer_Hallmark_GSEA_correlation_table[pancancer_Hallmark_GSEA_correlation_table > 0] = 0
}


## Hallmark_GSEA correlation heatmap
png(paste0("./5_Figures/Pancancer_plots/Cancer_clustering_ICR_vs_inverse_HallmarkPathways/Hallmark_GSEA_ICR_correlation_", display_correlations, "_", IPA_excluded, "_",
           MAPK_down_excluded, "_",
           pcor_pathways_excluded, "_", cor_cutoff, "_", positive_to_zero, "_in_minimum_number_of_cancers_is_", minimum_number_of_cancers,"_excluded_cancers_are_", 
           paste(Cancer_skip, collapse ="_"), ".png"),res=600,height= 16, width=16, unit="in")
heatmap.3 (pancancer_Hallmark_GSEA_correlation_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation between \nICR and Hallmark pathways (GSEA) \n", positive_to_zero),
           cex.main = 0.7,
           cexCol = 1.2,
           cexRow = 1.2,
           margins=c(12, 29))

title(sub = list(paste0(gsub("_", " ", str_to_title(display_correlations)), ". Pathways that correlate < ", cor_cutoff, 
                        " with ICR for at least ", minimum_number_of_cancers  ," cancers are included."), cex = 1), outer = FALSE, line = -1)
dev.off()


dir.create("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation", showWarnings = FALSE)
save(pancancer_Hallmark_GSEA_correlation_table, 
     file = paste0("./4_Analysis/TCGA_Assembler/Pan_Cancer/Correlation/", "Correlation_Hallmark_", display_correlations, "_", IPA_excluded, "_",  pcor_pathways_excluded, "_",
                   "excluded_cancers_are_", paste(Cancer_skip, collapse = "_minimum_n_cancers_corr_<", cor_cutoff, "_" ,minimum_number_of_cancers), ".Rdata"))
end.time <- Sys.time ()
time <- end.time - start.time
print (paste0("Between start script and completion script: ", time))
