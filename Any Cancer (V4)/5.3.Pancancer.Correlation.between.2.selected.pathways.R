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

required.packages <- c("ggplot2")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                 # Specify which genes will be correlated
assay.platform = "gene_RNAseq"
pathway1 = "MAPK UP GENES"
pathway2 = "[IPA] ERK MAPK Signaling"
  
# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)

if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Signature_Analysis/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Signature_Analysis/Correlation_plots"), showWarnings = FALSE)

N.sets = length(CancerTYPES)

load(file = "./4_Analysis/TCGA_Assembler/ACC/Correlation/Correlation_matrixes_Bindea_xCell_Hallmark_pearson_ACC.Rdata") 

pathway_correlation_df = data.frame(Correlation_values = 0, Cancertype = NA, Group = NA)

i=1
for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), sep = "\n", append = TRUE)
    next}
  
  load(file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/Correlation_matrixes_Bindea_xCell_Hallmark_pearson_", Cancer, ".Rdata"))
  
  pathway_correlation_df = rbind(pathway_correlation_df, data.frame(Correlation_values = Hallmark_GSEA_cor[pathway1, pathway2], Cancertype = Cancer, Group = "Tumor"))
}

pathway_correlation_df = pathway_correlation_df[-1,]
  
png(paste0("./5_Figures/Signature_Analysis/Correlation_plots/", pathway1, "_vs_", pathway2, "correlation_boxplot" ,".png"), res=600,height=7,width=4,unit="in")

boxplot_correlation = ggplot(data = pathway_correlation_df, aes(x= Group, y= Correlation_values)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_jitter() +
  ylim(-0.5,1) +
  #labs(fill = "Cancer type") +
  #scale_y_continuous("Mean correlation") +
  scale_fill_manual(my.palette) + 
  ggtitle(paste0("Mean correlation between \n", pathway1, " and \n", pathway2, " \nin each individual cancer")) + 
  theme(plot.title = element_text(size=11, face = "bold")) +
  guides(fill=FALSE) + 
  theme(axis.title = element_text(size = 13)) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2)) +
  theme(axis.line = element_line(color= "black", size = 0.4))
print(boxplot_correlation)
dev.off()

