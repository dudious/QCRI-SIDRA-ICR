#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Correlation", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                     # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                           # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages <- c("stringr", "ggplot2")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "TCGA_Assembler"                                                                                       # Specify download method (this information to be used when saving the file)
                                                                                        
Log_file = paste0("./1_Log_Files/5.2_Pancancer_ICR_Correlation_boxplot/5.2_Pancancer_ICR_Correlation_boxplot", 
                  "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
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
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/5.2_Pancancer_ICR_Correlation_boxplot"), showWarnings = FALSE)
cat("This is a log file for creating a pancancer correlation plot",
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

ICR_correlation_table = data.frame(Correlation_values = 0, Cancertype = NA)

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
  
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Correlation/Correlation_matrix_ICR_", test, "_", Cancer, ".Rdata"))
  
  correlation_triangle = ICR_cor[upper.tri(ICR_cor,diag = FALSE)]
  Cancertype = rep(Cancer, length(correlation_triangle))
  ICR_correlation_table = rbind(ICR_correlation_table, data.frame(Correlation_values = correlation_triangle, Cancertype = Cancertype))
}
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))
ICR_cluster_assignment_ICR_High = ICR_cluster_assignment_allcancers[ICR_cluster_assignment_allcancers$HML_cluster == "ICR High",]
mean_ICR_perCancer = aggregate(ICR_cluster_assignment_ICR_High$ICRscore, list(ICR_cluster_assignment_ICR_High$Cancer), mean)
mean_ICR_perCancer = mean_ICR_perCancer[order(mean_ICR_perCancer$x, decreasing = FALSE),]
Cancer_order = mean_ICR_perCancer$Group.1

ICR_correlation_table = ICR_correlation_table[-c(1),]

ICR_correlation_table = ICR_correlation_table[order(match(ICR_correlation_table$Cancertype,Cancer_order)),]
ICR_correlation_table$Cancertype = factor(ICR_correlation_table$Cancertype,levels = Cancer_order)

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

Cancer_color_table = Cancer_color_table[order(match(Cancer_color_table$Group.1, Cancer_order)),]
     
png(paste0("./5_Figures/Correlation_plots/ICR_Correlation_plots/",
           download.method, "/ICR_", test, "_Correlation_boxplot_pancancer.png"), res=600,height=6,width=14,unit="in")

boxplot_correlation = ggplot(data = ICR_correlation_table, aes(x= Cancertype, y= Correlation_values, fill = Cancertype)) + 
  geom_boxplot(outlier.colour = NA) +
  labs(fill = "Cancer type") +
  scale_y_continuous("Mean correlation") +
  scale_fill_manual(values = Cancer_color_table$color) + 
  ggtitle("Pearson correlation between expression of ICR genes across cancers") + 
  theme(plot.title = element_text(size=14, face = "bold")) +
  guides(fill=FALSE) + 
  theme(axis.title = element_text(size = 13)) + 
  theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "grey", size = 0.2), panel.grid.minor.y = element_line(colour = "grey", size = 0.2)) +
  theme(axis.line = element_line(color= "black", size = 0.4))
print(boxplot_correlation)
dev.off()

