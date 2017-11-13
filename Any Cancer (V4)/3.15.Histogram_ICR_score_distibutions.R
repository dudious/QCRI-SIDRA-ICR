#################################################################
###
### This script makes histograms of the ICR scores to visualize 
### the distribution of ICR scores across samples in each
### individual cancer.
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c()
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/3.15_Histogram_ICR/3.15_Histogram_ICR",                                                # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
Score_type = "ICRscore"                                                                                                 # Score_type can be "ICRscore" or "Scaled_ICRscore"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create folders and log file
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/3.15_Histogram_ICR"), showWarnings = FALSE)
cat("This is a log file for creating Boxplots with mean ICR score across ICRs",                                         # Set-up logfile
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
    "Creating boxplots",
    file = Log_file,
    append = FALSE, sep= "\n")

start.time.process.all = Sys.time()
msg = paste0("Processing", "\n")
cat(msg)

i= 1
for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(!file.exists(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
                         assay.platform, "_", "normalized.Rdata"))) 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  
  cat(paste0("Making histogram for ", Cancer), file = Log_file, sep = "\n", append = TRUE)
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  dir.create("./5_Figures", showWarnings = FALSE)
  dir.create("./5_Figures/ICR_distribution_plots", showWarnings = FALSE)
  dir.create(paste0("./5_Figures/ICR_distribution_plots/ICR_score_histograms/"), showWarnings = FALSE)
  dir.create(paste0("./5_Figures/ICR_distribution_plots/ICR_score_histograms/", download.method), showWarnings = FALSE)
  
  png(paste0("./5_Figures/ICR_distribution_plots/ICR_score_histograms/", download.method, "/Histogram_", Score_type , "_",Cancer,".png"), res=600,height=6,width=15,unit="in")
  #histogram_ICR2 =  ggplot(data = table_cluster_assignment,
                    #aes(x = table_cluster_assignment$Scaled_ICRscore)) + geom_histogram() + facet_wrap(~row.names(table_cluster_assignment), ncol = 1, scales= "free_x")
  #print(histogram_ICR2)
  
  histogram_ICR = qplot(table_cluster_assignment$ICRscore,
                        bins = 60,
                        main = paste0("Histogram ", Score_type," in ", Cancer),
                        xlab = Score_type,
                        fill = I("#EB6D7D"),
                        col = I("black")) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=16,face="bold"),
          title = element_text(size=17))
  print(histogram_ICR)
  
  dev.off()
}

