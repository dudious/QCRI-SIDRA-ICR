#################################################################
###
###
### Data input:
###
### Output :
### 
#################################################################


## Remark: The forest plot function has issues with zero and/or inf number values.
## Beware of this for troubleshooting!

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/heatmap.3.R"))

required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

## Set Parameters
CancerTYPES = c("ALL")                                                                                                   # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Assembler_Panca_Normalized"                                                                                       # Specify download method (this information to be used when saving the file)
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/", download.method, "/5.1_Pancancer_Correlation_matrix_Signatures/5.1_Pancancer_Correlation_matrix_signatures", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations)
load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Survival_Analysis/Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

# Add number of patients per cancer to All_survival_analysis data
for (i in 1:length(TCGA.cancersets$cancerType)){
  Cancer = TCGA.cancersets$cancerType[i]
  All_survival_analysis_data$N[All_survival_analysis_data$Cancertype == Cancer] = sum(ICR_cluster_assignment_allcancers$Cancer == Cancer)
}

if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Remove LAML and TGCT(has zero and inf values)
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "LAML"),]
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "THYM"),]
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "TGCT"),]
All_survival_analysis_data$p_value = round(All_survival_analysis_data$p_value, digits = 3)
HR.table = All_survival_analysis_data

#HR.table <- HR.table[-which(is.infinite(HR.table$Upper)),]
#HR.table <- HR.table[-which(HR.table$Lower==0),]

n.cancers = length(CancerTYPES)
x = n.cancers + 2

HR.table = HR.table[order(HR.table$HR),]

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR.table$HR[1:n.cancers]), NA),
    lower = c(NA,HR.table$CI_lower[c(1:n.cancers)], NA),
    upper = c(NA,HR.table$CI_upper[c(1:n.cancers)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")

tabletext<-cbind(
  c("Cancer", as.character(HR.table$Cancertype)[c(1:n.cancers)]),
  c("N", HR.table$N[c(1:n.cancers)]),
  c("p-value", HR.table$p_value[c(1:n.cancers)]),
  c("HR",      HR.table$HR[c(1:n.cancers)]))

dev.new()
forestplot(tabletext,
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,n.cancers),TRUE,rep(FALSE,n.cancers),TRUE,FALSE),
           #clip=c(0.1,2.5),
           xlog=TRUE,
           boxsize = .25,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 10), cex = 1))

