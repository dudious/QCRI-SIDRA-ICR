#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
### Manual adjustment of min and max
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize")
ibiopak("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,"/5.5.1.Pancancer_Bindea_Heatmap/5.5.1.Pancancer.Bindea.Heatmap_", # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("ICR clusters", "Cancers", "Proliferation cluster")
Legend = c("ICR Low","ICR Med","ICR High")
Cutoff_HR = 1
Proliferation_classification = "P_Pancancer"

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
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/5.5.1.Pancancer_Bindea_Heatmap"), showWarnings = FALSE)

cat("This is a log file for creating heatmap bindea enrichment plots",
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

load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", "Survival_analysis_High_vs_Low_Groups", 
            "HML_classification", ".Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]

load("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata")
load(paste0("./4_Analysis/", download.method, "/ACC/Signature_Enrichment/GSEA_ACC_Proliferation_metagene.Rdata"))

Enrichment_score_all = Proliferation.metagene.enrichment.score

for (i in 2:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer == "LAML"){next}
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, "_Proliferation_metagene.Rdata"))
  Enrichment_score_all = cbind(Enrichment_score_all, Proliferation.metagene.enrichment.score)
}

Enrichment_score_all.z.score = Enrichment_score_all 
for(j in 1: nrow(Enrichment_score_all.z.score))  {
  Enrichment_score_all.z.score[j,] = (Enrichment_score_all[j,]-mean(Enrichment_score_all[j,]))/sd(Enrichment_score_all[j,]) # z-score the enrichment matrix
}

if(Proliferation_classification == "P_Per_Cancer"){
  annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer", "Proliferation_cluster_name")]
  
}

if(Proliferation_classification == "P_Pancancer"){
  ## Bindea cluster (hierarchical)
  sHc = hclust(ddist <- dist(t(Enrichment_score_all.z.score)), method = "ward.D2")
  
  plot(sHc,labels=FALSE)
  
  ## Annotation for plotting
  annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer") ,drop=FALSE]
  # Pancancer Proliferation classification
  annotation$proliferation_cluster = cutree(sHc,k = 2)[match(rownames(annotation),names(cutree(sHc,k = 2)))]
  annotation = annotation[colnames(Enrichment_score_all.z.score),]
  annotation$proliferation_score = Enrichment_score_all.z.score[1,]
  proliferation_cluster_means = aggregate(proliferation_score~proliferation_cluster,data=annotation,FUN=mean)
  proliferation_cluster_means = proliferation_cluster_means[order(proliferation_cluster_means$proliferation_score),]
  proliferation_cluster_means$proliferation_cluster_name = c("Proliferation Low","Proliferation High")
  annotation$proliferation_cluster_name = proliferation_cluster_means$proliferation_cluster_name[match(annotation$proliferation_cluster,
                                                                                                       proliferation_cluster_means$proliferation_cluster)]
}

Cancer_order = CancerTYPES[-which(CancerTYPES =="LAML")]
ICR_order = c("ICR Low","ICR Medium","ICR High")
Proliferation_order = c("Proliferation Low", "Proliferation High")

#annotation = annotation[order(match(annotation$Cancer, Cancer_order)),]
#annotation = annotation[order(match(annotation$HML_cluster,ICR_order)),]

annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"

annotation$Cancer.col = Cancer_color_table$color[match(annotation$Cancer, Cancer_color_table$Group.1)]

annotation$Proliferation.col[annotation$proliferation_cluster_name == "Proliferation High"] = "#B97474"
annotation$Proliferation.col[annotation$proliferation_cluster_name == "Proliferation Low"] = "#ADDFAD"

annotation.blot = as.matrix(annotation[,c("HML_cluster.col","Cancer.col", "Proliferation.col"), drop = FALSE])
#annotation.blot = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
#Expression.data = Expression.data[colnames(annotation.blot),]

#Enrichment_score_all.z.score = Enrichment_score_all.z.score[,rownames(annotation.blot), drop = FALSE]
Enrichment_score_all.z.score = rbind(Enrichment_score_all.z.score, rep(0, 9475))

## Plot prep
Legend = CancerTYPES[-which(CancerTYPES == "LAML")]
rownames(Cancer_color_table) = Cancer_color_table$Group.1
Cancer_color_table = Cancer_color_table[CancerTYPES[-which(CancerTYPES == "LAML")],]
Legend_colors = c(Cancer_color_table$color)

Legend2 = c("ICR Low","ICR Med","ICR High")
Legend_colors2 = c("blue","green","red")

Legend3 = c("Proliferation High", "Proliferation Low")
Legend_colors3 = c("#B97474", "#ADDFAD")

### Plotting
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Proliferation_Clusters_determined_PanCancer", "_Heatmap_RNASeq_Pancancer_HML.png"),
    res = 600, height = 8, width = 15, unit = "in")
heatmap.3((as.matrix(Enrichment_score_all.z.score)),
          main= paste0("Pancancer enrichment scores \nssGSEA/Proliferation metagene Lance Miller"),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          Colv= as.dendrogram(sHc),
          dendrogram = "column",
          #Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.3,
          cexRow = 1.3,
          margins = c(13, 30))

title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                        "Oncogenic pathway enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
legend("bottomleft",legend = Legend,
       col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("top", legend = Legend2,
       col = Legend_colors2, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("topright", legend = Legend3,
       col = Legend_colors3, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)

dev.off()

if(Proliferation_classification == "P_Pancancer"){
  ICR_cluster_assignment_allcancers$Proliferation_cluster_Pancancer = annotation$proliferation_cluster_name[match(rownames(ICR_cluster_assignment_allcancers),
                                                                                                                  rownames(annotation.blot))]
  save(Cancer_color_table, ICR_cluster_assignment_allcancers, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))
}

