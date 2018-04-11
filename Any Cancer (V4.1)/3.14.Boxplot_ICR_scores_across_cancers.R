#################################################################
###
### This script makes boxplots of ICR scores in all cancers.
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
download.method = "Assembler_Panca_Normalized"                                                                                      # Specify download method "TCGA_Assembler", "Pancancer_matrix" (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/", download.method, "/3.14_Boxplot_ICR/3.14_Boxplot_ICR",                                                    # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
ICR_classification_k = "HML_classification"
basis_ordering = "Mean ICR High"
subset = "all"                                                                                                          # Subset can be "ICR High", "ICR Low", "All"
my.palette = colorRampPalette(c("#BEB9DA", "#FFD82F", "#BC7FBC", "#666666", "#387EB7", 
                                "#FEB462", "#E72A89", "#F781BF", "#B3B3B3", "#396BAF", 
                                "#FDCDAC", "#FFFD99", "#FC9998", "#B2DE68", "#A6CEE3",
                                "#F12B7E", "#4DAE4B", "#D96002", "#FB7F72", "#E4211E",
                                "#FC8D62", "#DECAE4", "#F1E1CC", "#D9D9D9", "#E6AB03",
                                "#FFFECC", "#7FB1D3", "#FFED6F", "#E5C494", "#A5761D",
                                "#34A02C", "#8DD3C7"))

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
load("./4_Analysis/Pancancer_matrix/Pan_Cancer/Clustering/Cancer_color_table.Rdata")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

# Create folders and log file
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method ,"/3.14_Boxplot_ICR"), showWarnings = FALSE)
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

for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "LAML") 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  
  cat(paste0("Making boxplot for ", Cancer), file = Log_file, sep = "\n", append = TRUE)
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  table_cluster_assignment$Cancer = Cancer
  assign(paste0(Cancer, "_table_cluster_assignment"), table_cluster_assignment)
}

ICR_cluster_assignment_allcancers = ACC_table_cluster_assignment

for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "LAML") 
  {next}
  if (Cancer == "ACC") {next}
  ICR_cluster_assignment_allcancers = rbind(get(paste0(Cancer, "_table_cluster_assignment")), ICR_cluster_assignment_allcancers )
}

#load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))
#ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[ICR_cluster_assignment_allcancers$HML_cluster == subset,]

dir.create("./5_Figures", showWarnings = FALSE)
dir.create("./5_Figures/ICR_distribution_plots", showWarnings = FALSE)
dir.create(paste0("./5_Figures/ICR_distribution_plots/ICR_boxplot_across_cancers/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/ICR_distribution_plots/ICR_boxplot_across_cancers/", download.method), showWarnings = FALSE)

ICR_cluster_assignment_allcancers$HML_cluster = as.factor(ICR_cluster_assignment_allcancers$HML_cluster)
ICR_cluster_assignment_allcancers$HML_cluster = relevel(ICR_cluster_assignment_allcancers$HML_cluster, "ICR High")
ICR_cluster_assignment_allcancers$HML_cluster = relevel(ICR_cluster_assignment_allcancers$HML_cluster, "ICR Medium")
ICR_cluster_assignment_allcancers$HML_cluster = relevel(ICR_cluster_assignment_allcancers$HML_cluster, "ICR Low")

if(basis_ordering == "Mean ICR"){
  mean_ICR_perCancer = aggregate(ICR_cluster_assignment_allcancers$ICRscore, list(ICR_cluster_assignment_allcancers$Cancer), mean)
  mean_ICR_perCancer = mean_ICR_perCancer[order(mean_ICR_perCancer$x, decreasing = FALSE),]
  Cancer_order = mean_ICR_perCancer$Group.1
}

if(basis_ordering == "Mean ICR High"){
  ICR_cluster_assignment_ICR_High = ICR_cluster_assignment_allcancers[ICR_cluster_assignment_allcancers$HML_cluster == "ICR High",]
  mean_ICR_perCancer = aggregate(ICR_cluster_assignment_ICR_High$ICRscore, list(ICR_cluster_assignment_ICR_High$Cancer), mean)
  mean_ICR_perCancer = mean_ICR_perCancer[order(mean_ICR_perCancer$x, decreasing = FALSE),]
  Cancer_order = mean_ICR_perCancer$Group.1
}

#ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[order(match(ICR_cluster_assignment_allcancers$Cancer, Cancer_order)),]
ICR_cluster_assignment_allcancers$Cancer = factor(ICR_cluster_assignment_allcancers$Cancer,levels = Cancer_order)

Cancer_color_table = Cancer_color_table[order(match(Cancer_color_table$Group.1, Cancer_order)),]

png(paste0("./5_Figures/ICR_distribution_plots/ICR_boxplot_across_cancers/", 
           download.method, "/ICR_distribution_boxplot_ordered_by_", basis_ordering, "_", ICR_classification_k,
           "_colour_by_cancer_across_cancers.png"), res=600,height=6,width=15,unit="in")
boxplot_ICR = ggplot(data = ICR_cluster_assignment_allcancers, aes(x= HML_cluster, y= ICRscore, fill = Cancer)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = Cancer_color_table$color) +
  scale_y_continuous("ICR score") +
  facet_grid(. ~ Cancer) +
  ggtitle(paste0("ICR scores in ", subset, " according to ", ICR_classification_k, " across cancers ordered by ", basis_ordering)) +
  theme(plot.title = element_text(size=14, face = "bold")) + theme(axis.title = element_text(size = 15, face = "bold", colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_text(angle = 90, colour = "black"),
        strip.text.x = element_text(colour = "black")) +
  theme(axis.line = element_line(color= "black", size = 0.4)) +
  guides(fill=FALSE)
print(boxplot_ICR)
dev.off()

png(paste0("./5_Figures/ICR_distribution_plots/ICR_boxplot_across_cancers/", 
           download.method, "/ICR_distribution_boxplot_ordered_by_", basis_ordering, "_", ICR_classification_k,
           "_colour_by_ICR_cluster_across_cancers.png"), res=600,height=6,width=15,unit="in")
boxplot_ICR = ggplot(data = ICR_cluster_assignment_allcancers, aes(x= HML_cluster, y= ICRscore, fill = HML_cluster)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = c("blue", "green", "red")) +
  scale_y_continuous("ICR score") +
  facet_grid(. ~ Cancer, switch = "x") +
  ggtitle(paste0("ICR scores in ", subset, " according to ", ICR_classification_k, " across cancers ordered by ", basis_ordering)) +
  xlab("Cancer types") +
  theme(plot.title = element_text(size=14, face = "bold")) + theme(axis.title = element_text(size = 15, face = "bold", colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_blank(),
        strip.text.x = element_text(colour = "black", size = 10),
        strip.placement = "outside",
        strip.background = element_blank()) +
  theme(axis.line = element_line(color= "black")) +
  guides(fill=FALSE)
print(boxplot_ICR)
dev.off()

save(ICR_cluster_assignment_allcancers, Cancer_color_table, file =paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))


## Extra plot

png(paste0("./5_Figures/ICR_distribution_plots/ICR_boxplot_across_cancers/", 
           download.method, "/ICR_distribution_boxplot_ordered_by_", basis_ordering, "_Mean_ICR",
           "_colour_by_cancer_across_cancers.png"), res=600,height=6,width=15,unit="in")
boxplot_ICR = ggplot(data = ICR_cluster_assignment_allcancers, aes(x= Cancer, y= ICRscore, fill = Cancer)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = Cancer_color_table$color) +
  scale_y_continuous("ICR score") +
  ggtitle(paste0("Mean ICR scores in ", subset, " across cancers ordered by ", basis_ordering)) +
  theme(plot.title = element_text(size=14, face = "bold")) + theme(axis.title = element_text(size = 15, face = "bold", colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_text(angle = 90, colour = "black")) +
  theme(axis.line = element_line(color= "black", size = 0.4)) +
  guides(fill=FALSE)
print(boxplot_ICR)
dev.off()
