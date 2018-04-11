#################################################################
###
### This script compares ICR characteristics of ICR enabled vs. ICR
### disabled cancers.
###
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("ggpubr")
ipak(required.packages)

# Set parameters
download.method = "Assembler_Panca_Normalized"  
Cutoff_HR = 1

# Load data
Cluster_assignment_analysis = read.csv(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Clustering/Cluster_assignment_analysis.csv"), stringsAsFactors = FALSE)
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))
All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Cancer_order = as.character(All_survival_analysis_data$Cancertype[-which(All_survival_analysis_data$Cancertype == "LAML")])
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]

# ICR mean
mean_ICR_perCancer = aggregate(ICR_cluster_assignment_allcancers$ICRscore, list(ICR_cluster_assignment_allcancers$Cancer), mean)
#mean_ICR_perCancer = mean_ICR_perCancer[order(mean_ICR_perCancer$x, decreasing = FALSE),]
colnames(mean_ICR_perCancer) = c("Cancer", "Mean ICR score")
mean_ICR_perCancer$`ICR ID` = NA
mean_ICR_perCancer$`ICR ID`[which(mean_ICR_perCancer$Cancer %in% ICR_enabled_cancers)] = "ICR enabled"
mean_ICR_perCancer$`ICR ID`[which(mean_ICR_perCancer$Cancer %in% ICR_disabled_cancers)] = "ICR disabled"
mean_ICR_perCancer = mean_ICR_perCancer[order(match(mean_ICR_perCancer$Cancer, Cancer_order)),]
mean_ICR_perCancer$Cancer = factor(as.character(mean_ICR_perCancer$Cancer), levels = Cancer_order)

Cancer_color_table = Cancer_color_table[order(match(Cancer_color_table$Group.1, Cancer_order)),]

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Compare_ICRscore_ICRE_ICRD"), showWarnings = FALSE)
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Compare_ICRscore_ICRE_ICRD/Mean_ICR_score_boxplot_ICR_EnabledvsDisabled.png"), res=600,height=6,width=6,unit="in")
ggplot(mean_ICR_perCancer, aes(x= `ICR ID`, y= `Mean ICR score`)) + 
  stat_compare_means(method = "t.test") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(0.3), aes(colour = Cancer)) +
  scale_color_manual(values = Cancer_color_table$color) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.line = element_line(color= "black", size = 0.4),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text = element_text(size = 12, colour = "black"))
dev.off()

t.test(`Mean ICR score`~`ICR ID`, data = mean_ICR_perCancer)

# mean ICR score ICR High
Cluster_assignment_analysis$ICR_ID = NA
Cluster_assignment_analysis$ICR_ID[which(Cluster_assignment_analysis$Cancertype %in% ICR_enabled_cancers)] = "ICR enabled"
Cluster_assignment_analysis$ICR_ID[which(Cluster_assignment_analysis$Cancertype %in% ICR_disabled_cancers)] = "ICR disabled"

mean_ICR_High_perCancer = Cluster_assignment_analysis[-which(Cluster_assignment_analysis$Cancertype == "LAML"), c("Cancertype", "mean.ICR.scores.HML", "ICR_ID")]
colnames(mean_ICR_High_perCancer) = c("Cancer", "mean.ICR.scores.HML", "ICR ID")
mean_ICR_High_perCancer$`Mean ICR score ICR High` = gsub(".*/ " ,"", mean_ICR_High_perCancer$mean.ICR.scores.HML)
mean_ICR_High_perCancer$`Mean ICR score ICR High` = as.numeric(mean_ICR_High_perCancer$`Mean ICR score ICR High`)
mean_ICR_High_perCancer$Cancer = factor(as.character(mean_ICR_perCancer$Cancer), levels = Cancer_order)

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Compare_ICRscore_ICRE_ICRD/Mean_ICR_score_ICRHigh_boxplot_ICR_EnabledvsDisabled.png"), res=600,height=6,width=6,unit="in")
ggplot(mean_ICR_High_perCancer, aes(x = `ICR ID`, y = `Mean ICR score ICR High`)) + 
  stat_compare_means(method = "t.test") +
  scale_color_manual(values = Cancer_color_table$color) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(0.3), aes(colour = Cancer)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.line = element_line(color= "black", size = 0.4),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text = element_text(size = 12, colour = "black"))
dev.off()

t.test(`Mean ICR score ICR High`~`ICR ID`, data = mean_ICR_High_perCancer)

# delta and ratio ICR High versus Low
ICR_High_vs_Low = Cluster_assignment_analysis[-which(Cluster_assignment_analysis$Cancertype == "LAML"), c("Cancertype", "delta.HL.in.HML.classification",
                                                                                                                "ratio.HL.in.HML.classification","ICR_ID")]
colnames(ICR_High_vs_Low) = c("Cancer", "Delta ICR score ICR High vs ICR Low", "Ratio ICR score ICR High vs ICR Low", "ICR ID")
ICR_High_vs_Low$Cancer = factor(as.character(ICR_High_vs_Low$Cancer), levels = Cancer_order)

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Compare_ICRscore_ICRE_ICRD/Delta_ICR_score_ICRHigh_vs_Low_boxplot_ICR_EnabledvsDisabled.png"), res=600,height=6,width=6,unit="in")
ggplot(ICR_High_vs_Low, aes(x = `ICR ID`, y = `Delta ICR score ICR High vs ICR Low`)) + 
  stat_compare_means(method = "t.test") +
  scale_color_manual(values = Cancer_color_table$color) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(0.3), aes(colour = Cancer)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.line = element_line(color= "black", size = 0.4),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text = element_text(size = 12, colour = "black"))
dev.off()

t.test(`Delta ICR score ICR High vs ICR Low`~`ICR ID`, data = ICR_High_vs_Low)

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Compare_ICRscore_ICRE_ICRD/Ratio_ICR_score_ICRHigh_vs_Low_boxplot_ICR_EnabledvsDisabled.png"), res=600,height=6,width=6,unit="in")
ggplot(ICR_High_vs_Low, aes(x = `ICR ID`, y = `Ratio ICR score ICR High vs ICR Low`)) + 
  stat_compare_means(method = "t.test") +
  scale_color_manual(values = Cancer_color_table$color) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(0.3), aes(colour = Cancer)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.line = element_line(color= "black", size = 0.4),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text = element_text(size = 12, colour = "black"))
dev.off()

t.test(`Ratio ICR score ICR High vs ICR Low`~`ICR ID`, data = ICR_High_vs_Low)

