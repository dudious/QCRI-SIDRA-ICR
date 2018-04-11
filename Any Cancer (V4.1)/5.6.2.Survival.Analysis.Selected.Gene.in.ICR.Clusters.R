####################################################################
###
### This Script calculates 
### 
### Input data:
### ("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/")
### Output data are saved as Rdata file:
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/ggkm.R"))

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                     # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
Log_file = paste0("./1_Log_Files/", download.method ,"/5.6.Pancancer_Survival_Analysis/5.6.Pancancer_Survival_Analysis_Log_File_",                              # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
ICR_k = "HML_classification"                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Cutoff_HR = 1
Gene_of_interest = "HLA-E"
ICR_ability = "ICR_All"
ICR_cluster = "ICR_Low"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", "Survival_analysis_High_vs_Low_Groups", 
            ICR_k, ".Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))

if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

# Create folders
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)  

load(paste0("./3_DataProcessing/", download.method, "/ACC/RNASeqData/ACC_gene_RNAseq_normalized_TP_filtered.Rdata"))
Expression.data.all = filtered.norm.RNAseqData

N.sets = length(CancerTYPES)
for (i in 2:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer == "LAML"){next}
  if(Cancer == "SKCM"){
    load(paste0("./3_DataProcessing/Assembler_Panca_Normalized/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  }else{load(paste0("./3_DataProcessing/Assembler_Panca_Normalized/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))}
  filtered.norm.RNAseqData = filtered.norm.RNAseqData[rownames(Expression.data.all),]
  Expression.data.all = cbind(Expression.data.all, filtered.norm.RNAseqData)
}
Expression.data.all = log(Expression.data.all +1, 2)

ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])
ICR_enabled_samples = rownames(ICR_cluster_assignment_allcancers)[ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers]
ICR_disabled_samples = rownames(ICR_cluster_assignment_allcancers)[ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers]
ICR_All_samples = rownames(ICR_cluster_assignment_allcancers)
samples_of_interest = get(paste0(ICR_ability, "_samples"))


ICR_High_samples = rownames(ICR_cluster_assignment_allcancers)[ICR_cluster_assignment_allcancers$HML_cluster == "ICR High"]
ICR_Medium_samples = rownames(ICR_cluster_assignment_allcancers)[ICR_cluster_assignment_allcancers$HML_cluster == "ICR Medium"]
ICR_Low_samples = rownames(ICR_cluster_assignment_allcancers)[ICR_cluster_assignment_allcancers$HML_cluster == "ICR Low"]
ICR_cluster_of_interest = get(paste0(ICR_cluster, "_samples"))

Expression.data.subset = Expression.data.all[, which(colnames(Expression.data.all) %in% samples_of_interest & colnames(Expression.data.all) %in% ICR_cluster_of_interest)]

Expression.data.subset = Expression.data.subset[Gene_of_interest, , drop = FALSE]

Expression.data.subset.z.score = Expression.data.subset 
for(j in 1: nrow(Expression.data.subset.z.score))  {
  Expression.data.subset.z.score[j,] = (Expression.data.subset[j,]-mean(Expression.data.subset[j,]))/sd(Expression.data.subset[j,]) # z-score the enrichment matrix
}

sHc = hclust(ddist <- dist(t(Expression.data.subset.z.score)), method = "ward.D2")
plot(sHc,labels=FALSE)
ICR_cluster_assignment_allcancers$gene_cluster = cutree(sHc,k = 3)[match(rownames(ICR_cluster_assignment_allcancers),names(cutree(sHc,k = 3)))]
ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[colnames(Expression.data.subset),]
ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[colnames(Expression.data.subset.z.score),]
ICR_cluster_assignment_allcancers$gene_score = Expression.data.subset.z.score[1,]
gene_cluster_means = aggregate(gene_score~gene_cluster,data=ICR_cluster_assignment_allcancers,FUN=mean)
gene_cluster_means = gene_cluster_means[order(gene_cluster_means$gene_score),]
gene_cluster_means$gene_cluster_name = c(paste0(Gene_of_interest, " Low"), paste0(Gene_of_interest, " Medium"),paste0(Gene_of_interest, " High"))
ICR_cluster_assignment_allcancers$gene_cluster_name = gene_cluster_means$gene_cluster_name[match(ICR_cluster_assignment_allcancers$gene_cluster,
                                                                                                 gene_cluster_means$gene_cluster)]

ICR_survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/ACC/SurvivalData/updatedsurvivaldata.csv"))
for (i in 2:length(CancerTYPES)){
  Cancer = CancerTYPES[i]
  if(Cancer == "LAML"){next}
  survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/SurvivalData/updatedsurvivaldata.csv"))
  ICR_survival_data = rbind(ICR_survival_data, survival_data)
}

if(ICR_ability == "ICR_All"){
  Survival_data = ICR_survival_data
}
if(ICR_ability == "ICR_enabled"){
  Survival_data = ICR_survival_data[which(ICR_survival_data$bcr_patient_barcode %in% substring(ICR_enabled_samples, 1, 12)),]
}
if(ICR_ability == "ICR_disabled"){
  Survival_data = ICR_survival_data[which(ICR_survival_data$bcr_patient_barcode %in% substring(ICR_disabled_samples, 1, 12)),]
}

Survival_data = Survival_data[which(Survival_data$bcr_patient_barcode %in% substring(ICR_cluster_of_interest, 1, 12)),]

# Create folders to save the data
dir.create(paste0("./4_Analysis/",download.method,"/Pan_Cancer"),showWarnings = FALSE)
dir.create(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis"),showWarnings = FALSE)
dir.create(paste0("./5_Figures"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Survival_Plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Survival_Plots/Genes_of_interest"), showWarnings = FALSE)


Survival_data$gene_cluster = ICR_cluster_assignment_allcancers$gene_cluster_name[match(Survival_data$bcr_patient_barcode,
                                                                                       substring(rownames(ICR_cluster_assignment_allcancers), 1, 12))]
Survival_data = Survival_data[!is.na(Survival_data$gene_cluster),]
Survival_data$gene_cluster = factor(Survival_data$gene_cluster, levels = c(paste0(Gene_of_interest, " High"), 
                                                                           paste0(Gene_of_interest, " Medium"),
                                                                           paste0(Gene_of_interest, " Low")))
Highest_group = paste0(Gene_of_interest, " High")

Y = Surv_cutoff_years * 365
TS.Alive = Survival_data[Survival_data$vital_status == "Alive", c("vital_status", "days_to_last_followup", "gene_cluster")]
colnames(TS.Alive) = c("Status","Time", "Group")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_data[Survival_data$vital_status == "Dead", c("vital_status", "days_to_death", "gene_cluster")]
colnames(TS.Dead) = c("Status","Time", "Group")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

#TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
TS.Surv[,"Group"] = as.factor(TS.Surv[,"Group"])

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c(paste0(Gene_of_interest, " High"),
                                                                                              paste0(Gene_of_interest, " Low")))
mHR.extract = extract.coxph(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)

PLOT_P = round(p[2],3)
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3)
PLOT_CI1 = CI[2,1]
PLOT_CI2 = CI[2,2]

# plots
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Survival_Plots/Genes_of_interest/", "Kaplan_Meier_",
           Gene_of_interest, "_", ICR_ability, "_", ICR_cluster, "_Pancancer.png"),res=600,height=6,width=9,unit="in")                                                                                           # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Group"]),
     ystrataname = NULL,
     main= paste0("Survival curve across ", Gene_of_interest, " expression groups in ",  ICR_ability, " cancers in ", ICR_cluster),
     xlabs = "Time in months",
     cbPalette = cbPalette,
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)
dev.off()

#####