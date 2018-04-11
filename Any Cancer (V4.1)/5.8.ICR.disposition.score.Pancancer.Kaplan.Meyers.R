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
Pathways = "ALL"
CancerTYPES = "ALL"
Pathway_skip = ""                                                                                                      # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
Log_file = paste0("./1_Log_Files/", download.method ,"/5.8.Pancancer_Survival_Analysis/5.8.Pancancer_Survival_Analysis_Log_File_",                              # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
ICR_k = "HML_classification"                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Cutoff_HR = 1
ICR_subset = "ICR_disabled"                                                                                                     # "ICR_All", "ICR_enabled", "ICR_disabled" or Type_of_Cancer
ICR_benefit_score_cutoff = 0
ICR_benefit_score_subset = "ICR_benefit_score_Low"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", "Survival_analysis_High_vs_Low_Groups", 
            ICR_k, ".Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefitv1.Rdata"))

if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

# Create folders
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)  
dir.create(paste0("./1_Log_Files/", download.method ,"/5.8.Pancancer_Survival_Analysis/"), showWarnings = FALSE)
cat("This is a log file for Survival Analysis of ",                                                                     # Set-up logfile
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
    paste0("Pathway_skip = ", Pathway_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    "",
    "Clustering",
    file = Log_file,
    append = FALSE, sep= "\n")

ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR)])
ICR_enabled_cancer_samples = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers)]
ICR_disabled_cancer_samples = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers)]

ICR_survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/ACC/SurvivalData/updatedsurvivaldata.csv"))
for (i in 2:length(CancerTYPES)){
  Cancer = CancerTYPES[i]
  if(Cancer == "LAML"){next}
  survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/SurvivalData/updatedsurvivaldata.csv"))
  ICR_survival_data = rbind(ICR_survival_data, survival_data)
}

if(ICR_subset == "ICR_All"){
  Survival_data1 = ICR_survival_data
}
if(ICR_subset == "ICR_enabled"){
  Survival_data1 = ICR_survival_data[which(ICR_survival_data$bcr_patient_barcode %in% substring(ICR_enabled_cancer_samples, 1, 12)),]
}
if(ICR_subset == "ICR_disabled"){
  Survival_data1 = ICR_survival_data[which(ICR_survival_data$bcr_patient_barcode %in% substring(ICR_disabled_cancer_samples, 1, 12)),]
}

if(ICR_benefit_score_subset == "ICR_benefit_score_Low"){
  samples_to_include = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_score < ICR_benefit_score_cutoff)]
  patients_to_include = substring(samples_to_include, 1, 12)
  Survival_data = Survival_data1[which(Survival_data1$bcr_patient_barcode %in% patients_to_include),]
}

if(ICR_benefit_score_subset == "ICR_benefit_score_High"){
  samples_to_include = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_score >= ICR_benefit_score_cutoff)]
  patients_to_include = substring(samples_to_include, 1, 12)
  Survival_data = Survival_data1[which(Survival_data1$bcr_patient_barcode %in% patients_to_include),]
}

# Create folders to save the data
dir.create(paste0("./4_Analysis/",download.method,"/Pan_Cancer"),showWarnings = FALSE)
dir.create(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis"),showWarnings = FALSE)
dir.create(paste0("./5_Figures"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Survival_Plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Survival_Plots/ICR_benefit_score_v1"), showWarnings = FALSE)

Survival_data$ICR_cluster = Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster[match(Survival_data$bcr_patient_barcode,substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers),1,12))]
Survival_data = Survival_data[!is.na(Survival_data$ICR_cluster),]
Survival_data$ICR_cluster = factor(Survival_data$ICR_cluster, levels = c("ICR High", "ICR Medium", "ICR Low")) 
Highest_ICR_group = "ICR High"

Y = Surv_cutoff_years * 365
TS.Alive = Survival_data[Survival_data$vital_status == "Alive", c("vital_status", "days_to_last_followup", "ICR_cluster")]
colnames(TS.Alive) = c("Status","Time", "Group")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_data[Survival_data$vital_status == "Dead", c("vital_status", "days_to_death", "ICR_cluster")]
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
mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR High", "ICR Low"))
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
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Survival_Plots/ICR_benefit_score_v1/",
           "/Kaplan_Meier_ICR_clusters_", "ICR_benefit_score_v1" ,"_", ICR_benefit_score_subset, "_in_", ICR_subset, "_Pancancer.png"),res=600,height=6,width=8,unit="in")                                                                                           # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Group"]),
     ystrataname = NULL,
     main= paste0("Survival curve across ICR groups in all ", ICR_subset, " and ", ICR_benefit_score_subset," samples"),
     xlabs = "Time in months",
     cbPalette = cbPalette,
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)
dev.off()
