# Setup environment
rm(list=ls())
## dependencies
## install java for xlsx export
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("plyr","ggplot2","reshape","survival","texreg")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (plyr)
library(survival)
library (reshape)
library (ggplot2)
library (texreg)
source ("~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline//")

# Set Parameters
Codepath          = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
download.method   = "TCGA_Assembler"   
Filtersamples     = "UnFiltered"    # altervatives : Filtered , UnFiltered
Surv.cutoff.years = 10            # SET cut-off
Gene              = "PLs.ALL"
GeneIsSignature   = "TRUE"
GeneSignature     = c("PLCB1","PLCB2","PLCB3","PLCB4","PLCD1","PLCD3","PLCD4","PLCE1","PLCG1","PLCG2",
                      "PLCH1","PLCH2","PLCL1","PLCL2","PLCZ1","PLD1","PLD2")
  
  #c("PLCB1", "PLCB2", "PLCB3", "PLCB4", "PLCG1", "PLCG2", "PLCD1", "PLCD3", "PLCD4", "PLCE1", "PLCH1", "PLCH2", "PLCL1", "PLCL2")

Division          = "Tertile"
Cancersets        = "ALL"

# DO ALL
TCGA.cancersets <- read.csv (paste0(Codepath,"Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$cancerType))
}
N.sets = length(Cancersets)
stats.table <- data.frame(stringsAsFactors = FALSE,Cancer = Cancersets,p.value = NA, HR = NA, Lower = NA, Upper = NA,N = NA, UpTert = NA, LowTert = NA)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("ACC","LAML","FPPP","READ","LIHC","UVM")) {next}

# Load data files
Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancerset,"/RNASeqData")
if(Cancerset == "SKCM"){
    load(paste0(Cancer_path, "/",Cancerset, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    load(paste0(Cancer_path, "/", Cancerset, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
}
# Calculate GeneSignature Score
if (GeneIsSignature == "TRUE"){
  RNASeq.Subset <- log(filtered.norm.RNAseqData[GeneSignature,],2)
  Gene.expression <- colMeans(RNASeq.Subset) 
} else {
  Gene.expression <- log(filtered.norm.RNAseqData[Gene,],2)
}
rm(filtered.norm.RNAseqData)
Clinical.data <- read.csv (paste0("./3_DataProcessing/TCGA_Assembler/",Cancerset,"/SurvivalData/updatedsurvivaldata.csv"),header=TRUE,stringsAsFactors = FALSE)[,-1]
Clinical.data[Clinical.data=="[Not Available]"] <- NA
Clinical.data[Clinical.data=="[Not Applicable]"] <- NA
if(length(which(is.na(Clinical.data$vital_status)))>0){Clinical.data <- Clinical.data[-which(is.na(Clinical.data$vital_status)),]}
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
dir.create(paste0("./5_Figures/Kaplan_Meier_Plots/HiLoGenes/",Gene),showWarnings = FALSE)

# split into quartiles
if (Division == "Quantile") {
  Gene.LO <- unique(substring(names(Gene.expression[Gene.expression<=quantile(Gene.expression,c(.25,.75))["25%"]]),1,12))
  Gene.HI <- unique(substring(names(Gene.expression[Gene.expression>=quantile(Gene.expression,c(.25,.75))["75%"]]),1,12))
}
if (Division == "Tertile") {
  Gene.LO <- unique(substring(names(Gene.expression[Gene.expression<=quantile(Gene.expression,c(.33333,.66666))["33.333%"]]),1,12))
  Gene.HI <- unique(substring(names(Gene.expression[Gene.expression>=quantile(Gene.expression,c(.33333,.66666))["66.666%"]]),1,12))
}

# Select relevant clinical data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude.post == "No")     # remove excluded patients
} else if  (Filtersamples=="UnFiltered"){
  Clinical.data.subset <- Clinical.data
}
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","days_to_death","days_to_last_followup")]  # select relevant data

# Add quantile to clinical data
Clinical.data.subset.TS$Quantile <- "Medium"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% Gene.HI] <- "High"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% Gene.LO] <- "Low"

# High vs Low
Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Quantile %in% c("High","Low"),]
cbPalette <- c("#FF0000","#0000FF")

# time / event object creation
Y <- Surv.cutoff.years * 365
TS.Alive <- subset(Clinical.data.subset.TS[,c(1,3,4)],Clinical.data.subset.TS$vital_status == "Alive")
colnames(TS.Alive) <- c("Status","Time","Group")
TS.Alive$Time <- as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] <- Y
TS.Dead <- subset(Clinical.data.subset.TS[,c(1,2,4)],Clinical.data.subset.TS$vital_status == "Dead")
colnames(TS.Dead) <- c("Status","Time","Group")
TS.Dead$Time <- as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] <- Y                                                                        
TS.Surv <- rbind (TS.Dead,TS.Alive)
TS.Surv$Time <- as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)             

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# stats 
print(Cancerset)
mdiff <- survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval <- pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

TS.Surv$Group <- relevel(as.factor(TS.Surv$Group),"High")
mHR <- coxph(formula = msurv ~ Group ,data = TS.Surv)
mHR.extract <-extract(mHR, include.aic = TRUE,
                      include.rsquared = TRUE, include.maxrs=TRUE,
                      include.events = TRUE, include.nobs = TRUE,
                      include.missings = TRUE, include.zph = TRUE)
HRtxt <- paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta <- coef(mHR)
se   <- sqrt(diag(mHR$var))
p    <- 1 - pchisq((beta/se)^2, 1)
CI   <- confint(mHR)
CI   <- round(exp(CI),2)
print(pvaltxt)
print(HRtxt)
print(CI)

PLOT_HR  = signif(exp(mHR.extract@coef),3)
PLOT_P   = signif(pval, 3)
PLOT_CI1 = CI[1,1]
PLOT_CI2 = CI[1,2]

result <- c(Cancerset,signif(pval, 3),signif(exp(mHR.extract@coef),3),CI[1,],length(Gene.expression),length(Gene.HI),length(Gene.LO))
stats.table[stats.table$Cancer == Cancerset,] <- result

# plots
png(paste0("./5_Figures/Kaplan_Meier_Plots/HiLoGenes/",Gene,"/ggplot.KM.",Gene,".HiVsLo.",Division,".TCGA.",Cancerset,".",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=c("High","Low") ,
     ystrataname="Legend",
     main=paste0("KM Plot for ",Cancerset,"-",Gene," RNASeq HiVsLo ",Division),
     xlabs = "Time in months",
     cbPalette = cbPalette,
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2
)
dev.off()

}

write.csv(stats.table,file = paste0("./5_Figures/Kaplan_Meier_Plots/HiLoGenes/",Gene,"/stats.",Gene,".",Division,".csv"))
          