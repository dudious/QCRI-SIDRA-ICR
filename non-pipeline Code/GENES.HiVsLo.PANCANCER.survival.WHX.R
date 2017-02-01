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
#source ("F:/DropBox Wouter/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")
#setwd("F:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")

# Set Parameters
Filtersamples     = "Filtered"    # altervatives : Filtered , UnFiltered
Surv.cutoff.years = 10            # SET cut-off
Gene              = "BRGS"
GeneIsSignature   = "TRUE"
GeneSignature     = c("TACSTD2","DSP","JUP","DST","FLG","PPL","PKP3")
Division          = "Tertile"
Cancersets        = "ALL"

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
stats.table <- data.frame(stringsAsFactors = FALSE,Cancer = Cancersets,p.value = NA, HR = NA, Lower = NA, Upper = NA,N = NA, UpTert = NA, LowTert = NA)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("ACC","LAML","FPPP","LIHC","UVM")) {next}

# Load data files
load(paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.BIOLINKS.Selected.NORMALIZED.TP_FILTERED_LOG2.RData"))
# Calculate GeneSignature Score
if (GeneIsSignature == "TRUE"){
  RNASeq.Subset <- RNASeq.NORM.TP_Log2[GeneSignature,]
  Gene.expression <- colMeans(RNASeq.Subset) 
}
else {
  Gene.expression <- (RNASeq.NORM.TP_Log2[Gene,])
}
rm(RNASeq.NORM.TP_Log2)
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_BIOLINKS_subset_clinicaldata.csv"),header=TRUE)
if(length(which(is.na(Clinical.data$X)))>0){Clinical.data <- Clinical.data[-which(is.na(Clinical.data$X)),]}
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
dir.create(paste0("./4 FIGURES/KM curves/",Gene),showWarnings = FALSE)

# split into quartiles
if (Division == "Quantile") {
Gene.LO <- Gene.expression[Gene.expression<quantile(Gene.expression)["25%"]]
Gene.HI <- Gene.expression[Gene.expression>quantile(Gene.expression)["75%"]]
}
if (Division == "Tertile") {
  Gene.LO <- Gene.expression[Gene.expression<quantile(Gene.expression,c(.33333,.66666))["33.333%"]]
  Gene.HI <- Gene.expression[Gene.expression>quantile(Gene.expression,c(.33333,.66666))["66.666%"]]
}

# Select relevant clinical data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude.post == "No")     # remove excluded patients
} else if  (Filtersamples=="UnFiltered")
{Clinical.data.subset <- Clinical.data
}
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","days_to_death","days_to_last_follow_up")]  # select relevant data

# Add quantile to clinical data
Clinical.data.subset.TS$Quantile <- "Medium"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% names(Gene.HI)] <- "High"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% names(Gene.LO)] <- "Low"

# High vs Low
Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Quantile %in% c("High","Low"),]
cbPalette <- c("#FF0000","#0000FF")

# time / event object creation
Y <- Surv.cutoff.years * 365
TS.Alive <- subset(Clinical.data.subset.TS[,c(1,3,4)],Clinical.data.subset.TS$vital_status == "alive")
colnames(TS.Alive) <- c("Status","Time","Group")
TS.Alive$Time <- as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] <- Y
TS.Dead <- subset(Clinical.data.subset.TS[,c(1,2,4)],Clinical.data.subset.TS$vital_status == "dead")
colnames(TS.Dead) <- c("Status","Time","Group")
TS.Dead$Time <- as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "alive"
TS.Dead$Time[TS.Dead$Time > Y] <- Y                                                                        
TS.Surv <- rbind (TS.Dead,TS.Alive)
TS.Surv$Time <- as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "dead"
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)             

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png(paste0("./4 FIGURES/KM curves/",Gene,"/ggplot.KM.",Gene,".HiVsLo.",Division,".TCGA.",Cancerset,".",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=c("High","Low") ,
     ystrataname="Legend",
     main=paste0("KM Plot for ",Cancerset,"-",Gene," RNASeq HiVsLo ",Division),
     xlabs = "Time in months",
     cbPalette = cbPalette
)
dev.off()
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

result <- c(Cancerset,signif(pval, 3),signif(exp(mHR.extract@coef),3),CI[1,],length(Gene.expression),length(Gene.HI),length(Gene.LO))
stats.table[stats.table$Cancer == Cancerset,] <- result
}

write.csv(stats.table,file = paste0("./4 FIGURES/KM curves/",Gene,"/stats.",Gene,".",Division,".csv"))
          