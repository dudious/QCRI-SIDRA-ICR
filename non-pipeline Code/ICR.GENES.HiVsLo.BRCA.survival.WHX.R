
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
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")

# Set Parameters
Filtersamples     = "Filtered"    # altervatives : Filtered , UnFiltered
Surv.cutoff.years = 10            # SET cut-off
Gene              = "ICRscore"
Division          = "Quantile"

# Load data files
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
imm.score <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BRCA.DBGS3.csv")
RNASeq.NORM_Log2.table <- as.data.frame(t(RNASeq.NORM_Log2))
RNASeq.NORM_Log2.table$ICRscore <- imm.score$unscaled.IS[match(rownames(RNASeq.NORM_Log2.table),imm.score$X)]
RNASeq.NORM_Log2 <-t(RNASeq.NORM_Log2.table)
rm(RNASeq.NORM_Log2.table)
Gene.expression <- (RNASeq.NORM_Log2[Gene,])
rm(RNASeq.NORM_Log2)
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.BSF2.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL

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
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","death_days_to","last_contact_days_to")]  # select relevant data

# Add quantile to clinical data
Clinical.data.subset.TS$Quantile <- "Medium"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% names(Gene.HI)] <- "High"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% names(Gene.LO)] <- "Low"

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

# plots
png(paste0("./4 FIGURES/KM curves/ggplot.KM.",Gene,".HiVsLo.",Division,".TCGA.BRCA.",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=c("High","Low") ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for ",Gene," RNASeq HiVsLo ",Division),
     xlabs = "Time in months",
     cbPalette = cbPalette
)
dev.off()

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
