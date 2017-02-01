
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
Surv.cutoff.years = 10            # SET cut-off
Gene              = "IS"
Division          = "Quantile"

# Load data files
load ("./2 DATA/SUBSETS/LM.Dataset/LM.Dataset.MA.subset.DBGS3.RData")
imm.score <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.LM.Data.csv")
MA.subset.table <- as.data.frame(t(MA.subset))
load ("./2 DATA/LM.BRCA/LM.Dataset.split.fixed.Rdata")
MA.subset <- MA.subset[-1955,]
MA.subset.table$Gene <- as.character(Gene.Meta.data$Symbol[match(rownames(MA.subset.table),Gene.Meta.data$Affy_Probe_ID)])
MA.subset.table.bygene <- aggregate(MA.subset.table[],list(MA.subset.table$Gene),FUN=mean)
rownames(MA.subset.table.bygene) <- MA.subset.table.bygene$Group.1
MA.subset.table.bygene$Group.1 <- NULL
MA.subset.table.bygene <- as.data.frame(t(MA.subset.table.bygene))
MA.subset.table.bygene <- MA.subset.table.bygene[-c(1955,1956),]
MA.subset.table.bygene$IS <- imm.score$unscaled.IS[match(rownames(MA.subset.table.bygene),imm.score$X)]

#select relevant expression data
Gene.expression.table <- MA.subset.table.bygene[,Gene,drop = FALSE]
Gene.expression <- Gene.expression.table[,1]
names(Gene.expression) <- rownames(Gene.expression.table)
Clinical.data.subset.TS <- Sample.Meta.data[,c("DMFS_10y_time","DMFS_10y_event")]  # select relevant data

# split into quartiles
if (Division == "Quantile") {
Gene.LO <- Gene.expression[Gene.expression<quantile(Gene.expression)["25%"]]
Gene.HI <- Gene.expression[Gene.expression>quantile(Gene.expression)["75%"]]
}
if (Division == "Tertile") {
  Gene.LO <- Gene.expression[Gene.expression<quantile(Gene.expression,c(.33333,.66666))["33.333%"]]
  Gene.HI <- Gene.expression[Gene.expression>quantile(Gene.expression,c(.33333,.66666))["66.666%"]]
}


# Add quantile to clinical data
Clinical.data.subset.TS$Quantile <- "Medium"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% names(Gene.HI)] <- "High"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% names(Gene.LO)] <- "Low"

# High vs Low
Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Quantile %in% c("High","Low"),]
cbPalette <- c("#FF0000","#0000FF")

# time / event object creation
Y <- Surv.cutoff.years * 365
TS.Alive <- subset(Clinical.data.subset.TS[,c(2,1,3)],Clinical.data.subset.TS$DMFS_10y_event == 0)
colnames(TS.Alive) <- c("Status","Time","Group")
TS.Alive$Time <- as.numeric(as.character(TS.Alive$Time))*365
TS.Alive$Time[TS.Alive$Time > Y] <- Y
TS.Dead <- subset(Clinical.data.subset.TS[,c(2,1,3)],Clinical.data.subset.TS$DMFS_10y_event == 1)
colnames(TS.Dead) <- c("Status","Time","Group")
TS.Dead$Time <- as.numeric(as.character(TS.Dead$Time))*365
TS.Dead$Status[which(TS.Dead$Time> Y)] = 0
TS.Dead$Time[TS.Dead$Time > Y] <- Y                                                                        
TS.Surv <- rbind (TS.Dead,TS.Alive)
TS.Surv$Time <- as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == 1
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)             

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png(paste0("./4 FIGURES/KM curves/ggplot.KM.",Gene,".HiVsLo.",Division,".LM.BRCA.",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=c("High","Low") ,
     ystrataname="Legend",
     main=paste0("KM Plot for ",Gene," RNASeq LM HiVsLo ",Division),
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

