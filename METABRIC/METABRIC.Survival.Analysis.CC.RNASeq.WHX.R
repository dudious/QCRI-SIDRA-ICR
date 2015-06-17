#################################################################
###
### This Script PLots kaplan Meier survival curves based on 
### Consensus Clustering grouping of METABRIC RNASeq Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/RNAseq/...
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/KM Curves
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("survival","reshape","ggplot2","plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
#required.packages.BioC <- c("y")
#missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
#source("http://bioconductor.org/biocLite.R")
#if(length(missing.packages)) biocLite(missing.packages)
library(survival)
library(reshape)
library(ggplot2)
library(plyr)

source ("./1 CODE/R tools/ggkm.R")

# Load data
Geneset <- "16G.I" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 4             # SET K here
DATASET <- 1       # SET METABRIC DATASET (1 or 2)
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.METABRIC.DATA.",DATASET,".k4.",Geneset,".reps1000/BRCA.METABRIC.DATA.",DATASET,".k4.",Geneset,".reps1000.k=",K,".consensusClass.csv"),header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames (Consensus.class) <- Consensus.class[,1]
load ("./2 DATA/METABRIC/FROM GABRIELE ZOPPOLI/GeneExpression_METABRIC2012_a.RData")                             # SET METABRIC a or b
if (!is.null(clinical_annotation$Sample.ID)) {rownames (clinical_annotation) <- clinical_annotation$Sample.ID}
 

#Clinical.data.subset <- Clinical.data
Clinical.data.subset.TS <- clinical_annotation[,c("T","last_follow_up_status")]  # select relevant data
#Clinical.data.subset.DFS <- Clinical.data.subset[,c("tumor_status","last_contact_days_to")]                 # select relevant data for Desease Free Survival

# Add Class to clinical data
Clinical.data.subset.TS <- merge(Clinical.data.subset.TS,Consensus.class["Group"],by="row.names",all.x=TRUE, all.y=FALSE)
row.names(Clinical.data.subset.TS) <- Clinical.data.subset.TS$Row.names
Clinical.data.subset.TS$Row.names <- NULL

# time / event object creation
TS.Alive <- subset(Clinical.data.subset.TS,Clinical.data.subset.TS$last_follow_up_status == "a")
colnames(TS.Alive) <- c("Time","Status","Group")
TS.Dead <- subset(Clinical.data.subset.TS,Clinical.data.subset.TS$last_follow_up_status %in% c("d","d-d.s.","d-o.c."))
colnames(TS.Dead) <- c("Time","Status","Group")
TS.Surv <- rbind (TS.Dead,TS.Alive)
TS.Surv$Time <- as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status != "a"

# survival curve
msurv <- Surv(TS.Surv$Time/7, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png(paste0("./4 FIGURES/KM curves/METABRIC.",DATASET,".plot.KM.RNASeq.excl.1000R.k=",K,".",Geneset,".png"),res=600,height=6,width=6,unit="in")     # set filename
plot(mfit,
     conf.int=F,
     mark.time=T,
     col=c("red","blue"),
     lwd=0.8,
     xlab="Time in weeks",
     ylab="Survival probability",
     main=paste0("Kaplan-Meier plot ",Geneset)
     )
dev.off()

png(paste0("./4 FIGURES/KM curves/METABRIC.",DATASET,".ggplot.KM.RNASeq.excl.1000R.k=",K,".",Geneset,".png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=52,
     ystratalabs=paste("Group",(1:K)),
     main=paste0("Kaplan-Meier Plot for ",Geneset," RNASeq selection"),
     xlabs = "Time in weeks",
     )
dev.off()

