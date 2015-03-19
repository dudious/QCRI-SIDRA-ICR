#################################################################
###
### This Script PLots kaplan Meier survival curves based on 
### Consensus Clustering grouping of MA Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/MA/...
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
Geneset <- "27G" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 3             # SET K here
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/BRCA.TCGA.MA.k4.",Geneset,".reps1000/BRCA.TCGA.MA.k4.",Geneset,".reps1000.k=",K,".consensusClass.csv"),header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
Clinical.data <- read.csv ("./3 ANALISYS/CLINICAL DATA/Agilent_subset_clinicaldata.csv",header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
#Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude == "No")                                # remove excluded patients
Clinical.data.subset <- Clinical.data
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","death_days_to","last_contact_days_to")]  # select relevant data
Clinical.data.subset.DFS <- Clinical.data.subset[,c("tumor_status","last_contact_days_to")]                 # select relevant data fro desease free survival

# Add Class to clinical data
Clinical.data.subset.TS <- merge(Clinical.data.subset.TS,Consensus.class["Group"],by="row.names",all.x=TRUE, all.y=FALSE)
row.names(Clinical.data.subset.TS) <- Clinical.data.subset.TS$Row.names
Clinical.data.subset.TS$Row.names <- NULL

# time / event object creation
TS.Alive <- subset(Clinical.data.subset.TS[,c(1,3,4)],Clinical.data.subset.TS$vital_status == "Alive")
colnames(TS.Alive) <- c("Status","Time","Group")
TS.Dead <- subset(Clinical.data.subset.TS[,c(1,2,4)],Clinical.data.subset.TS$vital_status == "Dead")
colnames(TS.Dead) <- c("Status","Time","Group")
TS.Surv <- rbind (TS.Dead,TS.Alive)
TS.Surv$Time <- as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"

# survival curve
msurv <- Surv(TS.Surv$Time/7, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots

png(paste0("./4 FIGURES/KM curves/plot.KM.MA.R1000.k=",K,".",Geneset,".png"),res=600,height=6,width=6,unit="in")     # set filename
plot(mfit,
     conf.int=F,
     mark.time=T,
     col=c("red","blue"),
     lwd=0.8,
     xlab="Time in weeks",
     ylab="Survival probability",
     main=paste0("Kaplan-Meier plot for MA ",Geneset)
     )
dev.off()

png(paste0("./4 FIGURES/KM curves/ggplot.KM.MA.R1000.k=",K,".",Geneset,".png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=52,
     ystratalabs=paste("Group",(1:K)),
     main=paste0("Kaplan-Meier Plot for ",Geneset," MA selection"),
     xlabs = "Time in weeks",
     )
dev.off()
