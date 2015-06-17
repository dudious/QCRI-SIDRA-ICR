#################################################################
###
### This Script PLots kaplan Meier survival curves based on 
### Consensus Clustering grouping of RNASeq Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/RNAseq/...
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/KM Curves
###
### Parameters : Geneset, K , for merged clusters change Cluster.names
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("survival","reshape","ggplot2","plyr","Rcpp","colorspace")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("reshape")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)
library(survival)
library(reshape)
library(ggplot2)
library(plyr)

source ("./1 CODE/R tools/ggkm.R")

# Load data
Geneset <- "ISGS" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 4            # SET K here
Y <- 10           # SET CUT-OFF here in years
Clusters.names <- rep(paste0("ICR",1:K))
Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LGG/LGG.TCGA.EDASeq.k7.ISGS.reps5000/LGG.TCGA.EDASeq.k7.ISGS.reps5000.k=4.consensusClass.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
Clinical.data <- read.csv ("./3 ANALISYS/CLINICAL DATA/TCGA.LGG.RNASeq_subset_clinicaldata.csv",header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude == "No")                              # remove excluded patients
#Clinical.data.subset <- Clinical.data
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","death_days_to","last_contact_days_to")]  # select relevant data
Clinical.data.subset.DFS <- Clinical.data.subset[,c("tumor_status","last_contact_days_to")]                 # select relevant data for Desease Free Survival

# Add Class to clinical data
Clinical.data.subset.TS <- merge(Clinical.data.subset.TS,Consensus.class["Group"],by="row.names",all.x=TRUE, all.y=FALSE)
row.names(Clinical.data.subset.TS) <- Clinical.data.subset.TS$Row.names
Clinical.data.subset.TS$Row.names <- NULL

# ICR4 vs ICR123
#levels (Clinical.data.subset.TS$Group) <- c(levels (Clinical.data.subset.TS$Group), "ICR123")
#Clinical.data.subset.TS$Group[which(Clinical.data.subset.TS$Group %in% c("ICR1","ICR2","ICR3"))] <- "ICR123"

# ICR4 vs ICR1
#Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Group %in% c("ICR1","ICR4"),]

# time / event object creation
Y <- Y*365
TS.Alive <- subset(Clinical.data.subset.TS[,c(1,3,4)],Clinical.data.subset.TS$vital_status == "Alive")#21
colnames(TS.Alive) <- c("Status","Time","Group")
TS.Alive$Time <- as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] <- Y#402
TS.Dead <- subset(Clinical.data.subset.TS[,c(1,2,4)],Clinical.data.subset.TS$vital_status == "Dead")#60
colnames(TS.Dead) <- c("Status","Time","Group")
TS.Dead$Time <- as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] <- Y        # set cut off point for survival#87
TS.Surv <- rbind (TS.Dead,TS.Alive)
TS.Surv$Time <- as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)  # remove patients with less then 1 day follow up time
#485 patients

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)#30.4 is avg 365/12
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png("./4 FIGURES/KM curves/ggplot.KM.TCGA.LGG.RNASeq.10Y.png",res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=Clusters.names ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for ",Geneset," RNASeq selection"),
     xlabs = "Time in months",
     )
dev.off()



