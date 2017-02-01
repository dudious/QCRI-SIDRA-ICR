#################################################################
###
### This Script PLots kaplan Meier survival curves based on 
### Consensus Clustering grouping of ",Cancerset," RNASeq Data
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
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")
#setwd("/mnt3/wouter/BREAST-QATAR/")
#Dependencies
required.packages <- c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("reshape")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)
library(reshape)
library(ggplot2)
library(plyr)
<<<<<<< HEAD
library(texreg)
library(survival)
=======
library(survival)
library(texreg)
>>>>>>> 5a74652a1b18b5d5e16e21f7fddb7033a8e8d4b5

source ("~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")
#source ("/mnt3/wouter/QCRI-SIDRA-ICR/R tools/ggkm.R")

# Set Parameters
DL.Method         = "BIOLINKS"     #Choose "ASSEMBLER" or "BIOLINKS"
sample.types      = "Selected"     #Alternatives TP , TP_TM , Selected
Cancersets        = "BRCA"     # SET Cancertype (include Filter type for BRCA.BSF of BRCA.PCF)
Filtersamples     = "Filtered" # altervatives : Filtered , UnFiltered
Geneset           = "DBGS3"     # SET GENESET and pruclustering filter 
K                 = 4                      # SET K
Surv.cutoff.years = 10     # SET cut-off
Km.type           = "1vs4"           # SET curve type  - altervatives :1vs2vs3vs4 4vs123 OR 1vs4

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}

# Load data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000/",
                                   Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
rownames(Consensus.class) <- Consensus.class$PatientID
if (DL.Method == "ASSEMBLER") {
  if (Cancerset == "BRCA"){
    if (Filter == "TRUE"){
      Cancerset <- paste0(Cancerset,".",BRCA.Filter)
    }
  }
  Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv")) 
}
if (DL.Method =="BIOLINKS") {
  Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_",DL.Method,"_subset_clinicaldata.csv")) 
}

if (length(which(is.na(Clinical.data[,1])))>0){Clinical.data <- Clinical.data[-which(is.na(Clinical.data[,1])),]}
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL

# Post-Clustering Filter Data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude.post == "No")     # remove excluded patients
} else if  (Filtersamples=="UnFiltered")
  {Clinical.data.subset <- Clinical.data
}

#Select data for survival analysis
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","days_to_death","days_to_last_follow_up")]  # select relevant data
#Clinical.data.subset.DFS <- Clinical.data.subset[,c("tumor_status","last_contact_days_to")]                 # select relevant data for Desease Free Survival

# Add Class to clinical data
Clinical.data.subset.TS <- merge(Clinical.data.subset.TS,Consensus.class["Group"],by="row.names",all.x=TRUE, all.y=FALSE)
row.names(Clinical.data.subset.TS) <- Clinical.data.subset.TS$Row.names
Clinical.data.subset.TS$Row.names <- NULL

if (Km.type =='4vs123') {
  # ICR4 vs ICR123
  levels (Clinical.data.subset.TS$Group) <- c(levels (Clinical.data.subset.TS$Group), "ICR123")
  Clinical.data.subset.TS$Group[which(Clinical.data.subset.TS$Group %in% c("ICR1","ICR2","ICR3"))] <- "ICR123"
  cbPalette <- c("#FF0000","#000000")
  
  } else if (Km.type =='1vs2vs3vs4') {
    cbPalette <- c("#0000FF","#00FF00","#FFA500","#FF0000")
  } else if (Km.type =='1vs4') {
    # ICR4 vs ICR1
    Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Group %in% c("ICR1","ICR4"),]
    cbPalette <- c("#0000FF","#FF0000")
  }
Clinical.data.subset.TS$Group <- droplevels(Clinical.data.subset.TS$Group) 
Clusters.names <- levels(Clinical.data.subset.TS$Group)

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
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)                                                                # remove patients with less then 1 day follow up time

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png(paste0("./4 FIGURES/KM curves/",DL.Method,"/ggplot.KM.",Km.type,".TCGA.",Cancerset,"-",Filtersamples,".RNASeq.",Geneset,".k=",K,".",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
#dev.new()
ggkm(mfit,
     timeby=12,
     ystratalabs=Clusters.names ,
     ystrataname="Legend",
     main=paste0("KM curve for ",Cancerset," - ",Geneset,"."),
     xlabs = "Time in months",
     cbPalette = cbPalette
     )
dev.off()
print (Cancerset)

mdiff <- survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval <- pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

TS.Surv$Group <- relevel(TS.Surv$Group,"ICR4")
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
}
