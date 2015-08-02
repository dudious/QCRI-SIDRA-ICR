#################################################################
###
### This Script PLots kaplan Meier survival curves based on 
### Consensus Clustering grouping of ",Cancerset," RNASeq Data
### From Hiseq and GA data combined , when ovelapping tech is available
### diffrentially clustered samples are remove  
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
#source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)
library(survival)
library(reshape)
library(ggplot2)
library(plyr)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")

# Set Parameters
Cancerset <- "COAD"           # SET Cancertype
Filtersamples <- "Filtered" # altervatives : Filtered , UnFiltered
Geneset <- "DBGS3.FLTR"            # SET GENESET and pruclustering filter 
K <- 4                        # SET K
Surv.cutoff.years <- 10       # SET cut-off
Km.type <- "1vs4"       # SET curve type  - altervatives :1vs2vs3vs4 4vs123 OR 1vs4


# Load data
Parent.Geneset <- substring(Geneset,1,5)
Cancerset.1 <- paste0(substring(Cancerset,1,4),"-hiseq")
Cancerset.2 <- paste0(substring(Cancerset,1,4),"-GA")
## load hiseq data
Consensus.class.1 <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset.1,"/",Cancerset.1,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000/",Cancerset.1,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class.1 <- Consensus.class.1[,-1]
rownames(Consensus.class.1) <- Consensus.class.1$PatientID
Clinical.data.1 <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset.1,".RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data.1) <- Clinical.data.1[,1]
Clinical.data.1[,1] <-NULL
## load GA data
Consensus.class.2 <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset.2,"/",Cancerset.2,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000/",Cancerset.2,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class.2 <- Consensus.class.2[,-1]
rownames(Consensus.class.2) <- Consensus.class.2$PatientID
Clinical.data.2 <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset.2,".RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data.2) <- Clinical.data.2[,1]
Clinical.data.2[,1] <-NULL
## Combine data
Consensus.class <- unique(as.data.frame(rbind(as.matrix(Consensus.class.1),as.matrix(Consensus.class.2))))
Clinical.data   <- unique(as.data.frame(rbind(as.matrix(Clinical.data.1),as.matrix(Clinical.data.2))))
overlap <- unique(Consensus.class$PatientID[which(duplicated(Consensus.class$PatientID))])
Consensus.class <- Consensus.class[-which(Consensus.class$PatientID %in% overlap),]     # remove overlapping samples
rownames(Consensus.class) <- Consensus.class$PatientID
Clinical.data$PatientID <- substring(rownames(Clinical.data),1,12)
#rownames(Clinical.data) <- NULL

# Select relevant clinical data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude.post == "No")     # remove excluded patients
} else if  (Filtersamples=="UnFiltered")
{Clinical.data.subset <- Clinical.data
}

Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","death_days_to","last_contact_days_to")]  # select relevant data
Clinical.data.subset.DFS <- Clinical.data.subset[,c("tumor_status","last_contact_days_to")]                 # select relevant data for Desease Free Survival

# Add Class to clinical data
Clinical.data.subset.TS <- merge(Clinical.data.subset.TS,Consensus.class["Group"],by="row.names",all.x=TRUE, all.y=FALSE)
row.names(Clinical.data.subset.TS) <- Clinical.data.subset.TS$Row.names
Clinical.data.subset.TS$Row.names <- NULL

# Seup for Curve Type
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
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)                                                                # remove patients with less then 1 day follow up time

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png(paste0("./4 FIGURES/KM curves/ggplot.KM.",Km.type,".TCGA.",Cancerset,"-",Filtersamples,".Combined-sources.RNASeq.",Geneset,".k=",K,".",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=Clusters.names ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for ",Parent.Geneset," RNASeq selection"),
     xlabs = "Time in months",
     cbPalette = cbPalette
)
dev.off()



