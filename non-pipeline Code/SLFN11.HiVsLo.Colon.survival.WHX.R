
# Setup environment
rm(list=ls())
## dependencies
## install java for xlsx export
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("xlsx","plyr","ggplot2","reshape","survival")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (xlsx) #xlsx needs java installed
library (plyr)
library(survival)
library (reshape)
library (ggplot2)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")

setwd("~/Dropbox/BREAST_QATAR/")

# Set Parameters
Filtersamples <- "Filtered" # altervatives : Filtered , UnFiltered
Surv.cutoff.years <- 10       # SET cut-off
Km.type <- "1vs2vs3vs4"       # SET curve type  - altervatives :1vs2vs3vs4 4vs123 OR 1vs4

# Load data files
load ("./2 DATA/TCGA RNAseq/RNASeq_COAD-hiseq_EDASeq/COAD-hiseq.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
SLFN11.hiseq <- as.data.frame(t(RNASeq.NORM_Log2["SLFN11",,drop=FALSE]))
load ("./2 DATA/TCGA RNAseq/RNASeq_COAD-GA_EDASeq/COAD-GA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
SLFN11.GA <- as.data.frame(t(RNASeq.NORM_Log2["SLFN11",,drop=FALSE]))
rm(RNASeq.NORM_Log2)
Clinical.data.1 <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.COAD-hiseq.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data.1) <- Clinical.data.1[,1]
Clinical.data.1[,1] <-NULL
Clinical.data.2 <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.COAD-GA.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data.2) <- Clinical.data.2[,1]
Clinical.data.2[,1] <-NULL
Clinical.data   <- unique(as.data.frame(rbind(as.matrix(Clinical.data.1),as.matrix(Clinical.data.2))))

# split into quartiles
SLFN11.hiseq.HI <- SLFN11.hiseq[SLFN11.hiseq$SLFN11<quantile(SLFN11.hiseq$SLFN11)["25%"],,drop=FALSE]
SLFN11.hiseq.LO <- SLFN11.hiseq[SLFN11.hiseq$SLFN11>quantile(SLFN11.hiseq$SLFN11)["75%"],,drop=FALSE]
SLFN11.GA.HI <- SLFN11.GA[SLFN11.GA$SLFN11<quantile(SLFN11.GA$SLFN11)["25%"],,drop=FALSE]
SLFN11.GA.LO <- SLFN11.GA[SLFN11.GA$SLFN11>quantile(SLFN11.GA$SLFN11)["75%"],,drop=FALSE]
SLFN11.HI <- rbind(SLFN11.hiseq.HI,SLFN11.GA.HI)
SLFN11.LO <- rbind(SLFN11.hiseq.LO,SLFN11.GA.LO)
rm(SLFN11.GA.HI,SLFN11.GA.LO,SLFN11.hiseq.HI,SLFN11.hiseq.LO)

# Select relevant clinical data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude.post == "No")     # remove excluded patients
} else if  (Filtersamples=="UnFiltered")
{Clinical.data.subset <- Clinical.data
}
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","death_days_to","last_contact_days_to")]  # select relevant data

# Add quantile to clinical data
Clinical.data.subset.TS$Quantile <- "Medium"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% rownames(SLFN11.HI)] <- "High"
Clinical.data.subset.TS$Quantile[rownames(Clinical.data.subset.TS) %in% rownames(SLFN11.LO)] <- "Low"

# High vs Low
Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Quantile %in% c("High","Low"),]
cbPalette <- c("#0000FF","#FF0000")

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
png(paste0("./4 FIGURES/KM curves/ggplot.KM.SLFN11.HiVsLo.TCGA.COAD.",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=c("High","Low") ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for SLFN11 RNASeq HiVsLo"),
     xlabs = "Time in months",
     cbPalette = cbPalette
)
dev.off()
