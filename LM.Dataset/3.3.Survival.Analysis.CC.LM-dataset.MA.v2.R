#################################################################
###
### This Script PLots kaplan Meier survival curves based on 
### Consensus Clustering grouping of ",Cancerset," MA Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/MA/...
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
required.packages <- c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
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
library(texreg)
source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")

# Set Parameters
Cancerset         <- "LM.Dataset"
Filtersamples     <- "Filtered" # altervatives : Filtered , UnFiltered
Geneset           <- "DBGS3"          # SET GENESET HERE !!!!!!!!!!!!!!
K                 <- 4                      # SET K here
Surv.cutoff.years <- 10     # SET cut-off here
Km.type           <- "4vs123"     # altervatives :1vs2vs3vs4 4vs123 OR 1vs4

# Load data
#Clusters.names <- rep(paste0("ICR",1:K))
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/",Cancerset,".MA.k7.",
                                   Geneset,".reps5000/",Cancerset,".MA.k7.",
                                   Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
rownames(Consensus.class) <- Consensus.class$PatientID
load ("./2 DATA/LM.BRCA/LM.Dataset.split.fixed.Rdata")

Clinical.data.subset.TS <- Sample.Meta.data[,c("DMFS_10y_time","DMFS_10y_event")]  # select relevant data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset.TS <- Sample.Meta.data[,c("DMFS_10y_time","DMFS_10y_event","PAM50")]
  #replace NA with "Unknown"
  levels(Clinical.data.subset.TS$PAM50)<-c(levels(Clinical.data.subset.TS$PAM50),"Unknown")
  Clinical.data.subset.TS[which(is.na(Clinical.data.subset.TS$PAM50)),"PAM50"] <- "Unknown"
  # exclude Unknown subtype N=26
  Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$PAM50!="Unknown",]
  # exclude normal-like N=257
  Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$PAM50!="Normal",]
  Clinical.data.subset.TS<-Clinical.data.subset.TS[,-ncol(Clinical.data.subset.TS)]
}       
                             

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
TS.Surv <- Clinical.data.subset.TS
colnames(TS.Surv) <- c("Time","Status","Group")
TS.Surv$Time <- TS.Surv$Time
TS.Surv$Time <- TS.Surv$Time * 365
TS.Surv <- subset(TS.Surv,TS.Surv$Time > 1)
TS.Surv$Status <- TS.Surv$Status == 1

# survival curve
msurv <- Surv(TS.Surv$Time/30.4, TS.Surv$Status)
mfit <- survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# plots
png(paste0("./4 FIGURES/KM curves/ggplot.KM.",Km.type,".",Cancerset,"-",Filtersamples,".MA.",Geneset,".k=",K,".",Surv.cutoff.years,"Y.v4.png"),res=600,height=6,width=6,unit="in")  # set filename
#dev.new()
ggkm(mfit,
     timeby=12,
     ystratalabs=Clusters.names ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for ",Geneset," MA selection"),
     xlabs = "Time in months",
     cbPalette = cbPalette
     )
dev.off()

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
