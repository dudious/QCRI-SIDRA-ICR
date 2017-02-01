
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
Filtersamples     = "Filtered" # altervatives : Filtered , UnFiltered
Surv.cutoff.years = 10       # SET cut-off
Km.type           = "1vs2vs3vs4"       # SET curve type  - altervatives :1vs2vs3vs4 4vs123 OR 1vs4
Gene              = "CXCL9"
ICR.filter        = "ICR"
Cancerset         = "BRCA.BSF2"
IMS.filter        = "ALL"
Geneset           = "DBGS3.FLTR"
matrix.type       = "NonSilent"

# Load data files
#combined mutation and CNV datamatix
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
#clustering Data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

# select data to plot
selected.data <- merged.matrix[,"CXCL9"]
selected.data[selected.data=="MUT"] <- NA
selected.data <- gsub("MUT;","",as.matrix(selected.data))
selected.data[selected.data=="HOMDEL"] <- "Deleted"
selected.data[selected.data=="AMP"] <- "Amplified"
selected.data <- as.matrix(selected.data)
selected.data[is.na(selected.data)] <- "Normal"

AMP.patients <- names(selected.data[selected.data[,1]=="Amplified",])
DEL.patients <- names(selected.data[selected.data[,1]=="Deleted",])
NOR.patients <- names(selected.data[selected.data[,1]=="Normal",])


Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.BSF2.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
#Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
#Consensus.class <- Consensus.class[,-1]
#colnames (Consensus.class) <- c("Patient_ID","Cluster")
#rownames(Consensus.class) <- Consensus.class[,1]

# Select relevant clinical data
if (Filtersamples=="Filtered"){     
  Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude.post == "No")     # remove excluded patients
} else if  (Filtersamples=="UnFiltered")
{Clinical.data.subset <- Clinical.data
}
if (ICR.filter != "ALL" ){
  selected.ICR <- rownames(Consensus.class[Consensus.class$Cluster %in% c("ICR1","ICR4"),])
  selected.ICR <- selected.ICR[which(selected.ICR %in% rownames(Clinical.data.subset))]
  Clinical.data.subset <- Clinical.data.subset[selected.ICR,]
}
Clinical.data.subset.TS <- Clinical.data.subset[,c("vital_status","death_days_to","last_contact_days_to")]  # select relevant data

# Add CNA and ICR cluster to clinical data
Clinical.data.subset.TS$CNA <- "Unknown"
Clinical.data.subset.TS$CNA[rownames(Clinical.data.subset.TS) %in% AMP.patients] <- "Amplified"
Clinical.data.subset.TS$CNA[rownames(Clinical.data.subset.TS) %in% DEL.patients] <- "Deleted"
Clinical.data.subset.TS$CNA[rownames(Clinical.data.subset.TS) %in% NOR.patients] <- "Normal"
Clinical.data.subset.TS$Cluster <- Consensus.class$Cluster[match(rownames(Clinical.data.subset.TS),rownames(Consensus.class))]

# CNA
Clinical.data.subset.TS$Grouping <- Clinical.data.subset.TS$CNA
Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Grouping %in% c("Amplified"),]
label = c("Amplified")
cbPalette <- c("red","blue")

if (ICR.filter != "ALL" ){
  Clinical.data.subset.TS$Grouping <- paste0 (Clinical.data.subset.TS$Cluster,"-",Clinical.data.subset.TS$CNA)
  Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Grouping %in% c("ICR1-Amplified","ICR1-Deleted","ICR4-Amplified","ICR4-Deleted"),]
  label = c("ICR1-Amplified","ICR1-Deleted","ICR4-Amplified","ICR4-Deleted")
  cbPalette <- c("#0000FF","#FF0000","#8b0000","#00008b")
}

# time / event object creation
Y <- Surv.cutoff.years * 365
TS.Alive <- subset(Clinical.data.subset.TS[,c(1,3,6)],Clinical.data.subset.TS$vital_status == "Alive")
colnames(TS.Alive) <- c("Status","Time","Group")
TS.Alive$Time <- as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] <- Y
TS.Dead <- subset(Clinical.data.subset.TS[,c(1,2,6)],Clinical.data.subset.TS$vital_status == "Dead")
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
png(paste0("./4 FIGURES/KM curves/ggplot.KM.",Gene,".CNA.",ICR.filter,".TCGA.BRCA.",Surv.cutoff.years,".AMP.v2.Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     #ystratalabs=label ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for ",Gene," RNASeq CNA ",ICR.filter),
     xlabs = "Time in months",
     cbPalette = cbPalette
)
dev.off()

mdiff <- survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval <- pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

TS.Surv$Group <- relevel(as.factor(TS.Surv$Group),"Normal")
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
