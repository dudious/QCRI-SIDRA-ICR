
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
Gene              = "TP53"
ICR.filter        = "ICR"

# Load data files
load ("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.All.DBGS3.FLTR.Mutation.Matrixes.NonSilent.Rdata")
MAPKX.MUT <- rownames(genes.mutations.selected[which(genes.mutations.selected[,"MAP3K1"]==1 | genes.mutations.selected[,"MAP2K4"]==1),])
genes.mutations.selected$MAPKX <- NA
genes.mutations.selected[MAPKX.MUT,"MAPKX"] <- 1
patients.MUT <- rownames(genes.mutations.selected[which(genes.mutations.selected[,Gene]==1),])
patients.WT  <- rownames(genes.mutations.selected[-which(genes.mutations.selected[,Gene]==1),])
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.BSF2.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

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

# Add MUTSTAT and ICR cluster to clinical data
Clinical.data.subset.TS$MUTSTAT <- "Unknown"
Clinical.data.subset.TS$MUTSTAT[rownames(Clinical.data.subset.TS) %in% patients.WT] <- "WT"
Clinical.data.subset.TS$MUTSTAT[rownames(Clinical.data.subset.TS) %in% patients.MUT] <- "MUT"
Clinical.data.subset.TS$Cluster <- Consensus.class$Cluster[match(rownames(Clinical.data.subset.TS),rownames(Consensus.class))]

# MUT vs WT
Clinical.data.subset.TS$Grouping <- Clinical.data.subset.TS$MUTSTAT
Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Grouping %in% c("MUT","WT"),]
label = c("Mutated","Wild-Type")
cbPalette <- c("#0000FF","#FF0000")

if (ICR.filter != "ALL" ){
  Clinical.data.subset.TS$Grouping <- paste0 (Clinical.data.subset.TS$Cluster,"-",Clinical.data.subset.TS$MUTSTAT)
  Clinical.data.subset.TS <- Clinical.data.subset.TS[Clinical.data.subset.TS$Grouping %in% c("ICR1-MUT","ICR1-WT","ICR4-MUT","ICR4-WT"),]
  label = c("ICR1-Mutated","ICR1-Wild-Type","ICR4-Mutated","ICR4-Wild-Type")
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
png(paste0("./4 FIGURES/KM curves/ggplot.KM.",Gene,".MUTvsWT.",ICR.filter,".TCGA.BRCA.",Surv.cutoff.years,"Y.png"),res=600,height=6,width=6,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs=label ,
     ystrataname="Legend",
     main=paste0("Kaplan-Meier Plot for ",Gene," RNASeq MUTvsWT ",ICR.filter),
     xlabs = "Time in months",
     cbPalette = cbPalette
)
dev.off()

mdiff <- survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval <- pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

#TS.Surv$Group <- relevel(as.factor(TS.Surv$Group),"WT")
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
