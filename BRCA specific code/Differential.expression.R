# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

#install missing packages
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("edgeR")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
#source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

#Load the packages
library(gplots)
library(edgeR)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

#load the functions

TimeUse <-function(func){ 
  time1 <- proc.time() 
  result <-func
  time2 <- proc.time()
  show(timeelapsed <- time2-time1)
}

edgeRDEA4GSEAOff <- function(geneCounts, groups, baseLine){
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  DEA <- DGEList(geneCounts, group = groups)
  DEA <- calcNormFactors(DEA)
  keep <- rowSums(cpm(DEA)>1)>=10
  DEA <- DEA[keep,]
  DEA$samples$lib.size <- colSums(DEA$counts)
  message("Estimating GLM Common Dispersion")
  TimeUse(DEA <- estimateGLMCommonDisp(DEA, design = design))
  DEA <- estimateGLMTrendedDisp(DEA, design = design)
  message("Estimating GLM Tagwise Dispersion")
  TimeUse(DEA <- estimateGLMTagwiseDisp(DEA, design = design))
  message("Fitting the GLM Model")
  TimeUse(fit <- glmFit(DEA, design = design))
  contrast <- rep(1, 2)
  blId <- which(colnames(design) == baseLine)
  contrast[blId] <- -1
  message("Maximum Likelihood Estimate")
  TimeUse(test <- glmLRT(fit, contrast = contrast))
  ans <- topTags(test, n = nrow(DEA$counts))
  ans <- ans$table
  diffGenes <- abs(ans[, "logFC"]) >= 1 & ans[, "FDR"] <= 0.01
  detags <- rownames(ans)[diffGenes]
  pdf("smearPlot.pdf")
  plotSmear(test,de.tags=detags)
  abline(h=c(-1,1),col="blue")
  dev.off()
  return(ans)
}

##############################################################################

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
Signif.tresh = 0.05
fch.tresh = 1
method = "edgeR"
IMS_filter = "Luminal"
Grouping = "ICR1v4"
RNASeq_filter = 0.25


#load MAPKMUT data
MAPK.MUT.state <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
#load and subset RNAseq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")
selectedGenes <- which(rowMeans(RNASeq.NORM) >= quantile(rowMeans(RNASeq.NORM), RNASeq_filter))
RNASeq.NORM <- RNASeq.NORM[selectedGenes,]
RNASeq.subset <- as.data.frame(t(RNASeq.NORM[,])) #select all / no subset used
#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
#Clustering data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#Filter luminal only
if (IMS_filter=="Luminal"){
  subtype <- ClinicalData.subset[rownames(RNASeq.subset),"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
  luminal <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A" | subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B",,drop=FALSE]
  RNASeq.subset <- RNASeq.subset[which(rownames(RNASeq.subset) %in% rownames(luminal)),]
}

##edgeR
  #MAPK MUT vs WT
  #Add MAPK status (MUT/WT) (non-silent)
  RNASeq.subset$MAPK <- MAPK.MUT.state$MAP2K4.MAP3K1[match(rownames(RNASeq.subset),MAPK.MUT.state$Sample)]
  RNASeq.subset <- RNASeq.subset[-which(is.na(RNASeq.subset$MAPK)),]
  
  groups <- rep("",nrow(RNASeq.subset))
  names(groups)<- rownames(RNASeq.subset)
  G1<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]=="MUT",])
  G2<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]=="WT",])
  groups[G1]<-"G1"
  groups[G2]<-"G2"
  groups<-as.factor(groups)
  
  RNASeq.subset <- t(RNASeq.subset[,-ncol(RNASeq.subset)])  #remove grouping column
  
  system.time(MAPKs.MUT.DEGs <- edgeRDEA4GSEAOff(RNASeq.subset[,names(groups)], groups = groups,baseLine = "G2"))
  MAPKs.MUT.DEGs<-MAPKs.MUT.DEGs[order((MAPKs.MUT.DEGs$logFC),decreasing = T),]
  
  RNASeq.subset <- t(RNASeq.subset)
                     
  #ICR1 vs ICR4
  #Add Cluster data
  RNASeq.subset$Cluster <- Consensus.class$Cluster[match(rownames(RNASeq.subset),rownames(Consensus.class))]
  RNASeq.subset <- RNASeq.subset[-which(is.na(RNASeq.subset$Cluster)),]
  
  groups <- rep("",nrow(RNASeq.subset))
  names(groups)<- rownames(RNASeq.subset)
  G1<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]=="ICR1",])
  G2<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]=="ICR4",])
  groups[G1]<-"G1"
  groups[G2]<-"G2"
  groups<-as.factor(groups)
  
  RNASeq.subset <- t(RNASeq.subset[,-ncol(RNASeq.subset)])  #remove grouping column
  
  system.time(ICR1v4.DEGs <- edgeRDEA4GSEAOff(RNASeq.subset[,names(groups)], groups = groups,baseLine = "G2"))
  ICR1v4.DEGs<-ICR1v4.DEGs[order((ICR1v4.DEGs$logFC),decreasing = T),]
  
  RNASeq.subset <- t(RNASeq.subset)
  
  save(MAPKs.MUT.DEGs,ICR1v4.DEGs,file=paste0("./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEG_",IMS_filter,".RDATA"))
  load(file=paste0("./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEG_",IMS_filter,".RDATA"))
  
  