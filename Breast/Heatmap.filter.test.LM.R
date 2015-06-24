# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
library("plyr")
source ("./1 CODE/R tools/heatmap.3.R")

# Parameters
K <- 4
Cancerset <- "BRCA"
Geneset <- "DBGS1.FLTR.LMF" # SET GENESET HERE !!!!!!!!!!!!!!
Parent.Geneset <- substring(Geneset,1,5)

# Load Data
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Parent.Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)
mode(RNASeq.subset) = "numeric"
RNASeq.subset <- cbind(RNASeq.subset,IS=rowMeans(RNASeq.subset))
RNASeq.subset <- RNASeq.subset[order(RNASeq.subset[,"IS"]),]

Consensus.class.LMF <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
rownames(Consensus.class.LMF) <- Consensus.class.LMF$X
Consensus.class.LMF$X <- NULL

Geneset <- "DBGS1" # SET GENESET HERE !!!!!!!!!!!!!!
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
rownames(Consensus.class) <- Consensus.class$X
Consensus.class$X <- NULL

#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL

#create clusterassignmnet table
Consensus.class <- Consensus.class[rownames(RNASeq.subset),]
Consensus.class.LMF <- Consensus.class.LMF[rownames(RNASeq.subset),]

Consensus.class$Group.LMF <- Consensus.class.LMF$Group[match(Consensus.class$PatientID,Consensus.class.LMF$PatientID)]
colnames(Consensus.class) <- c("PatientID","Cluster","Cluster.LMF")

#add clinical data to clusterassignment table - Patientcolors
patientcolors <- Consensus.class
patientcolors$Histology <- ClinicalData.subset$histological_type[match(patientcolors$PatientID,rownames(ClinicalData.subset))]
patientcolors$Subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(patientcolors$PatientID,rownames(ClinicalData.subset))]

# zeta score calculation (NOT USED)
RNASeq.subset <- RNASeq.subset[,-ncol(RNASeq.subset)]
RNASeq.subset.zeta <- (RNASeq.subset - rowMeans(RNASeq.subset))/apply(RNASeq.subset, 1, sd)
quantile(RNASeq.subset.zeta, 0.15); quantile(RNASeq.subset.zeta, 0.85)
RNASeq.subset.zeta[RNASeq.subset.zeta <= quantile(RNASeq.subset.zeta, 0.05)] <- quantile(RNASeq.subset.zeta, 0.05)
RNASeq.subset.zeta[RNASeq.subset.zeta >= quantile(RNASeq.subset.zeta, 0.95)] <- quantile(RNASeq.subset.zeta, 0.95)

# Heatmap.3
patientcolors <- patientcolors[,-1]
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to cluster
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000" #red
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500" #orange
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00" #green
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF" #blue

levels (patientcolors$Cluster.LMF) <- c(levels (patientcolors$Cluster.LMF),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to Cluster.LMF
patientcolors$Cluster.LMF[patientcolors$Cluster.LMF=="ICR4"] <- "#FF0000"
patientcolors$Cluster.LMF[patientcolors$Cluster.LMF=="ICR3"] <- "#FFA500"
patientcolors$Cluster.LMF[patientcolors$Cluster.LMF=="ICR2"] <- "#00FF00"
patientcolors$Cluster.LMF[patientcolors$Cluster.LMF=="ICR1"] <- "#0000FF"
patientcolors$Cluster.LMF[is.na(patientcolors$Cluster.LMF)] <- "#000000" #black

levels (patientcolors$Histology) <- c(levels (patientcolors$Histology),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to Histology
patientcolors$Histology[patientcolors$Histology=="Infiltrating Ductal Carcinoma"] <- "#FF0000"
patientcolors$Histology[patientcolors$Histology=="Infiltrating Lobular Carcinoma"] <- "#FFA500"
patientcolors$Histology[patientcolors$Histology=="Mixed Histology (please specify)"] <- "#00FF00"
patientcolors$Histology[patientcolors$Histology=="Other  specify"] <- "#0000FF"
patientcolors$Histology[patientcolors$Histology=="Metaplastic Carcinoma"] <- "#0000FF"
patientcolors$Histology[patientcolors$Histology=="Mucinous Carcinoma"] <- "#0000FF"
patientcolors$Histology[patientcolors$Histology=="Infiltrating Carcinoma NOS"] <- "#0000FF"
patientcolors$Histology[patientcolors$Histology=="Medullary Carcinoma"] <- "#0000FF"
patientcolors$Histology[patientcolors$Histology=="[Not Available]"] <- NA
patientcolors$Histology[is.na(patientcolors$Histology)] <- "#000000"

levels (patientcolors$Subtype) <- c(levels (patientcolors$Subtype),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to Subtype
patientcolors$Subtype[patientcolors$Subtype=="Luminal A"] <- "#FF0000"
patientcolors$Subtype[patientcolors$Subtype=="Luminal B"] <- "#FFA500"
patientcolors$Subtype[patientcolors$Subtype=="HER2-enriched"] <- "#00FF00"
patientcolors$Subtype[patientcolors$Subtype=="Basal-like"] <- "#0000FF"
patientcolors$Subtype[patientcolors$Subtype=="Normal-like"] <- NA
patientcolors$Subtype[is.na(patientcolors$Subtype)] <- "#000000"

colnames(patientcolors) <- c("No Filter","LM Filter","Histology","Subtype")


patientcolors <- as.matrix(patientcolors)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100)))


png("./4 FIGURES/Heatmaps/heatmap.3/Heatmap.3.RNASeq.TCGA.BRCA.DBGS1.FLTR.LM.filtertest.png",res=600,height=6,width=12,unit="in")     # set filename
heatmap.3(t(RNASeq.subset),
          main = paste0("Heatmap RNASeq - ",Geneset," sel., K=",K),
          col=my.palette,
          breaks=my.colors,
          Colv=FALSE,
          ColSideColors=patientcolors,
          labCol=FALSE,
          scale="row",
          cexRow=1.3,cexCol=0.1,margins=c(2,7)
)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()
