# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
library("plyr")
source ("./1 CODE/R tools/heatmap.3.R")

# Parameters
K <- 4
Cancerset <- "BRCA"
Geneset <- "DBGS1.FLTR.NMS" # SET GENESET HERE !!!!!!!!!!!!!!
Parent.Geneset <- substring(Geneset,1,5)

# Load Data
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Parent.Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)
mode(RNASeq.subset) = "numeric"
RNASeq.subset <- cbind(RNASeq.subset,IS=rowMeans(RNASeq.subset))
RNASeq.subset <- RNASeq.subset[order(RNASeq.subset[,"IS"]),]

Consensus.class.NMS <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
rownames(Consensus.class.NMS) <- Consensus.class.NMS$X
Consensus.class.NMS$X <- NULL

Geneset <- "DBGS1.FLTR.N" # SET GENESET HERE !!!!!!!!!!!!!!
Consensus.class.N <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
rownames(Consensus.class.N) <- Consensus.class.N$X
Consensus.class.N$X <- NULL

Geneset <- "DBGS1" # SET GENESET HERE !!!!!!!!!!!!!!
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
rownames(Consensus.class) <- Consensus.class$X
Consensus.class$X <- NULL


#create clusterassignmnet table
Consensus.class <- Consensus.class[rownames(RNASeq.subset),]
Consensus.class.N <- Consensus.class.N[rownames(RNASeq.subset),]
Consensus.class.NMS <- Consensus.class.NMS[rownames(RNASeq.subset),]

Consensus.class$Group.N <- Consensus.class.N$Group[match(Consensus.class$PatientID,Consensus.class.N$PatientID)]
Consensus.class$Group.NMS <- Consensus.class.NMS$Group[match(Consensus.class$PatientID,Consensus.class.NMS$PatientID)]
colnames(Consensus.class) <- c("PatientID","Cluster","Cluster.N","Cluster.NMS")

# zeta score calculation (NOT USED)
RNASeq.subset <- RNASeq.subset[,-ncol(RNASeq.subset)]
RNASeq.subset.zeta <- (RNASeq.subset - rowMeans(RNASeq.subset))/apply(RNASeq.subset, 1, sd)
quantile(RNASeq.subset.zeta, 0.15); quantile(RNASeq.subset.zeta, 0.85)
RNASeq.subset.zeta[RNASeq.subset.zeta <= quantile(RNASeq.subset.zeta, 0.05)] <- quantile(RNASeq.subset.zeta, 0.05)
RNASeq.subset.zeta[RNASeq.subset.zeta >= quantile(RNASeq.subset.zeta, 0.95)] <- quantile(RNASeq.subset.zeta, 0.95)

# Heatmap.3
patientcolors <- Consensus.class[,-1]
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to cluster
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
levels (patientcolors$Cluster.N) <- c(levels (patientcolors$Cluster.N),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to Cluster.NMS
patientcolors$Cluster.N[patientcolors$Cluster.N=="ICR4"] <- "#FF0000"
patientcolors$Cluster.N[patientcolors$Cluster.N=="ICR3"] <- "#FFA500"
patientcolors$Cluster.N[patientcolors$Cluster.N=="ICR2"] <- "#00FF00"
patientcolors$Cluster.N[patientcolors$Cluster.N=="ICR1"] <- "#0000FF"
patientcolors$Cluster.N[is.na(patientcolors$Cluster.N)] <- "#000000"
levels (patientcolors$Cluster.NMS) <- c(levels (patientcolors$Cluster.NMS),c("#FF0000","#FFA500","#00FF00","#0000FF","#000000"))  #Aply color scheme to Cluster.NMS
patientcolors$Cluster.NMS[patientcolors$Cluster.NMS=="ICR4"] <- "#FF0000"
patientcolors$Cluster.NMS[patientcolors$Cluster.NMS=="ICR3"] <- "#FFA500"
patientcolors$Cluster.NMS[patientcolors$Cluster.NMS=="ICR2"] <- "#00FF00"
patientcolors$Cluster.NMS[patientcolors$Cluster.NMS=="ICR1"] <- "#0000FF"
patientcolors$Cluster.NMS[is.na(patientcolors$Cluster.NMS)] <- "#000000"


patientcolors <- as.matrix(patientcolors)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100)))


png("./4 FIGURES/Heatmaps/heatmap.3/Heatmap.3.RNASeq.TCGA.BRCA.DBGS1.filtertest.png",res=600,height=6,width=12,unit="in")     # set filename
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
