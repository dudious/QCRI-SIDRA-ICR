#################################################################
###
### This Script PLots Heatmaps based on 
### Master file containing RNAseq Data
### 
### Input data :
### ./3 ANALISYS/MASTER FILES/TCGA.LIHC.RNASeq_subset_ISGS.Master.csv
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/Heatmaps/heatmap.3
###
### Parameters to modify : 
### K , Geneset
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
source ("./1 CODE/R tools/heatmap.3.R")

# Parameters
K <- 4
Geneset <- "ISGS"

# Load Data
Master.data <- read.csv("./3 ANALISYS/MASTER FILES/TCGA.LIHC.RNASeq_subset_ISGS.Master.csv",header=TRUE) # select source data
rownames(Master.data) <- Master.data$X
Master.data$X <- NULL
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_INES_MA.txt") 
gene.list.selected <- as.character(gene.list[which(gene.list[,"Selected_by_DB"]==1),1])
RNASeq.subset <- Master.data[,gene.list.selected]
load ("./3 ANALISYS/CLUSTERING/RNAseq/LIHC/LIHC.TCGA.EDASeq.k7.ISGS.reps5000/ConsensusClusterObject.Rdata")

# Ordering within the clusters
consensusClusters <- as.factor(ConsensusClusterObject[[K]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
hhc <- ConsensusClusterObject[[K]]$consensusTree
sampleOrder <- consensusClusters[hhc$order]  ## get the order of cluser assignments based on the consensus tree
RNASeq.subset <-RNASeq.subset[names(sampleOrder),]
ConsensusClusterObject.oGE <- t(RNASeq.subset[complete.cases(RNASeq.subset),])

# zeta score calculation
ConsensusClusterObject.oGEz <- (ConsensusClusterObject.oGE - rowMeans(ConsensusClusterObject.oGE))/apply(ConsensusClusterObject.oGE, 1, sd)
quantile(ConsensusClusterObject.oGEz, 0.15); quantile(ConsensusClusterObject.oGEz, 0.85)
ConsensusClusterObject.oGEz[ConsensusClusterObject.oGEz <= quantile(ConsensusClusterObject.oGEz, 0.05)] <- quantile(ConsensusClusterObject.oGEz, 0.05)
ConsensusClusterObject.oGEz[ConsensusClusterObject.oGEz >= quantile(ConsensusClusterObject.oGEz, 0.95)] <- quantile(ConsensusClusterObject.oGEz, 0.95)

#ordering of the clusters
ConsensusClusterObject.oGEz <- as.data.frame(t(ConsensusClusterObject.oGEz))
ConsensusClusterObject.oGEz$Cluster <- Master.data$Cluster.ISGS.RNSeq[match(rownames(ConsensusClusterObject.oGEz),rownames(Master.data))]
ConsensusClusterObject.oGEz <- ConsensusClusterObject.oGEz[order(factor(ConsensusClusterObject.oGEz$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
Consensus.class <- ConsensusClusterObject.oGEz[,"Cluster",drop=FALSE]
ConsensusClusterObject.oGEz$Cluster <- NULL

# Heatmap.3
patientcolors <- Consensus.class
patientcolors$Mutation <- Master.data$P53_mutation[match(rownames(patientcolors),rownames(Master.data))]
patientcolors$Grade <- Master.data$tumor_grade[match(rownames(patientcolors),rownames(Master.data))]

levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to cluster
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
levels (patientcolors$Mutation) <- c(levels (patientcolors$Mutation),c("#FF0000","#0000FF"))  #Aply color scheme to mutation
patientcolors$Mutation[patientcolors$Mutation=="TRUE"] <- "#000000"
patientcolors$Mutation[patientcolors$Mutation=="FALSE"] <- "#a9a9a9"
levels (patientcolors$Grade) <- c(levels (patientcolors$Grade),c("#b3b3b3","#787878","#3d3d3d","#020202"))  #Aply color scheme to grade
patientcolors$Grade[patientcolors$Grade=="G1"] <- "#b3b3b3"
patientcolors$Grade[patientcolors$Grade=="G2"] <- "#787878"
patientcolors$Grade[patientcolors$Grade=="G3"] <- "#3d3d3d"
patientcolors$Grade[patientcolors$Grade=="G4"] <- "#020202"
patientcolors$Grade[patientcolors$Grade=="[Not Available]"] <- NA


patientcolors <- as.matrix(patientcolors)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 300)

png("./4 FIGURES/Heatmaps/heatmap.3/Heatmap.3.RNASeq.TCGA.LIHC.ISGS.png",res=600,height=6,width=6,unit="in")     # set filename
heatmap.3(t(ConsensusClusterObject.oGEz),
          main = paste0("Heatmap RNASeq - ",Geneset," sel., K=",K),
          col=my.palette,
          Colv=FALSE,
          ColSideColors=patientcolors,
          labCol=FALSE,
          cexRow=1.3,cexCol=0.1,margins=c(2,7)
          )
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

             
          
          
     
          



