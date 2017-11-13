#################################################################
###
### This Script PLots Heatmaps based on 
### Master file containing RNAseq Data
### 
### Input data :
### ./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_ISGS.Master.csv
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/Heatmaps/heatmap.3
###
### Parameters to modify : 
### K , Geneset, Cancerset
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
source ("./1 CODE/R tools/heatmap.3.R")

# Parameters
Cancerset <- "BRCA"
Geneset <- "DBGS1"
K <- 4

# Load Data
Master.data <- read.csv(paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.csv"),header=TRUE) # select source data
rownames(Master.data) <- Master.data$X
Master.data$X <- NULL
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_v2.2.txt")                                 # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,Geneset]==1),1])
RNASeq.subset <- Master.data[,gene.list.selected]
load (paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/ConsensusClusterObject.Rdata"))

# Ordering within the clusters
consensusClusters <- as.factor(ConsensusClusterObject[[K]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
hhc <- ConsensusClusterObject[[K]]$consensusTree
sampleOrder <- consensusClusters[hhc$order]  ## get the order of cluser assignments based on the consensus tree
#RNASeq.subset <-RNASeq.subset[names(sampleOrder),]
ConsensusClusterObject.oGE <- t(RNASeq.subset[complete.cases(RNASeq.subset),])

# zeta score calculation
ConsensusClusterObject.oGEz <- (ConsensusClusterObject.oGE - rowMeans(ConsensusClusterObject.oGE))/apply(ConsensusClusterObject.oGE, 1, sd)
quantile(ConsensusClusterObject.oGEz, 0.15); quantile(ConsensusClusterObject.oGEz, 0.85)
ConsensusClusterObject.oGEz[ConsensusClusterObject.oGEz <= quantile(ConsensusClusterObject.oGEz, 0.05)] <- quantile(ConsensusClusterObject.oGEz, 0.05)
ConsensusClusterObject.oGEz[ConsensusClusterObject.oGEz >= quantile(ConsensusClusterObject.oGEz, 0.95)] <- quantile(ConsensusClusterObject.oGEz, 0.95)

#ordering of the clusters
ConsensusClusterObject.oGEz <- as.data.frame(t(ConsensusClusterObject.oGEz))
order.command <- paste0("ConsensusClusterObject.oGEz$Cluster <- Master.data$Cluster.",Geneset,".RNSeq[match(rownames(ConsensusClusterObject.oGEz),rownames(Master.data))]")
eval(parse(text=order.command))
ConsensusClusterObject.oGEz <- ConsensusClusterObject.oGEz[order(factor(ConsensusClusterObject.oGEz$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
Consensus.class <- ConsensusClusterObject.oGEz[,"Cluster",drop=FALSE]
ConsensusClusterObject.oGEz$Cluster <- NULL

# Heatmap.3
patientcolors <- Consensus.class
patientcolors$Mutation <- Master.data$GOF_mutation[match(rownames(patientcolors),rownames(Master.data))]
patientcolors$Neo <- Master.data$history_neoadjuvant_treatment[match(rownames(patientcolors),rownames(Master.data))]

levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to cluster
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
levels (patientcolors$Mutation) <- c(levels (patientcolors$Mutation),c("#FF0000","#0000FF"))  #Aply color scheme to mutation
patientcolors$Mutation[patientcolors$Mutation=="TRUE"] <- "#000000"
patientcolors$Mutation[patientcolors$Mutation=="FALSE"] <- "#a9a9a9"
levels (patientcolors$Neo) <- c(levels (patientcolors$Neo),c("#000000","#a9a9a9"))  #Aply color scheme to grade
patientcolors$Neo[patientcolors$Neo=="No"] <- "#a9a9a9"
patientcolors$Neo[patientcolors$Neo=="Yes"] <- "#000000"
patientcolors$Neo[patientcolors$Neo=="[Not Available]"] <- NA


patientcolors <- as.matrix(patientcolors)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 300)

png(paste0("./4 FIGURES/Heatmaps/heatmap.3/Heatmap.3.RNASeq.TCGA.",Cancerset,".",Geneset,".png"),res=600,height=6,width=6,unit="in")     # set filename
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

             
          
          
     
          



