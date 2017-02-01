#################################################################
###
### This Script PLots Heatmaps based on 
### Consensus Clustering grouping of MA Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/MA/...
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/Heatmaps
###
### Parameters to modified : Geneset, K 
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)//BREAST_QATAR/")
#Dependencies
required.packages <- c("gplots","GMD")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library("GMD")

# Load Data
Geneset <- "DBGS1" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 3             # SET K here

Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/GBM/GBM.TCGA.MA.Affy.k7.ISGS.reps5000/GBM.TCGA.MA.Affy.k7.ISGS.reps5000.k=",K,".consensusClass.csv"),header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
load (paste0("./2 DATA/SUBSETS/ASSEMBLER/GBM/TCGA.GBM.MA.subset.",Geneset,".RData"))
MA.subset <- as.matrix(MA.subset)

## Code to reorder within cluster, expression data not used
load ("./3 ANALISYS/CLUSTERING/MA/GBM/GBM.TCGA.MA.k7.DBGS1.reps5000/ConsensusClusterObject.Rdata")
consensusClusters <- as.factor(ConsensusClusterObject[[K]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
hhc <- ConsensusClusterObject[[K]]$consensusTree
sampleOrder <- consensusClusters[hhc$order]  ## get the order of cluser assignments based on the consensus tree

ConsensusClusterObject.oGE <- t(MA.subset[names(sampleOrder),])

column_annotation <- matrix(" ", nrow = ncol(ConsensusClusterObject.oGE), ncol = 1)
column_annotation[, 1] <- as.character(consensusClusters[names(sampleOrder)])
unique(column_annotation[, 1])
ConsensusClusterObject.oGEz <- (ConsensusClusterObject.oGE - rowMeans(ConsensusClusterObject.oGE))/apply(ConsensusClusterObject.oGE, 1, sd)
quantile(ConsensusClusterObject.oGEz, 0.15); quantile(ConsensusClusterObject.oGEz, 0.85)
ConsensusClusterObject.oGEz[ConsensusClusterObject.oGEz <= quantile(ConsensusClusterObject.oGEz, 0.15)] <- quantile(ConsensusClusterObject.oGEz, 0.15)
ConsensusClusterObject.oGEz[ConsensusClusterObject.oGEz >= quantile(ConsensusClusterObject.oGEz, 0.85)] <- quantile(ConsensusClusterObject.oGEz, 0.85)

#Add cluster assignment to data
MA.subset <- merge (MA.subset,Consensus.class,by="row.names")
row.names(MA.subset) <- MA.subset$Row.names
MA.subset$Row.names <- NULL
MA.subset$PatientID <- NULL

#Rename ICR clusters
Cluster.order <- data.frame(Group=MA.subset[,ncol(MA.subset)], avg=rowMeans (MA.subset[,1:(ncol(MA.subset)-1)]))
Cluster.order <- aggregate(Cluster.order,by=list(Cluster.order$Group),FUN=mean)
Cluster.order <- cbind(Cluster.order[order(Cluster.order$avg),c(2,3)],ICR.name=c("ICR1","ICR2","ICR3"))#,"ICR4"))
Consensus.class$Group[Consensus.class$Group==Cluster.order[1,1]] <- as.character(Cluster.order[1,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[2,1]] <- as.character(Cluster.order[2,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[3,1]] <- as.character(Cluster.order[3,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[4,1]] <- as.character(Cluster.order[4,3])
write.csv (Consensus.class,file=paste0("./3 ANALISYS/CLUSTERING/MA/GBM/GBM.TCGA.MA.k7.DBGS1.reps5000/GBM.TCGA.MA.k7.DBGS1.reps5000.k=",K,".consensusClass.ICR.csv"))     

#Update Cluster names
MA.subset$Group <- NULL
MA.subset <- merge (MA.subset,Consensus.class,by="row.names")
row.names(MA.subset) <- MA.subset$Row.names
MA.subset$Row.names <- NULL
MA.subset$PatientID <- NULL

#ordeing within clusters
MA.subset   <- MA.subset[colnames(ConsensusClusterObject.oGEz),]

#ordering of the clusters
MA.subset <- MA.subset[order(factor(MA.subset$Group,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
MA.subset$Group <- NULL

#re-order the labels
Consensus.class <- Consensus.class[rownames(MA.subset),]

# Heatmap 2 (simple no extra annotations)
patientcolors <- Consensus.class
levels (patientcolors$Group) <- c(levels (patientcolors$Group),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
#patientcolors$Group[patientcolors$Group=="ICR4"] <- "#FF0000"
#patientcolors$Group[patientcolors$Group=="ICR3"] <- "#FFA500"
patientcolors$Group[patientcolors$Group=="ICR3"] <- "#FF0000"
patientcolors$Group[patientcolors$Group=="ICR2"] <- "#00FF00"
patientcolors$Group[patientcolors$Group=="ICR1"] <- "#0000FF"
#patientcolors$Group <- droplevels(patientcolors$Group)
patientcolors <- patientcolors$Group
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100)))
png("./4 FIGURES/Heatmaps/Heatmap.MA.TCGA.GBM.DBGS1.png",res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(MA.subset),
          main = paste0("Heatmap MA - ",Geneset," sel., K=",K),
          col=my.palette,                   #set color sheme RED High, GREEN low
          breaks=my.colors,                #Uncomment and comment key                 
          ColSideColors=patientcolors,      #set goup colors                 
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1.3,cexCol=0.1,margins=c(2,7),
          Colv=FALSE)                       #warning is normal
par(lend = 1)
#legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
legend("topright",legend = c("ICR3","ICR2","ICR1"),
       col = c("red","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()
######################

