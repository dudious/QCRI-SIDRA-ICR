#################################################################
###
### This Script PLots Heatmaps based on 
### Consensus Clustering grouping of ",Cancerset," RNASeq Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/...
### Data is saved :
### ./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/...
### Figures are saved :
### ./4 FIGURES/Heatmaps
###
### Parameters to modified : Geneset, K 
###
#################################################################

# Setup environment
rm(list=ls())
#setwd("~/Dropbox/BREAST_QATAR/")
setwd("/mnt3/wouter/BREAST-QATAR/")
#Dependencies
required.packages <- c("gplots","GMD")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")

# Set Parameters
DL.Method    = "BIOLINKS"     # Choose "ASSEMBLER" or "BIOLINKS"
sample.types = "Selected"     # Alternatives TP , TP_TM , Selected
Cancersets    = "ALL"         # SET Cancertype
BRCA.Filter  = "PCF"          # "PCF" or "BSF" Pancer or Breast specific
Geneset      = "DBGS3"        # SET GENESET and proclustering filter 
K <- 4                        # SET K here

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}

# Load Data
load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))  
RNASeq.subset <- as.matrix(RNASeq.subset)
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000.k=4.consensusClass.csv"),header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
load (paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000/ConsensusClusterObject.Rdata"))

## Code to reorder within cluster, expression data not used
consensusClusters <- as.factor(ConsensusClusterObject[[K]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
hhc <- ConsensusClusterObject[[K]]$consensusTree
sampleOrder <- consensusClusters[hhc$order]                                            ## get the order of cluser assignments based on the consensus tree
ConsensusClusterObject.oGE <- t(RNASeq.subset[names(sampleOrder),])

#Add cluster assignment to data
RNASeq.subset <- merge (RNASeq.subset,Consensus.class,by="row.names")
row.names(RNASeq.subset) <- RNASeq.subset$Row.names
RNASeq.subset$Row.names <- NULL
RNASeq.subset$PatientID <- NULL

#Rename ICR clusters
Cluster.order <- data.frame(Group=RNASeq.subset[,ncol(RNASeq.subset)], avg=rowMeans (RNASeq.subset[,1:(ncol(RNASeq.subset)-1)]))
Cluster.order <- aggregate(Cluster.order,by=list(Cluster.order$Group),FUN=mean)
Cluster.order <- cbind(Cluster.order[order(Cluster.order$avg),c(2,3)],ICR.name=c("ICR1","ICR2","ICR3","ICR4"))
Consensus.class$Group[Consensus.class$Group==Cluster.order[1,1]] <- as.character(Cluster.order[1,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[2,1]] <- as.character(Cluster.order[2,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[3,1]] <- as.character(Cluster.order[3,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[4,1]] <- as.character(Cluster.order[4,3])
write.csv (Consensus.class,paste0(file="./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000.k=4.consensusClass.ICR.csv"))       

#Update Cluster names
RNASeq.subset$Group <- NULL
RNASeq.subset <- merge (RNASeq.subset,Consensus.class,by="row.names")
row.names(RNASeq.subset) <- RNASeq.subset$Row.names
RNASeq.subset$Row.names <- NULL
RNASeq.subset$PatientID <- NULL

#ordeing within clusters (comment out if no reordering within cluster is required)
# RNASeq.subset   <- RNASeq.subset[colnames(ConsensusClusterObject.oGE),]

#ordering of the clusters
RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$Group,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
RNASeq.subset$Group <- NULL

#re-order the labels
Consensus.class <- Consensus.class[rownames(RNASeq.subset),]

# Heatmap 2 (simple no extra annotations)
patientcolors <- Consensus.class
levels (patientcolors$Group) <- c(levels (patientcolors$Group),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
patientcolors$Group[patientcolors$Group=="ICR4"] <- "#FF0000"
patientcolors$Group[patientcolors$Group=="ICR3"] <- "#FFA500"
patientcolors$Group[patientcolors$Group=="ICR2"] <- "#00FF00"
patientcolors$Group[patientcolors$Group=="ICR1"] <- "#0000FF"
#patientcolors$Group <- droplevels(patientcolors$Group)
patientcolors <- patientcolors$Group
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100)))
png(paste0("./4 FIGURES/Heatmaps/ICR_expression/Heatmap.RNASeq.TCGA.",DL.Method,".",Cancerset,".",Geneset,".png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(RNASeq.subset),
          main = paste0("Heatmap RNASeq: ",Cancerset,"-",Geneset),
          col=my.palette,                   #set color sheme RED High, GREEN low
          breaks=my.colors,                                 
          ColSideColors=patientcolors,      #set goup colors                 
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1.3,cexCol=0.1,cex.main=0.5,
          margins=c(2,7),
          Colv=FALSE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()
print(paste0("Clusterassigment and heatmap for ",Cancerset," done."))
}
