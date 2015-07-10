#################################################################
###
### This Script PLots Heatmaps based on 
### Consensus Clustering clustering of RNASeq Data and mutation data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/RNAseq/...
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/Heatmaps
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
# Dependencies
required.packages <- c("gplots","plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library("plyr")

## Parameters
Cancerset <- "BLCA"
Geneset = "DBGS3.FLTR"
matrix.type = "Any"       # Alterantives "Any" , "Missense"
plot.type = "high"         # Alterantives "low" , "high" , "373genes"  (NOT WORKING)

# Load Data
load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

# select data to plot
if (plot.type == "low"){allmuts.mutatedgenes <- genes.mutations.low}
if (plot.type == "high"){allmuts.mutatedgenes <- genes.mutations.high}
if (plot.type == "373genes"){allmuts.mutatedgenes <- genes.mutations.373genes}
allmuts.mutatedgenes[is.na(allmuts.mutatedgenes)] = 0

#merge data
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,Consensus.class,by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL

#ordeing rows (patients)
allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg

#calculate frequency (mean)
allmuts.mutatedgenes.mean <- ddply(allmuts.mutatedgenes,.(Cluster),colwise(mean))                             # Calculate frequency (mean)

#generate numeric mutation matrix
allmuts.mutatedgenes$Cluster <- NULL
allmuts.mutatedgenes$Patient_ID <- NULL
allmuts.mutatedgenes.mean$Cluster <- NULL
allmuts.mutatedgenes <- as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean <- as.matrix(allmuts.mutatedgenes.mean)
mode(allmuts.mutatedgenes) <- "numeric"
mode(allmuts.mutatedgenes.mean) <- "numeric"

#ordering columns (genes)
allmuts.mutatedgenes.sd <- as.data.frame(apply(allmuts.mutatedgenes.mean,2,sd))                             # Calculate frequency SD 
colnames (allmuts.mutatedgenes.sd) <- c("SD")
allmuts.mutatedgenes.sd <- allmuts.mutatedgenes.sd[order(allmuts.mutatedgenes.sd$SD),,drop = FALSE]

allmuts.mutatedgenes.mean <- as.data.frame(allmuts.mutatedgenes.mean[,rownames(allmuts.mutatedgenes.sd)])   # order mutaion frequency by SD
#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,rownames(allmuts.mutatedgenes.sd)])            # order mutation count by SD freq
allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,order(colnames(allmuts.mutatedgenes))])         # order mutation count alphabeticaly
Consensus.class<-as.data.frame(Consensus.class[rownames(allmuts.mutatedgenes),])                      # sort cluster asignments like mutation matrix

#enforce numeric mutation matrix
allmuts.mutatedgenes = as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean = as.matrix(allmuts.mutatedgenes.mean)
mode(allmuts.mutatedgenes)="numeric"
mode(allmuts.mutatedgenes.mean)="numeric"

# Heatmap
patientcolors <- Consensus.class
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
#patientcolors$Cluster <- droplevels(patientcolors$cluster)
patientcolors <- as.character(patientcolors$Cluster)

my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".reordered_alpha.png"),res=600,height=6,width=25,unit="in")     # set filename
heatmap.2(allmuts.mutatedgenes,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=patientcolors,                        # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=3,
          margins=c(10,2),
          labRow=FALSE,
          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
          )
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()


# Heatmap 2
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors <- c(seq(0,0.01,length=100),seq(0.01,0.05,length=100),seq(0.05,1,length=100))
patientcolors <- c("#0000FF","#00FF00","#FFA500","#FF0000")
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",Geneset,".Mean.Mutation.HeatMap.",matrix.type,".",plot.type,".png"),res=600,height=5,width=50,unit="in")     # set filename
par(mar=c(1,1,1,1))
heatmap.2(allmuts.mutatedgenes.mean,
          main = "HeatMap-MutatedGenes frequency",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          breaks=my.colors,                                   # set manual color gradient of color scheme
          RowSideColors=patientcolors,                        # set goup colors
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          #density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=3,
          margins=c(15,15),
          labRow=FALSE,
          Colv=FALSE,Rowv=FALSE                               # reorder row/columns by dendogram
          )
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1)
dev.off()



