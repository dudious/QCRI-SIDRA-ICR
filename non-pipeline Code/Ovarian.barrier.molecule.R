# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

#load extra data
#RNASeq data
load("./2 DATA/TCGA RNAseq/RNASeq_OV_EDASeq/OV.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
Barrier.genes <- c("DST","DSC3","DSP","JUP","TACSTD2","PKP3","PPL","FLG")
TH.1.genes <- c("CXCR3","GNLY","CCL4","CCL3","CD4","IFNG","IRF1","CD8A","CCL5","CD3E","GZMA","GZMB","CCR5","CXCL11","GBP1","CXCL10","CXCL9")
Barrier.genes %in% rownames(RNASeq.NORM_Log2)
TH.1.genes %in% rownames(RNASeq.NORM_Log2)
#clustering data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/OV/OV.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/OV.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#Th.1 rowmeans
RNASeq.th1.subset <- t(RNASeq.NORM_Log2[TH.1.genes,])
RNASeq.th1.means <- data.frame(z.score=rowMeans(RNASeq.th1.subset))
RNASeq.th1.means <- RNASeq.th1.means[order(RNASeq.th1.means$z.score),,drop=FALSE]
RNASeq.th1.means$cluster <- Consensus.class$Cluster[match(rownames(RNASeq.th1.means),rownames(Consensus.class))]
RNASeq.th1.means <- RNASeq.th1.means[-which(is.na(RNASeq.th1.means$cluster)),]

#Barrier genes heatmap
RNASeq.barrier.subset <- t(RNASeq.NORM_Log2[Barrier.genes,])

#lookup cluster
cluster <- Consensus.class[rownames(RNASeq.barrier.subset),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#drop unclusterd samples
cluster <- cluster[-which(is.na(cluster)),,drop=FALSE]
RNASeq.barrier.subset <- RNASeq.barrier.subset[rownames(cluster),]

#color matrix
clustercolors <- as.character(cluster$Cluster)
color.matrix <- as.matrix(rbind (clustercolors))
rownames(color.matrix) <- c("Cluster")

#barrier heatmap
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))

dev.new(with=6,height=6)
heatmap.3(t(RNASeq.barrier.subset),
          main = "barrier genes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          ColSideColors=t(color.matrix),                         # set goup colors
          breaks = my.colors,
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1,cexCol=0.1,margins=c(9,9),
          Colv=TRUE,Rowv=TRUE                           # reorder row/columns by dendogram
)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)

#Barrier - TH1 ordered genes heatmap
RNASeq.combined.subset <- t(RNASeq.NORM_Log2[c(Barrier.genes,TH.1.genes),])

#lookup cluster
cluster <- Consensus.class[rownames(RNASeq.combined.subset),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#drop unclusterd samples
cluster <- cluster[-which(is.na(cluster)),,drop=FALSE]
RNASeq.combined.subset <- RNASeq.combined.subset[rownames(cluster),]

#color matrix
clustercolors <- as.character(cluster$Cluster)
color.matrix <- as.matrix(rbind (clustercolors))
rownames(color.matrix) <- c("Cluster")

#combined heatmap
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))

dev.new(with=6,height=6)
heatmap.3(t(RNASeq.combined.subset),
          main = "barrier/TH1 genes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          ColSideColors=t(color.matrix),                         # set goup colors
          breaks = my.colors,
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1,cexCol=0.1,margins=c(9,9),
          Colv=TRUE,Rowv=TRUE                           # reorder row/columns by dendogram
)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)

#re-order heatmap
RNASeq.combined.subset <- RNASeq.combined.subset[rownames(RNASeq.th1.means),]
RNASeq.combined.subset <- RNASeq.combined.subset[,c(Barrier.genes,TH.1.genes)]
#lookup cluster
cluster <- Consensus.class[rownames(RNASeq.combined.subset),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#color matrix
clustercolors <- as.character(cluster$Cluster)
color.matrix <- as.matrix(rbind (clustercolors))
rownames(color.matrix) <- c("Cluster")

#combined reordered heatmap
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))

dev.new(with=6,height=6)
heatmap.3(t(RNASeq.combined.subset),
          main = "barrier/TH1 genes reordered",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          ColSideColors=t(color.matrix),                         # set goup colors
          breaks = my.colors,
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1,cexCol=0.1,margins=c(9,9),
          Colv=FALSE,Rowv=FALSE                           # reorder row/columns by dendogram
)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)

