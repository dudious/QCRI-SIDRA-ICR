# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
## dependencies
required.packages <- c("plotrix")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(plotrix)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

#parameters
Cancerset <- "COAD-HiSeq"     
Geneset = "DBGS3.FLTR"  # SET GENESET HERE !!!!!!!!!!!!!!

#gene ordere derived from clustering DEG MAPK MUT/WT in luminal
up.genes.order <- c("TAOK2","TP53","MAPK3","MAP3K1","MAPT","HSPA1A","FLNB","TAOK3","CRK","RPS6KA2",
                    "MAP2K4","DUSP5","CACNA1D","MAPK8","RASGRP1","CACNA1G")
down.genes.order <- c("CACNG6","CACNA1B","CACNA2D3","FASLG","RASGRF1","JUN","JUND","DUSP16","PPM1B",
                      "SOS1","FGF12","RASGRP2","PRKCB","MAP4K1","PTPN7","GADD45G","DDIT3","DUSP8",
                      "DUSP10","FGFR4","FGF14","FGF13","MAP2K6","DUSP2")
#RNASeq data
load(paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
RNASeq.subset.UP <- RNASeq.NORM_Log2[up.genes.order,]
RNASeq.subset.DOWN <- RNASeq.NORM_Log2[down.genes.order,]
#reverse DOWN expression
RNASeq.subset.DOWN <- -(RNASeq.subset.DOWN)

#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL

#clustering data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#rowmeans & order
means.UP <-  as.data.frame(colMeans(RNASeq.subset.UP))
colnames(means.UP) <- "UP"
means.DOWN <- as.data.frame(colMeans(RNASeq.subset.DOWN))
colnames(means.DOWN) <- "DOWN"
#combined z-score
means.BOTH <- means.UP
means.BOTH$DOWN <- means.DOWN$DOWN[match(rownames(means.BOTH),rownames(means.DOWN))]
colnames(means.BOTH) <- c("UP","DOWN")
means.BOTH$UP.rank<-rank(means.BOTH$UP)
means.BOTH$DOWN.rank<-rank(means.BOTH$DOWN)
means.BOTH$BOTH <- rowMeans(means.BOTH[,c(1,2)])
means.BOTH$BOTH.rank <- rowMeans(means.BOTH[,c(3,4)])
means.BOTH <- means.BOTH[order(means.BOTH$BOTH.rank),]
means.BOTH$Cluster <- Consensus.class$Cluster[match(rownames(means.BOTH),Consensus.class$Patient_ID)]
#Drop data without cluster assignment
means.BOTH<-means.BOTH[-which(is.na(means.BOTH$Cluster)),]
#Re-order according to combined z-score
RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(means.BOTH)]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(means.BOTH)]

#lookup cluster
cluster <- Consensus.class[colnames(RNASeq.subset.UP),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#cluster filter
cluster <- cluster[cluster$Cluster=="#FF0000"|cluster$Cluster=="#0000FF",,drop=FALSE]
means.BOTH <- means.BOTH[rownames(cluster),]
RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(cluster)]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(cluster)]

#color matrix
clustercolors <- as.character(cluster$Cluster)
MAPK.color.scale <-  color.scale(means.BOTH$BOTH.rank,extremes=c("#ffc1cb","#ff284d"),na.color=NA)
color.matrix <- as.matrix(rbind (clustercolors,MAPK.color.scale))
rownames(color.matrix) <- c("Cluster","MAPK.color.scale")

#reorder rows (genes)
RNASeq.subset.UP <- RNASeq.subset.UP[up.genes.order,]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[down.genes.order,]

#heatmap UP
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))

dev.new(width=6, height=6)
#png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.TCGA.",filter.label,"SIGN.UP.rankmixorder.ICR1vs4.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.3(RNASeq.subset.UP,
          main = paste0("SIGN.UP.MUTvsWT.in.",Cancerset),
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
          #margins=c(2,8),
          Colv=FALSE,Rowv=FALSE                          # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR1","","Mutated","Wild Type"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","blue","white","#8b0000","grey"),lty= 1,lwd = 5,cex = 0.3)
#dev.off()
dim(RNASeq.subset.UP)

#heatmap DOWN
my.palette <- colorRampPalette(c("red", "yellow", "blue"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))
dev.new(width=10, height=10)
#png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.TCGA.",filter.label,"SIGN.DOWN.rankmixorder.ICR1vs4.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.3(RNASeq.subset.DOWN,
          main = paste0("SIGN.DOWN.MUTvsWT.in.",Cancerset),
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
          #margins=c(2,8),
          Colv=FALSE,Rowv=FALSE# reorder row/columns by dendogram
        
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR1","","Mutated","Wild Type"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","blue","white","#8b0000","grey"),lty= 1,lwd = 5,cex = 0.3)
#dev.off()
dim(RNASeq.subset.DOWN)



