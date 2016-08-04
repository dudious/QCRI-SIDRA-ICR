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
filter.samples = "" #"Luminal" OR "Basal-Like" OR "HER2-enriched"
#gene ordere derived from clustering DEG MAPK MUT/WT in luminal
up.genes.order <- c("TAOK2","TP53","MAPK3","MAP3K1","MAPT","HSPA1A","FLNB","TAOK3","CRK","RPS6KA2",
                    "MAP2K4","DUSP5","CACNA1D","MAPK8","RASGRP1","CACNA1G")
down.genes.order <- c("CACNG6","CACNA1B","CACNA2D3","FASLG","RASGRF1","JUN","JUND","DUSP16","PPM1B",
                      "SOS1","FGF12","RASGRP2","PRKCB","MAP4K1","PTPN7","GADD45G","DDIT3","DUSP8",
                      "DUSP10","FGFR4","FGF14","FGF13","MAP2K6","DUSP2")

#load MAPK pathway and MAPKDEG data
MAPK.PW.genes <- read.csv ("./3 ANALISYS/MAPK/MAPK.pathway.genes.csv",header = FALSE)
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.rdata")

#select significant MAPK.PW genes in TCGA ICR DEG  or MAPK DEG
###### UPREGULATED Luminal MUT
DEGinMUTvsWT.lum.MAPK.PW<-DEGs
DEGinMUTvsWT.lum.MAPK.PW<-DEGs[rownames(DEGs)%in%MAPK.PW.genes$V2,] #filter for MAPK genes
DEGinMUTvsWT.lum.MAPK.PW.UP<-DEGinMUTvsWT.lum.MAPK.PW[DEGinMUTvsWT.lum.MAPK.PW$logFC>0,]
DEGinMUTvsWT.lum.MAPK.PW.UP<-DEGinMUTvsWT.lum.MAPK.PW.UP[order(DEGinMUTvsWT.lum.MAPK.PW.UP$FDR),]
DEGinMUTvsWT.lum.MAPK.PW.UP<-DEGinMUTvsWT.lum.MAPK.PW.UP[DEGinMUTvsWT.lum.MAPK.PW.UP$FDR<0.05,]

#######DOWNREGULATED in Luminal MUT
DEGinMUTvsWT.lum.MAPK.PW<-DEGs
DEGinMUTvsWT.lum.MAPK.PW<-DEGs[rownames(DEGs)%in%MAPK.PW.genes$V2,] #filter for MAPK genes
DEGinMUTvsWT.lum.MAPK.PW.DOWN<-DEGinMUTvsWT.lum.MAPK.PW[DEGinMUTvsWT.lum.MAPK.PW$logFC<0,]
DEGinMUTvsWT.lum.MAPK.PW.DOWN<-DEGinMUTvsWT.lum.MAPK.PW.DOWN[order(DEGinMUTvsWT.lum.MAPK.PW.DOWN$FDR),]
DEGinMUTvsWT.lum.MAPK.PW.DOWN<-DEGinMUTvsWT.lum.MAPK.PW.DOWN[DEGinMUTvsWT.lum.MAPK.PW.DOWN$FDR<0.05,]

#RNASeq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
RNASeq.subset.UP <- RNASeq.NORM_Log2[rownames(DEGinMUTvsWT.lum.MAPK.PW.UP),]
RNASeq.subset.DOWN <- RNASeq.NORM_Log2[rownames(DEGinMUTvsWT.lum.MAPK.PW.DOWN),]
#reverse DOWN expression
RNASeq.subset.DOWN <- -(RNASeq.subset.DOWN)

#MAPKMUT data
MAPKMUT <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
rownames(MAPKMUT) <- MAPKMUT$Sample
MAPKMUT$X <- NULL
MAPKMUT$Sample <- NULL

#Load master file
BRCA.master <- read.csv("./3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF2.RNASeq_subset_DBGS3.FLTR.Master.Summary.csv",stringsAsFactors = FALSE)
rownames(BRCA.master) <- BRCA.master$X
BRCA.master$X <- NULL

#clinical data
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
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
#means.BOTH$DOWN <- -(means.BOTH$DOWN)
means.BOTH$UP.rank<-rank(means.BOTH$UP)
means.BOTH$DOWN.rank<-rank(means.BOTH$DOWN)
means.BOTH$BOTH <- rowMeans(means.BOTH[,c(1,2)])
means.BOTH$BOTH.rank <- rowMeans(means.BOTH[,c(3,4)])
means.BOTH <- means.BOTH[order(means.BOTH$BOTH.rank),]

#Re-order according to combined z-score
RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(means.BOTH)]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(means.BOTH)]

#lookup MUTstat
mutstat <- MAPKMUT[colnames(RNASeq.subset.UP),"MAP2K4.MAP3K1",drop=FALSE]
levels (mutstat$MAP2K4.MAP3K1) <- c(levels (mutstat$MAP2K4.MAP3K1),c("#8b0000","grey"))
mutstat<-mutstat[!is.na(mutstat$MAP2K4.MAP3K1),,drop=F]
mutstat$MAP2K4.MAP3K1[mutstat$MAP2K4.MAP3K1=="WT"] <- "grey"
mutstat$MAP2K4.MAP3K1[mutstat$MAP2K4.MAP3K1=="MUT"] <- "#8b0000"
mutstat$TP53 <- BRCA.master$TP53.NS.mutations[match(rownames(mutstat),rownames(BRCA.master))]
mutstat$TP53[mutstat$TP53=="WT"] <- "grey"
mutstat$TP53[mutstat$TP53=="MUT"] <- "#00468b"
#drop data without mutation data
RNASeq.subset.UP <- RNASeq.subset.UP[,which(colnames(RNASeq.subset.UP) %in% rownames(mutstat))]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,which(colnames(RNASeq.subset.DOWN) %in% rownames(mutstat))]

#lookup cluster
cluster <- Consensus.class[colnames(RNASeq.subset.UP),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#lookup subtype
subtype <- ClinicalData.subset[colnames(RNASeq.subset.UP),"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
subtype <- subtype [!is.na(subtype$TCGA.PAM50.RMethod.RNASeq),,drop=F]
levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"]     <- "#eaff00"    
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"]     <- "#00c0ff"     
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"]    <- "#da70d6"   
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#daa520"            
#drop normal-like
subtype<-subtype[-which(subtype$TCGA.PAM50.RMethod.RNASeq=="Normal-like"),,drop=FALSE]
#drop data without subtype data
RNASeq.subset.UP <- RNASeq.subset.UP[,which(colnames(RNASeq.subset.UP) %in% rownames(subtype))]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,which(colnames(RNASeq.subset.DOWN) %in% rownames(subtype))]
cluster <- cluster[which(rownames(cluster) %in% rownames(subtype)),,drop=FALSE]
mutstat <- mutstat[which(rownames(mutstat) %in% rownames(subtype)),,drop=FALSE]

#cluster filter
cluster <- cluster[cluster$Cluster=="#FF0000"|cluster$Cluster=="#0000FF",,drop=FALSE]
subtype <- subtype[rownames(cluster),,drop=FALSE]
mutstat <- mutstat[rownames(cluster),,drop=FALSE]
means.BOTH <- means.BOTH[rownames(cluster),,drop=FALSE]
RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(cluster)]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(cluster)]

#subtype filter
filter.label <-""
if (filter.samples == "Luminal") {
subtype <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="#eaff00"|subtype$TCGA.PAM50.RMethod.RNASeq=="#00c0ff",,drop=FALSE]
cluster <- cluster[rownames(subtype),,drop=FALSE]
mutstat <- mutstat[rownames(subtype),,drop=FALSE]
means.BOTH <- means.BOTH[rownames(subtype),,drop=FALSE]
RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(subtype)]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(subtype)]
filter.label <- "LUM."
}
if (filter.samples == "Basal-Like") {
  subtype <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="#da70d6",,drop=FALSE]
  cluster <- cluster[rownames(subtype),,drop=FALSE]
  mutstat <- mutstat[rownames(subtype),,drop=FALSE]
  means.BOTH <- means.BOTH[rownames(subtype),,drop=FALSE]
  RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(subtype)]
  RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(subtype)]
  filter.label <- "BASAL."
}
if (filter.samples == "HER2-enriched") {
  subtype <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="#daa520",,drop=FALSE]
  cluster <- cluster[rownames(subtype),,drop=FALSE]
  mutstat <- mutstat[rownames(subtype),,drop=FALSE]
  means.BOTH <- means.BOTH[rownames(subtype),,drop=FALSE]
  RNASeq.subset.UP <- RNASeq.subset.UP[,rownames(subtype)]
  RNASeq.subset.DOWN <- RNASeq.subset.DOWN[,rownames(subtype)]
  filter.label <- "HER2."
}
#color matrix
mapmutcolors <- as.character(mutstat$MAP2K4.MAP3K1)
tpmutcolors <- as.character(mutstat$TP53)
clustercolors <- as.character(cluster$Cluster)
subtypecolors <- as.character(subtype$TCGA.PAM50.RMethod.RNASeq)
MAPK.color.scale <-  color.scale(means.BOTH$BOTH.rank,extremes=c("#ffc1cb","#ff284d"),na.color=NA)
color.matrix <- as.matrix(rbind (clustercolors,subtypecolors,mapmutcolors,tpmutcolors,MAPK.color.scale))
rownames(color.matrix) <- c("Cluster","Subtype","MAPK-mutation","TP53-mutation","MAPK.color.scale")

#reorder rows (genes)
RNASeq.subset.UP <- RNASeq.subset.UP[up.genes.order,]
RNASeq.subset.DOWN <- RNASeq.subset.DOWN[down.genes.order,]

#heatmap UP
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))

#dev.new(width=6, height=6)
png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.TCGA.",filter.label,"SIGN.UP.rankmixorder.ICR1vs4.TP53.HC.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.3(RNASeq.subset.UP,
          main = "SIGN.UP.MUTvsWT",
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
dev.off()
dim(RNASeq.subset.UP)

#heatmap DOWN
#SAME as above
#dev.new(width=10, height=10)
png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.TCGA.",filter.label,"SIGN.DOWN.rankmixorder.ICR1vs4.TP53.HC.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.3(RNASeq.subset.DOWN,
          main = "SIGN.DOWN.rev.MUTvsWT",
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
dev.off()
dim(RNASeq.subset.DOWN)



