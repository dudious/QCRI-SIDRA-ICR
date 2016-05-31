# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
## dependencies
required.packages <- c("plotrix","data.table")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(plotrix)
library (data.table)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

#parameters
filter.samples = "Basal-Like" #"Luminal" OR "Basal-Like" OR "HER2-enriched"
genes.probes = "genes"

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

# Load expression data
load ("./2 DATA/LM.BRCA/LM.Dataset.split.Rdata")
Expression.Data<-Expression.Data[,-1955]

#LM UPREGULATED
LM.MAPK.genes.UP <- DEGinMUTvsWT.lum.MAPK.PW.UP
LM.MAPK.genes.UP$symbol <- rownames(LM.MAPK.genes.UP)
Available.probes.UP <- data.frame(Probe.ID=Gene.Meta.data$Affy_Probe_ID[Gene.Meta.data$Symbol %in% LM.MAPK.genes.UP$symbol],Symbol=Gene.Meta.data$Symbol[Gene.Meta.data$Symbol %in% LM.MAPK.genes.UP$symbol])
rownames(Available.probes.UP) <- Available.probes.UP$Probe.ID
MAPk.expresion.matrix.UP <- Expression.Data[rownames(Available.probes.UP),]

#LM DOWNREGULATED
LM.MAPK.genes.DOWN <- DEGinMUTvsWT.lum.MAPK.PW.DOWN
LM.MAPK.genes.DOWN$symbol <- rownames(LM.MAPK.genes.DOWN)
Available.probes.DOWN <- data.frame(Probe.ID=Gene.Meta.data$Affy_Probe_ID[Gene.Meta.data$Symbol %in% LM.MAPK.genes.DOWN$symbol],Symbol=Gene.Meta.data$Symbol[Gene.Meta.data$Symbol %in% LM.MAPK.genes.DOWN$symbol])
rownames(Available.probes.DOWN) <- Available.probes.DOWN$Probe.ID
MAPk.expresion.matrix.DOWN <- Expression.Data[rownames(Available.probes.DOWN),]

#colapse probes to genes
MAPk.expresion.matrix.UP.bygene <- as.matrix(MAPk.expresion.matrix.UP)
mode (MAPk.expresion.matrix.UP.bygene) <- "numeric"
MAPk.expresion.matrix.UP.bygene <- cbind.data.frame(MAPk.expresion.matrix.UP.bygene,Symbol=as.factor(Available.probes.UP[rownames(MAPk.expresion.matrix.UP.bygene),"Symbol"]))
MAPk.expresion.matrix.UP.bygene<-aggregate (.~Symbol, data=MAPk.expresion.matrix.UP.bygene,FUN=mean)
rownames(MAPk.expresion.matrix.UP.bygene) <- MAPk.expresion.matrix.UP.bygene$Symbol
MAPk.expresion.matrix.UP.bygene$Symbol <-NULL
MAPk.expresion.matrix.DOWN.bygene <- as.matrix(MAPk.expresion.matrix.DOWN)
mode (MAPk.expresion.matrix.DOWN.bygene) <- "numeric"
MAPk.expresion.matrix.DOWN.bygene <- cbind.data.frame(MAPk.expresion.matrix.DOWN.bygene,Symbol=as.factor(Available.probes.DOWN[rownames(MAPk.expresion.matrix.DOWN.bygene),"Symbol"]))
MAPk.expresion.matrix.DOWN.bygene<-aggregate (.~Symbol, data=MAPk.expresion.matrix.DOWN.bygene,mean)
rownames(MAPk.expresion.matrix.DOWN.bygene) <- MAPk.expresion.matrix.DOWN.bygene$Symbol
MAPk.expresion.matrix.DOWN.bygene$Symbol <-NULL

#load cluster data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/LM.Dataset/LM.Dataset.MA.k7.DBGS3.reps5000/LM.Dataset.MA.k7.DBGS3.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
rownames(Consensus.class) <-  Consensus.class$PatientID
Consensus.class$X <-NULL

#select data
heatmap.table <-data.frame(t(MAPk.expresion.matrix.UP.bygene))
heatmap.table$subtype <- Sample.Meta.data$PAM50_SUBTYPE[match(rownames(heatmap.table),rownames(Sample.Meta.data))]
heatmap.table$cluster <- Consensus.class$Group[match(rownames(heatmap.table),rownames(Consensus.class))]
heatmap.table <- heatmap.table[complete.cases(heatmap.table),]
#drop normal-like & Claudine low
heatmap.table<-heatmap.table[-which(heatmap.table$subtype=="Normal"),]
heatmap.table<-heatmap.table[-which(heatmap.table$subtype=="ClaudinLow"),]
table(heatmap.table$subtype)
#subtype filter
filter.label <- ""
if (filter.samples == "Luminal") {
  heatmap.table <- heatmap.table[heatmap.table$subtype == "LumA" | heatmap.table$subtype == "LumB",]
  filter.label <-"LUM."
}
if (filter.samples == "Basal-Like") {
  heatmap.table <- heatmap.table[heatmap.table$subtype == "Basal",]
  filter.label <- "BASAL."
}
if (filter.samples == "HER2-enriched") {
  heatmap.table <- heatmap.table[heatmap.table$subtype == "Her2",]
  filter.label <- "HER2."
}
#heatmap.table <- heatmap.table[heatmap.table$subtype == "LumB",]
heatmap.table <- heatmap.table[heatmap.table$cluster == "ICR1" | heatmap.table$cluster == "ICR4",]

#Combined Z-score
heatmap.matrix.UP <- as.matrix(heatmap.table[,1:(ncol(heatmap.table)-2)])
mode(heatmap.matrix.UP) <- "numeric"
heatmap.matrix.DOWN <- as.matrix(t(MAPk.expresion.matrix.DOWN.bygene[,rownames(heatmap.matrix.UP)]))
mode(heatmap.matrix.DOWN) <- "numeric"
heatmap.matrix.DOWN <- -(heatmap.matrix.DOWN)
means.UP <- as.data.frame(rowMeans(heatmap.matrix.UP))
means.DOWN <- as.data.frame(rowMeans(heatmap.matrix.DOWN))
colnames(means.UP)<- "UP"
colnames(means.DOWN)<- "DOWN"
means.BOTH <- means.UP
means.BOTH$DOWN <- means.DOWN$DOWN[match(rownames(means.BOTH),rownames(means.DOWN))]
colnames(means.BOTH) <- c("UP","DOWN")
means.BOTH$UP.rank<-rank(means.BOTH$UP)
means.BOTH$DOWN.rank<-rank(means.BOTH$DOWN)
means.BOTH$BOTH <- rowMeans(means.BOTH[,c(1,2)])
means.BOTH$BOTH.rank <- rowMeans(means.BOTH[,c(3,4)])
means.BOTH <- means.BOTH[order(means.BOTH$BOTH.rank),]

#reorder columns (samples)
heatmap.meta <- heatmap.table[,(ncol(heatmap.table)-1):ncol(heatmap.table)]
heatmap.meta$Z.score <- means.BOTH$BOTH.rank[match(rownames(heatmap.meta),rownames(means.BOTH))]
heatmap.meta <- heatmap.meta[order(heatmap.meta$Z.score),]
#heatmap.meta <- heatmap.meta[sample.order,] #order like other script
heatmap.matrix.UP <- heatmap.matrix.UP[rownames(heatmap.meta),]
heatmap.matrix.DOWN <- heatmap.matrix.DOWN[rownames(heatmap.meta),]


#creat color legend
cluster.colors <- heatmap.meta$cluster
levels (cluster.colors) <- c(levels (cluster.colors),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
cluster.colors[cluster.colors=="ICR4"] <- "#FF0000"
cluster.colors[cluster.colors=="ICR3"] <- "#FFA500"
cluster.colors[cluster.colors=="ICR2"] <- "#00FF00"
cluster.colors[cluster.colors=="ICR1"] <- "#0000FF"
cluster.colors <- as.character(cluster.colors)

subtype.colors <-heatmap.meta$subtype
levels (subtype.colors) <- c(levels (subtype.colors),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","#000000"))
subtype.colors[subtype.colors=="LumA"] <- "#eaff00"
subtype.colors[subtype.colors=="LumB"] <- "#00c0ff"
subtype.colors[subtype.colors=="Basal"] <- "#da70d6"
subtype.colors[subtype.colors=="Her2"] <- "#daa520"
subtype.colors[subtype.colors=="Normal"] <- "#d3d3d3"
subtype.colors[subtype.colors=="ClaudinLow"] <- "#000000"


subtype.colors <- as.character(subtype.colors)

MAPK.color.scale <-  color.scale(heatmap.meta$Z.score,alpha=1,extremes=c("#ffc1cb","#ff284d"),na.color=NA)

meta.matrix <- as.matrix(rbind(cluster.colors,subtype.colors,MAPK.color.scale))
meta.matrix<- t(meta.matrix)

#reorder rows (genes)
up.genes.order <- c("TAOK2","TP53","MAPK3","MAP3K1","MAPT","HSPA1A","FLNB","TAOK3","CRK","RPS6KA2","MAP2K4","DUSP5","CACNA1D","MAPK8","RASGRP1","CACNA1G")
down.genes.order <- c("CACNG6","CACNA1B","CACNA2D3","FASLG","RASGRF1","JUN","JUND","DUSP16","PPM1B","SOS1","FGF12","RASGRP2","PRKCB","MAP4K1","PTPN7","GADD45G","DDIT3","DUSP8","DUSP10","FGFR4","FGF14","FGF13","MAP2K6","DUSP2")
heatmap.matrix.UP <- heatmap.matrix.UP[,up.genes.order]
heatmap.matrix.DOWN <- heatmap.matrix.DOWN[,down.genes.order[which(down.genes.order %in% colnames(heatmap.matrix.DOWN))]]

#LM UPREGULATED heatmap
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors <- c(seq(-10,-1,length=100),seq(-1,1,length=100),seq(1,10,length=100))
#dev.new(width=6, height=6)
png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.LM.",filter.label,"SIGN.UP.rankmixorder.ICR1vs4.bygene.png"),res=600,height=10,width=14,unit="in")
heatmap.3((t(heatmap.matrix.UP)),
          main = "SIGN.UP.MUTvsWT",
          ColSideColors = meta.matrix,
          col=my.palette, 
          breaks=my.colors, 
          trace = "none",
          scale ="row",
          labCol = FALSE,
          margins=c(2,10),
          Colv = FALSE,Rowv = FALSE,
          cexRow=1,cexCol=0.1
          )
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR1"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","blue"),lty= 1,lwd = 5,cex = 0.6)
dev.off()


#LM DOWNREGULATED heatmap
my.palette <- colorRampPalette(c("red", "yellow", "blue"))(n = 299)
my.colors <- c(seq(-10,-1,length=100),seq(-1,1,length=100),seq(1,10,length=100))
#dev.new(width=6, height=6)
png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.LM.",filter.label,"SIGN.DOWN.rankmixorder.ICR1vs4.bygene.png"),res=600,height=10,width=14,unit="in")
heatmap.3((t(heatmap.matrix.DOWN)),
          main = "SIGN.DOWN.rev.MUTvsWT",
          ColSideColors = meta.matrix,
          col=my.palette, 
          breaks=my.colors, 
          trace = "none",
          scale ="row",
          labCol = FALSE,
          margins=c(2,10),
          Colv = FALSE,Rowv = FALSE,
          cexRow=1,cexCol=0.1
)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR1"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","blue"),lty= 1,lwd = 5,cex = 0.6)
dev.off()



