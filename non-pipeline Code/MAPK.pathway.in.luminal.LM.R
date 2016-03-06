
# Setup environment
#rm(list=ls())
## dependencies
## install java for xlsx export
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("xlsx","plyr","ggplot2","reshape","plotrix")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (plyr)
library (reshape)
library (ggplot2)
library(mygene)
library(plotrix)

setwd("~/Dropbox/BREAST_QATAR/")

# Load expression data
load ("./2 DATA/LM.BRCA/LM.Dataset.split.Rdata")
Expression.Data<-Expression.Data[,-1955]
#MAPK.genes <- read.csv ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/genes.from.MAPK.pathway.csv") FROM our pathway analysis
#colnames(MAPK.genes) <-"entrez" 
#gene.names <- queryMany(MAPK.genes$entrez, scopes="entrezgene", fields=c("symbol", "go"), species="human")
#gene.table <- as.data.frame(gene.names)
#MAPK.genes$symbol <- gene.table$symbol[match(MAPK.genes$entrez,gene.table$query)]
#MAPK.genes <- unique(MAPK.genes)

#Genes selection
#MAPK.genes <- read.csv ("./3 ANALISYS/MAPK/MAPK.pathway.genes.csv",header = FALSE,stringsAsFactors = FALSE)
#MAPK.genes <- MAPK.genes[,2,drop=FALSE]
#colnames(MAPK.genes) <- "symbol"

#MAPK.genes <- DEGinICR4vs1.MAPK.PW
#MAPK.genes$symbol <- rownames(MAPK.genes)

#LM UPREGULATED

MAPK.genes <- DEGinMUTvsWT.lum.MAPK.PW.UP
MAPK.genes$symbol <- rownames(MAPK.genes)

Available.probes <- data.frame(Probe.ID=Gene.Meta.data$Affy_Probe_ID[Gene.Meta.data$Symbol %in% MAPK.genes$symbol],Symbol=Gene.Meta.data$Symbol[Gene.Meta.data$Symbol %in% MAPK.genes$symbol])
rownames(Available.probes) <- Available.probes$Probe.ID
MAPk.expresion.matrix <- Expression.Data[rownames(Available.probes),]


#load cluster data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/LM.Dataset/LM.Dataset.MA.k7.DBGS3.reps5000/LM.Dataset.MA.k7.DBGS3.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data

#select data
heatmap.table <-as.data.frame(t(MAPk.expresion.matrix))
heatmap.table$subtype <- Sample.Meta.data$PAM50_SUBTYPE[match(rownames(heatmap.table),rownames(Sample.Meta.data))]
heatmap.table$cluster <- Consensus.class$Group[match(rownames(heatmap.table),Consensus.class$PatientID)]
heatmap.table <- heatmap.table[complete.cases(heatmap.table),]

#Filter
heatmap.table <- heatmap.table[heatmap.table$subtype == "LumA" | heatmap.table$subtype == "LumB",]
#heatmap.table <- heatmap.table[heatmap.table$subtype == "LumB",]
heatmap.table <- heatmap.table[heatmap.table$cluster == "ICR1" | heatmap.table$cluster == "ICR4",]


#split data and create color mapping
heatmap.matrix <- as.matrix(heatmap.table[,1:(ncol(heatmap.table)-2)])
mode(heatmap.matrix) <- "numeric"
heatmap.meta <- heatmap.table[,(ncol(heatmap.table)-1):ncol(heatmap.table)]
heatmap.matrix <- heatmap.matrix[row.names(heatmap.meta),]
means.UP <- as.data.frame(rowMeans(heatmap.matrix))
colnames(means.UP)<- "means.UP"
heatmap.meta <- cbind(heatmap.meta,means.UP$means.UP)
heatmap.meta <- heatmap.meta[order(heatmap.meta$means),]
heatmap.genes <- Available.probes[colnames(heatmap.matrix),]

#combined z-score
heatmap.meta <- heatmap.meta[rownames(means.BOTH),]

heatmap.matrix <- heatmap.matrix[row.names(heatmap.meta),]

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

MAPK.color.scale <-  color.scale(heatmap.meta$means,cs1 = c(0,1),cs2=c(0,0),cs3=c(0,0),alpha=1,extremes=NA,na.color=NA)
MAPK.color.scale <-  color.scale(means.BOTH$BOTH.rank,cs1 = c(0,1),cs2=c(0,0),cs3=c(0,0),alpha=1,extremes=NA,na.color=NA)

heatmap.genes$label <- paste0(heatmap.genes$Symbol," - ",heatmap.genes$Probe.ID)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

meta.matrix <- as.matrix(rbind(cluster.colors,subtype.colors,MAPK.color.scale))
meta.matrix<- t(meta.matrix)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors <- c(seq(-10,-1,length=100),seq(-1,1,length=100),seq(1,10,length=100))
#dev.new(width=6, height=6)
png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.LM.LUM.SIGN.UP.rankmixorder.ICR1vs4.png"),res=600,height=10,width=14,unit="in")
heatmap.3((t(heatmap.matrix)),
          main = "SIGN.UP.MUTvsWT",
          ColSideColors = meta.matrix,
          col=my.palette, 
          breaks=my.colors, 
          trace = "none",
          scale ="row",
          labCol = FALSE,
          labRow = heatmap.genes$label,
          margins=c(2,10),
          Colv = FALSE
          )
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR1"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","blue"),lty= 1,lwd = 5,cex = 0.6)
dev.off()
#correlation between ICR score and mapscore
immunescore<-read.csv("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BRCA.LMDATA.csv",row.names = 1)
immunescore<-immunescore[rownames(heatmap.meta),,drop=FALSE]
cor.test(immunescore$unscaled.score,heatmap.meta$means)

#LM DOWNREGULATED

MAPK.genes <- DEGinMUTvsWT.lum.MAPK.PW.DOWN
MAPK.genes$symbol <- rownames(MAPK.genes)

Available.probes <- data.frame(Probe.ID=Gene.Meta.data$Affy_Probe_ID[Gene.Meta.data$Symbol %in% MAPK.genes$symbol],Symbol=Gene.Meta.data$Symbol[Gene.Meta.data$Symbol %in% MAPK.genes$symbol])
rownames(Available.probes) <- Available.probes$Probe.ID
MAPk.expresion.matrix <- Expression.Data[rownames(Available.probes),]

#reverse expression
MAPk.expresion.matrix <- as.matrix(MAPk.expresion.matrix)
mode (MAPk.expresion.matrix) <- "numeric"
MAPk.expresion.matrix <- -(MAPk.expresion.matrix)

#load cluster data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/LM.Dataset/LM.Dataset.MA.k7.DBGS3.reps5000/LM.Dataset.MA.k7.DBGS3.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data

#select data
heatmap.table <-as.data.frame(t(MAPk.expresion.matrix))
heatmap.table$subtype <- Sample.Meta.data$PAM50_SUBTYPE[match(rownames(heatmap.table),rownames(Sample.Meta.data))]
heatmap.table$cluster <- Consensus.class$Group[match(rownames(heatmap.table),Consensus.class$PatientID)]
heatmap.table <- heatmap.table[complete.cases(heatmap.table),]

#Filter
heatmap.table <- heatmap.table[heatmap.table$subtype == "LumA" | heatmap.table$subtype == "LumB",]
#heatmap.table <- heatmap.table[heatmap.table$subtype == "LumB",]
heatmap.table <- heatmap.table[heatmap.table$cluster == "ICR1" | heatmap.table$cluster == "ICR4",]


#split data and create color mapping
heatmap.matrix <- as.matrix(heatmap.table[,1:(ncol(heatmap.table)-2)])
mode(heatmap.matrix) <- "numeric"
heatmap.meta <- heatmap.table[,(ncol(heatmap.table)-1):ncol(heatmap.table)]
heatmap.matrix <- heatmap.matrix[row.names(heatmap.meta),]
means.DOWN <- as.data.frame(rowMeans(heatmap.matrix))
colnames(means.DOWN)<- "means.DOWN"
heatmap.meta <- cbind(heatmap.meta,means.DOWN$means.DOWN)
heatmap.meta <- heatmap.meta[order(heatmap.meta$means),]
heatmap.genes <- Available.probes[colnames(heatmap.matrix),]

#combined z-score
heatmap.meta <- heatmap.meta[rownames(means.BOTH),]

heatmap.matrix <- heatmap.matrix[row.names(heatmap.meta),]

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

MAPK.color.scale <-  color.scale(heatmap.meta$means,cs1 = c(0,1),cs2=c(0,0),cs3=c(0,0),alpha=1,extremes=NA,na.color=NA)
MAPK.color.scale <-  color.scale(means.BOTH$BOTH.rank,cs1 = c(0,1),cs2=c(0,0),cs3=c(0,0),alpha=1,extremes=NA,na.color=NA)

heatmap.genes$label <- paste0(heatmap.genes$Symbol," - ",heatmap.genes$Probe.ID)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

meta.matrix <- as.matrix(rbind(cluster.colors,subtype.colors,MAPK.color.scale))
meta.matrix<- t(meta.matrix)
my.palette <- colorRampPalette(c("red", "yellow", "blue"))(n = 299)
my.colors <- c(seq(-10,-1,length=100),seq(-1,1,length=100),seq(1,10,length=100))
#dev.new(width=6, height=6)
png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.LM.LUM.SIGN.DOWN.rankmixorder.ICR1vs4.png"),res=600,height=10,width=14,unit="in")
heatmap.3((t(heatmap.matrix)),
          main = "SIGN.DOWN.rev.MUTvsWT",
          ColSideColors = meta.matrix,
          col=my.palette, 
          breaks=my.colors, 
          trace = "none",
          scale ="row",
          labCol = FALSE,
          labRow = heatmap.genes$label,
          margins=c(2,10),
          Colv = FALSE
)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR1"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","blue"),lty= 1,lwd = 5,cex = 0.6)
dev.off()

#combined z-score
means.BOTH <- means.UP
means.BOTH$DOWN <- means.DOWN$means.DOWN[match(rownames(means.BOTH),rownames(means.DOWN))]
colnames(means.BOTH) <- c("UP","DOWN")
#means.BOTH$DOWN <- -(means.BOTH$DOWN)
means.BOTH$UP.rank<-rank(means.BOTH$UP)
means.BOTH$DOWN.rank<-rank(means.BOTH$DOWN)
means.BOTH$BOTH <- rowMeans(means.BOTH[,c(1,2)])
means.BOTH$BOTH.rank <- rowMeans(means.BOTH[,c(3,4)])
means.BOTH <- means.BOTH[order(means.BOTH$BOTH.rank),]
