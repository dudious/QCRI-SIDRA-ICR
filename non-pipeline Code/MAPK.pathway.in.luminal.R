
# Setup environment
rm(list=ls())
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
library(limma)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")

setwd("~/Dropbox/BREAST_QATAR/")


# Load expression data
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
#MAPK.genes <- read.csv ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/genes.from.MAPK.pathway.csv") FROM our pathway analysis
#colnames(MAPK.genes) <-"entrez" 
#gene.names <- queryMany(MAPK.genes$entrez, scopes="entrezgene", fields=c("symbol", "go"), species="human")
#gene.table <- as.data.frame(gene.names)
#MAPK.genes$symbol <- gene.table$symbol[match(MAPK.genes$entrez,gene.table$query)]
#MAPK.genes <- unique(MAPK.genes)
MAPK.genes <- read.csv ("./3 ANALISYS/MAPK/MAPK.pathway.genes.csv",header = FALSE,stringsAsFactors = FALSE)
MAPK.genes <- MAPK.genes[,2,drop=FALSE]
colnames(MAPK.genes) <- "symbol"

Available.genes <- MAPK.genes$symbol[MAPK.genes$symbol %in% rownames(RNASeq.NORM_Log2)]
MAPk.expresion.matrix <- RNASeq.NORM_Log2[Available.genes,]
rm(RNASeq.NORM_Log2)

#load clinical
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.BSF2.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL

#load MAPKMUT data
MAPMUT <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")

#load cluster data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data

#load immunoscore
immunescores <- read.csv("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BRCA.DBGS3.csv",row.names = 1)

#select data
heatmap.table <-as.data.frame(t(MAPk.expresion.matrix))
heatmap.table$subtype <- Clinical.data$TCGA.PAM50.RMethod.RNASeq[match(rownames(heatmap.table),rownames(Clinical.data))]
heatmap.table$cluster <- Consensus.class$Group[match(rownames(heatmap.table),Consensus.class$PatientID)]
heatmap.table$MAPKs <- MAPMUT$MAP2K4.MAP3K1[match(rownames(heatmap.table),MAPMUT$Sample)]
heatmap.table$IS <- immunescores$scaled.IS[match(rownames(heatmap.table),rownames(immunescores))]
heatmap.table <- heatmap.table[complete.cases(heatmap.table),]

#Filter
#heatmap.table <- heatmap.table[heatmap.table$subtype == "Luminal A" | heatmap.table$subtype == "Luminal B",]
#heatmap.table <- heatmap.table[heatmap.table$subtype == "Luminal A" ,]
#heatmap.table <- heatmap.table[heatmap.table$subtype == "Basal-like",]
heatmap.table <- heatmap.table[heatmap.table$cluster == "ICR1" | heatmap.table$cluster == "ICR4",]


#split data and create color mapping
heatmap.matrix <- as.matrix(heatmap.table[,1:(ncol(heatmap.table)-4)])
mode(heatmap.matrix) <- "numeric"
heatmap.meta <- heatmap.table[,(ncol(heatmap.table)-3):ncol(heatmap.table)]
means <- rowMeans(heatmap.matrix)
heatmap.meta <- cbind(heatmap.meta,means)
heatmap.meta <- heatmap.meta[order(heatmap.meta$means),]
heatmap.matrix <- heatmap.matrix[row.names(heatmap.meta),]

cluster.colors <- heatmap.meta$cluster
levels (cluster.colors) <- c(levels (cluster.colors),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
cluster.colors[cluster.colors=="ICR4"] <- "#FF0000"
cluster.colors[cluster.colors=="ICR3"] <- "#FFA500"
cluster.colors[cluster.colors=="ICR2"] <- "#00FF00"
cluster.colors[cluster.colors=="ICR1"] <- "#0000FF"
cluster.colors <- as.character(cluster.colors)

subtype.colors <-heatmap.meta$subtype
levels (subtype.colors) <- c(levels (subtype.colors),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
subtype.colors[subtype.colors=="Luminal A"] <- "#eaff00"
subtype.colors[subtype.colors=="Luminal B"] <- "#00c0ff"
subtype.colors[subtype.colors=="Basal-like"] <- "#da70d6"
subtype.colors[subtype.colors=="HER2-enriched"] <- "#daa520"
subtype.colors[subtype.colors=="Normal-like"] <- "#d3d3d3"
subtype.colors <- as.character(subtype.colors)

mutation.colors <-heatmap.meta$MAPKs
levels (mutation.colors) <- c(levels (mutation.colors),c("#8b0000","grey"))
mutation.colors[mutation.colors=="MUT"] <- "#8b0000"
mutation.colors[mutation.colors=="WT"] <- "grey"
mutation.colors <- as.character(mutation.colors)

MAPK.color.scale <-  color.scale(heatmap.meta$means,cs1 = c(0,1),cs2=c(0,0),cs3=c(0,0),alpha=1,extremes=NA,na.color=NA)

library(gplots)
#dev.new(width=15, height=15)
#heatmap.2((t(heatmap.matrix)),
#          ColSideColors = mutation.colors,
#          trace = "none",
#          scale ="row",
#          labCol = FALSE,
#          Colv = TRUE)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

meta.matrix <- as.matrix(rbind(cluster.colors,subtype.colors,mutation.colors,MAPK.color.scale))
meta.matrix<- t(meta.matrix)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors <- c(seq(-5,-1,length=100),seq(-1,1,length=100),seq(1,5,length=100))
dev.new(width=6, height=6)
heatmap.3((t(heatmap.matrix)),
          ColSideColors = meta.matrix,
          col=my.palette, 
          breaks=my.colors, 
          trace = "none",
          scale ="row",
          labCol = FALSE,
          Colv = FALSE
          )

#correlation between ICR score and mapscore
cor.test(heatmap.meta$IS,heatmap.meta$means)


                 
