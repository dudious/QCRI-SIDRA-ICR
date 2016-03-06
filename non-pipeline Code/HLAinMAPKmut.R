
# Setup environment
rm(list=ls())
## dependencies
## install java for xlsx export
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("xlsx","plyr","ggplot2","reshape","survival")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (xlsx) #xlsx needs java installed
library (plyr)
library(survival)
library (reshape)
library (ggplot2)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/ggkm.R")

setwd("~/Dropbox/BREAST_QATAR/")


# Load data files
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
HLA.genes <- RNASeq.NORM_Log2[grepl("HLA",rownames(RNASeq.NORM_Log2)),]
HLA.expression.matrix <- as.data.frame(t(RNASeq.NORM_Log2["PTX3",,drop=FALSE]))
rm(RNASeq.NORM_Log2)
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.BSF2.RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL

#load MAPKMUT data
MAPMUT <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")

#load cluster data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data

#select data
heatmap.table <-as.data.frame(t(HLA.genes))
heatmap.table$subtype <- Clinical.data$TCGA.PAM50.RMethod.RNASeq[match(rownames(heatmap.table),rownames(Clinical.data))]
heatmap.table$cluster <- Consensus.class$Group[match(rownames(heatmap.table),Consensus.class$PatientID)]
heatmap.table$MAPKs <- MAPMUT$MAP2K4.MAP3K1[match(rownames(heatmap.table),MAPMUT$Sample)]
heatmap.table <- heatmap.table[complete.cases(heatmap.table),]

#Filter
heatmap.table <- heatmap.table[heatmap.table$subtype == "Luminal A" | heatmap.table$subtype == "Luminal B",]
heatmap.table <- heatmap.table[heatmap.table$cluster == "ICR1" | heatmap.table$cluster == "ICR4",]


#split data and create color mapping
heatmap.matrix <- as.matrix(heatmap.table[,1:(ncol(heatmap.table)-3)])
mode(heatmap.matrix) <- "numeric"
heatmap.meta <- heatmap.table[,(ncol(heatmap.table)-2):ncol(heatmap.table)]
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


library(gplots)
dev.new(width=15, height=15)
heatmap.2((t(heatmap.matrix)),
          ColSideColors = mutation.colors,
          trace = "none",
          scale ="row",
          labCol = FALSE,
          Colv = TRUE)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

meta.matrix <- as.matrix(rbind(cluster.colors,subtype.colors,mutation.colors))
meta.matrix<- t(meta.matrix)
dev.new(width=6, height=6)
heatmap.3((t(heatmap.matrix)),
          ColSideColors = meta.matrix,
          trace = "none",
          scale ="row",
          labCol = FALSE,
          Colv = FALSE)

