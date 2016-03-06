#################################################################
###
### This Script creates a subset of the LM data 
### for a selected set  genes.(Gene_selection_XXX.txt)
### source data :
### "./2 DATA/SUBSETS/Gene_selection_v2.2.txt"
### ~/Dropbox/Data_LM/LM.Dataset.split.Rdata"
### 
### Results are saved in
### ./2 DATA/SUBSETS/
### File to use :
### 
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

required.packages <- c("clue")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages( "./1 CODE/R tools/clue_0.3-49.zip", repos=NULL,type= "win.binary")

required.packages.BioC <- c("ConsensusClusterPlus")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

library(ConsensusClusterPlus)
library(clue)
source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/stefanofunctions.R")

# Parameters
Cancerset <- "LM.Dataset"
Geneset <- "DBGS3"

# load data
load ("./2 DATA/LM.BRCA/LM.Dataset.split.Rdata")                                
DEGs <- read.csv ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.entrez.csv")
rownames(DEGs) <- DEGs$X
DEGs$X <- NULL
MAPK.genes <- data.frame(Entrez=c(5599, 776, 4137, 4214, 3303, 9344, 2317, 7157, 1398, 6196, 8913, 6416, 51347, 5595, 1847, 10125)) #16 from DEG in Luminal MUTATED
MAPK.genes.clusterDEG <- data.frame(Entrez=c(783, 4776, 786, 27330, 59285, 627, 4908, 9479, 1846, 55799, 3306, 8399, 776, 4137,     #26 from DEG luminal ICR
                                             27092, 778, 9254, 7786, 2261, 2260, 8912, 2259, 8913, 2255, 2257, 4915))

MAPK.genes$Name <- rownames(DEGs)[match(MAPK.genes$Entrez,DEGs$entrez)]
MAPK.genes.clusterDEG$Name <- rownames(DEGs)[match(MAPK.genes.clusterDEG$Entrez,DEGs$entrez)]


# check availabilety of the genes in the dataset
available.genes.MA <- MAPK.genes$Name[which(MAPK.genes$Name %in% Gene.Meta.data$Symbol)]
unavailable.genes.MA <- MAPK.genes$Name[-which(MAPK.genes$Name %in% Gene.Meta.data$Symbol)]
c("B7-H","B7-H1","B7H1","PD-L1","PDL1") %in% Gene.Meta.data$Symbol                                # check for CD274 alternative names

## Subset data
MA.probes.subset <- Gene.Meta.data[Gene.Meta.data$Symbol %in% available.genes.MA,]
MA.subset <- t(Expression.Data[rownames(Expression.Data) %in% MA.probes.subset$Affy_Probe_ID,])
mode(MA.subset) <- "numeric"

## report
print (paste0("Geneset Selected : MAPK DEG genes in luminal",Geneset))
print ("Genes selected for MA : ")
print (available.genes.MA)
print ("Genes missing for MA :")
print (unavailable.genes.MA)

## filter luminal patients
Luminal.samples <- rownames(Sample.Meta.data[Sample.Meta.data$PAM50_SUBTYPE %in% c("LumA","LumB"),])
MA.subset.luminal <- MA.subset[Luminal.samples,]

## add ICR cluster
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/",Cancerset,".MA.k7.",
                                   Geneset,".reps5000/",Cancerset,".MA.k7.",
                                   Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
rownames(Consensus.class) <- Consensus.class$PatientID
MA.subset.luminal <- as.data.frame(MA.subset.luminal)
MA.subset.luminal$Cluster <- Consensus.class$Group[match(rownames(MA.subset.luminal),Consensus.class$PatientID)]
MA.subset.luminal <- MA.subset.luminal[MA.subset.luminal$Cluster %in% c("ICR1","ICR4"),]

## Stats ICR4 vs ICR1 Of these genes in LM data
stats <- data.frame(probes = colnames(MA.subset.luminal[,-ncol(MA.subset.luminal)]),p.value=NA)
for (i in 1:(ncol(MA.subset.luminal)-1)) {
  test <- kruskal.test(MA.subset.luminal[,i]~MA.subset.luminal$Cluster)
  stats$mean.ICR1[i] <- mean(MA.subset.luminal[MA.subset.luminal$Cluster == "ICR1",i])
  stats$mean.ICR4[i] <- mean(MA.subset.luminal[MA.subset.luminal$Cluster == "ICR4",i])
  stats$logFC[i] <- log(stats$mean.ICR4[i]/stats$mean.ICR1[i])
  stats$p.value[i] <- test$p.value 
}

stats$gene <- Gene.Meta.data$Symbol[match(stats$probes,Gene.Meta.data$Affy_Probe_ID)]
stats <- stats[order(stats$p.value),]
stats$p.value.adjusted <- p.adjust(stats$p.value, method = "bonferroni", n = length(stats$p.value))
stats <- stats[,c("gene","probes","p.value","p.value.adjusted","mean.ICR4","mean.ICR1","logFC")]

write.csv(stats,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/stats_LM_datset_MAPK_pathway.csv")

DEGs.filtered <- DEGs[available.genes.MA,]
DEGs.filtered <- DEGs.filtered[order(DEGs.filtered$FDR),]

write.csv(DEGs.filtered,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/stats_TCGA_MAPK_pathway.csv")

## TCGA original data
Cancerset <- "BRCA"
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
#MAPK pathway
RNASeq.NORM_Log2.filtered <- as.data.frame(t(RNASeq.NORM_Log2[available.genes.MA,]))
#Cluster & Subtype
Cancerset <- "BRCA.BSF2"
Geneset <- "DBGS3.FLTR"
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
rownames(Consensus.class) <- Consensus.class$PatientID
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL

RNASeq.NORM_Log2.filtered$Cluster <- Consensus.class$Group[match(rownames(RNASeq.NORM_Log2.filtered),Consensus.class$PatientID)] 
RNASeq.NORM_Log2.filtered$Subtype <- Clinical.data$TCGA.PAM50.RMethod.RNASeq[match(rownames(RNASeq.NORM_Log2.filtered),rownames(Clinical.data))]
RNASeq.NORM_Log2.filtered <- RNASeq.NORM_Log2.filtered[RNASeq.NORM_Log2.filtered$Cluster %in% c("ICR1","ICR4"),]
RNASeq.NORM_Log2.filtered <- RNASeq.NORM_Log2.filtered[RNASeq.NORM_Log2.filtered$Subtype %in% c("Luminal A","Lumninal B"),]

## Stats
stats <- data.frame(genes = colnames(RNASeq.NORM_Log2.filtered[,-c((ncol(RNASeq.NORM_Log2.filtered)-1),ncol(RNASeq.NORM_Log2.filtered))]),p.value=NA)
for (i in 1:(ncol(RNASeq.NORM_Log2.filtered)-2)) {
  test <- kruskal.test(RNASeq.NORM_Log2.filtered[,i]~RNASeq.NORM_Log2.filtered$Cluster)
  stats$mean.ICR1[i] <- mean(RNASeq.NORM_Log2.filtered[RNASeq.NORM_Log2.filtered$Cluster == "ICR1",i])
  stats$mean.ICR4[i] <- mean(RNASeq.NORM_Log2.filtered[RNASeq.NORM_Log2.filtered$Cluster == "ICR4",i])
  stats$logFC[i] <- log(stats$mean.ICR4[i]/stats$mean.ICR1[i])
  stats$p.value[i] <- test$p.value 
}

stats <- stats[order(stats$p.value),]
stats$p.value.adjusted <- p.adjust(stats$p.value, method = "bonferroni", n = length(stats$p.value))
stats <- stats[,c("genes","mean.ICR4","mean.ICR1","logFC","p.value","p.value.adjusted")]

write.csv(stats,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/statsNEW_TCGA_MAPK_pathway.csv")


#heatmaps
library("gplots")
source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

heatmap.genes <- MAPK.genes[MAPK.genes$Name %in% c("CACNA1D","FLNB","MAPT","MAP2K4","MAPK8","HSPA1A","TAOK2"),] #adjusted p value <0.01 in TCGA & LM dataset (ICR4vs1,Luminal)
heatmap.probes <- c("210108_at","207998_s_at","203266_s_at","208614_s_at","203929_s_at","204986_s_at","210671_x_at","200799_at","200800_s_at")
TCGA.heatmap <- RNASeq.NORM_Log2.filtered[,heatmap.genes$Name]
TCGA.heatmap.meta <- RNASeq.NORM_Log2.filtered[,"Cluster",drop=FALSE]
LM.heatmap <- MA.subset.luminal[,heatmap.probes]
LM.heatmap.meta <- MA.subset.luminal[,"Cluster",drop=FALSE]
MAPKS.state <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
TCGA.heatmap.meta$MAPKs <- MAPKS.state$MAP2K4.MAP3K1[match(rownames(TCGA.heatmap.meta),MAPKS.state$Sample)]
TCGA.heatmap.meta <- TCGA.heatmap.meta[order(TCGA.heatmap.meta$Cluster,TCGA.heatmap.meta$MAPKs),]
LM.heatmap.meta <- LM.heatmap.meta[order(LM.heatmap.meta$Cluster),,drop=FALSE]

#TCGA

color.matrix <- as.matrix(TCGA.heatmap.meta[which(complete.cases(TCGA.heatmap.meta)),])
TCGA.heatmap <- TCGA.heatmap[rownames(color.matrix),]
color.matrix[color.matrix[,"Cluster"]=="ICR4","Cluster"] <- "#FF0000"
color.matrix[color.matrix[,"Cluster"]=="ICR1","Cluster"] <- "#0000FF"
color.matrix[color.matrix[,"MAPKs"]=="MUT","MAPKs"] <-"#8b0000"
color.matrix[color.matrix[,"MAPKs"]=="WT","MAPKs"] <- "#ffff00"
color.matrix <- t(color.matrix)
colnames(color.matrix) <-NULL
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
png(paste0("./4 FIGURES/Heatmaps/MAPK.pathway.TCGA.Tree.png"),res=600,height=9,width=18,unit="in")     # set filename
heatmap.3(TCGA.heatmap,
          main = "HeatMap : RNASeq expression - MAPK.pathway - TCGA/Luminal - ICR1 vs ICR4",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=color.matrix,                         # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="col", 
          density.info="none",
          trace="none",
          #labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2,
          margins=c(10,2),
          labRow=FALSE,
          Colv=TRUE, Rowv=TRUE                              # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topleft",legend = c("ICR 1",  "ICR 4"  ,""     ,"Mutated","Wild-type"),
                    col = c("#0000FF","#FF0000","white","#8b0000","#ffff00"),lty= 1,lwd = 5,cex = 1.3)
dev.off()

#LM

color.matrix <- as.matrix(LM.heatmap.meta[which(complete.cases(LM.heatmap.meta)),])
LM.heatmap <- LM.heatmap[rownames(LM.heatmap.meta),]
color.matrix[color.matrix[,1]=="ICR4",1] <- "#FF0000"
color.matrix[color.matrix[,1]=="ICR1",1] <- "#0000FF"
color.matrix <- t(color.matrix)
colnames(color.matrix) <-NULL
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
png(paste0("./4 FIGURES/Heatmaps/MAPK.pathway.LM.png"),res=600,height=9,width=18,unit="in")     # set filename
heatmap.3(LM.heatmap,
          main = "HeatMap : MA expression - MAPK.pathway - LM/Luminal - ICR1 vs ICR4",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=color.matrix,                         # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="col", 
          density.info="none",
          trace="none",
          #labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2,
          margins=c(15,5),
          labRow=FALSE,
          Colv=TRUE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topleft",legend = c("ICR 1",  "ICR 4"  ,""     ,"Mutated","Wild-type"),
       col = c("#0000FF","#FF0000","white","#8b0000","#ffff00"),lty= 1,lwd = 5,cex = 1.3)
dev.off()

