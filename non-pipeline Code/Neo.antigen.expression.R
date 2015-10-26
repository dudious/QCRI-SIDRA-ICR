# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

#install missing packages
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("edgeR")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

#Load the packages
library(gplots)
library(edgeR)

#load the functions

TimeUse <-function(func){ 
  time1 <- proc.time() 
  result <-func
  time2 <- proc.time()
  show(timeelapsed <- time2-time1)
}

edgeRDEA4GSEAOff <- function(geneCounts, groups, baseLine){
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  DEA <- DGEList(geneCounts, group = groups)
  DEA <- calcNormFactors(DEA)
  keep <- rowSums(cpm(DEA)>1)>=10
  DEA <- DEA[keep,]
  DEA$samples$lib.size <- colSums(DEA$counts)
  message("Estimating GLM Common Dispersion")
  TimeUse(DEA <- estimateGLMCommonDisp(DEA, design = design))
  message("Estimating GLM Tagwise Dispersion")
  TimeUse(DEA <- estimateGLMTagwiseDisp(DEA, design = design))
  message("Fitting the GLM Model")
  TimeUse(fit <- glmFit(DEA, design = design))
  contrast <- rep(1, 2)
  blId <- which(colnames(design) == baseLine)
  contrast[blId] <- -1
  message("Maximum Likelihood Estimate")
  TimeUse(test <- glmLRT(fit, contrast = contrast))
  ans <- topTags(test, n = nrow(DEA$counts))
  ans <- ans$table
  diffGenes <- abs(ans[, "logFC"]) >= 1 & ans[, "FDR"] <= 0.01
  detags <- rownames(ans)[diffGenes]
  pdf("smearPlot.pdf")
  plotSmear(test,de.tags=detags)
  abline(h=c(-1,1),col="blue")
  dev.off()
  return(ans)
}

##############################################################################



Signif.tresh = 0.05
method = "edgeR"

#load neo-antigen data
neo.antigens <- read.csv("./2 DATA/Rooneys.neo.antigens.csv")
neo.antigens <- neo.antigens [neo.antigens$Cancer=="BRCA",]
neo.antigens.genes <- c(as.character(neo.antigens$Gene),c("MAGEC2","MAGEC3"))

#subset RNAseq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")
neo.antigens.genes.available <- neo.antigens.genes[which(neo.antigens.genes %in% rownames(RNASeq.NORM))]
RNASeq.subset <- as.data.frame(t(RNASeq.NORM[neo.antigens.genes.available,]))

#Add cluster assignment
cluster.assignment <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv")
rownames(cluster.assignment) <-  cluster.assignment$PatientID
cluster.assignment$X <- NULL 
colnames(cluster.assignment) <- c("PatientID","Cluster")
RNASeq.subset$Cluster <- cluster.assignment$Cluster[match(rownames(RNASeq.subset),cluster.assignment$PatientID)]
RNASeq.subset <- RNASeq.subset[-which(is.na(RNASeq.subset$Cluster)),]

#ordering of the clusters
RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
RNASeq.subset$Cluster <- NULL

#re-order the labels
cluster.assignment <- cluster.assignment[rownames(RNASeq.subset),]

# Heatmap 2 (simple no extra annotations)
RNASeq.subset <- as.matrix(RNASeq.subset)
mode(RNASeq.subset) <- "numeric"
patientcolors <- cluster.assignment
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
#patientcolors$Cluster <- droplevels(patientcolors$Cluster)
patientcolors <- as.character(patientcolors$Cluster)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-2,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))
png(paste0("./4 FIGURES/Heatmaps/Heatmap.Neo-antigens.RNASeq.TCGA.BRCA-BSF2.DBGS3.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(RNASeq.subset),
          main = paste0("Heatmap RNASeq - neo-antigens"),
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
          cexRow=0.8,cexCol=0.1,margins=c(2,7),
          Colv=FALSE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

if (method == "kruskal"){

# filter for significant genes
RNASeq.subset <- cbind(RNASeq.subset,cluster.assignment$Cluster)
p.value.matrix <- as.data.frame(matrix(vector(),(ncol(RNASeq.subset)-1),2))
colnames(p.value.matrix) <- c("Gene","p.value")
p.value.matrix$Gene <- colnames(RNASeq.subset)[-length(colnames(RNASeq.subset))]
rownames(p.value.matrix) <- p.value.matrix$Gene
p.value.matrix$Gene <- NULL
for (i in 1:(ncol(RNASeq.subset)-1)) {
  Gene = colnames(RNASeq.subset)[i]
  test = kruskal.test(RNASeq.subset[,i]~RNASeq.subset[,ncol(RNASeq.subset)])
  p.value.matrix[Gene,1] = test$p.value
  
}
signif.genes <- p.value.matrix[p.value.matrix$p.value<Signif.tresh,,drop=FALSE]
RNASeq.subset <- as.data.frame(RNASeq.subset[,rownames(signif.genes)])

#ordering of the clusters
RNASeq.subset$Cluster <- cluster.assignment$Cluster[match(rownames(RNASeq.subset),cluster.assignment$PatientID)]
RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
RNASeq.subset$Cluster <- NULL

#re-order the labels
cluster.assignment <- cluster.assignment[rownames(RNASeq.subset),]

# Heatmap 2 (simple no extra annotations)
RNASeq.subset <- as.matrix(RNASeq.subset)
mode(RNASeq.subset) <- "numeric"
patientcolors <- cluster.assignment
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
#patientcolors$Cluster <- droplevels(patientcolors$Cluster)
patientcolors <- as.character(patientcolors$Cluster)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-2,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))
png(paste0("./4 FIGURES/Heatmaps/Heatmap.Neo-antigens.kruskal.RNASeq.TCGA.BRCA-BSF2.DBGS3.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(RNASeq.subset),
          main = paste0("Heatmap RNASeq - neo-antigens"),
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
          cexRow=0.8,cexCol=0.1,margins=c(2,7),
          Colv=TRUE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

}

##edgeR
if (method == "edgeR"){
  RNASeq.subset <- cbind(RNASeq.subset,cluster.assignment$Cluster)
  groups <- rep("",nrow(RNASeq.subset))
  names(groups)<- rownames(RNASeq.subset)

#ICR4 vs ICR123
G1<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]==4,])
groups[G1]<-"G1"
G2<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]!=4,])
groups[G2]<-"G2"
groups<-as.factor(groups)

RNASeq.subset <- t(RNASeq.subset[,-ncol(RNASeq.subset)])

system.time(DEGs <- edgeRDEA4GSEAOff(RNASeq.subset[,names(groups)], groups = groups,  
                                               baseLine = "G2"))

DEGs<-DEGs[order(DEGs$logFC,decreasing = T),]
write.csv(DEGs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEGinICR4vs123BRCA_BSF2_neoantignes.csv")
signif.genes <- DEGs[DEGs$FDR<Signif.tresh,,drop=FALSE]
save(DEGs,signif.genes,RNASeq.subset,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEGinICR4vs123_neoantignes.RDATA")

#reconstitute heatmap matrix
load(file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEGinICR4vs123_neoantignes.RDATA")
RNASeq.subset <- as.data.frame(RNASeq.subset[rownames(signif.genes),])
RNASeq.subset <- as.data.frame(t(RNASeq.subset))
RNASeq.subset$Cluster <- cluster.assignment$Cluster[match(rownames(RNASeq.subset),cluster.assignment$PatientID)]

#ordering of the clusters
RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
RNASeq.subset$Cluster <- NULL

#re-order the labels
cluster.assignment <- cluster.assignment[rownames(RNASeq.subset),]

# Heatmap 2 (simple no extra annotations)
RNASeq.subset <- as.matrix(RNASeq.subset)
mode(RNASeq.subset) <- "numeric"
patientcolors <- cluster.assignment
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
#patientcolors$Cluster <- droplevels(patientcolors$Cluster)
patientcolors <- as.character(patientcolors$Cluster)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-2,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))
png(paste0("./4 FIGURES/Heatmaps/Heatmap.Neo-antigens.EdgeR.RNASeq.TCGA.BRCA-BSF2.DBGS3.DENDO.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(RNASeq.subset),
          main = paste0("Heatmap RNASeq - neo-antigens"),
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
          cexRow=0.8,cexCol=0.1,margins=c(2,7),
          Colv=TRUE,Rowv=TRUE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

}