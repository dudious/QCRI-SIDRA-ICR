#################################################################
###
### This Script PLots Heatmaps based on 
### Consensus Clustering grouping of RNASeq Data
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
#Dependencies
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")

# Load Data
Geneset <- "ISGS" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 4             # SET K here
Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LIHC.TCGA.EDASeq.k4.ISGS.reps2000/LIHC.TCGA.EDASeq.k4.ISGS.reps2000.k=4.consensusClass.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
load (paste0("./2 DATA/SUBSETS/LIHC/TCGA.LIHC.RNASeq.subset.",Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)

#ordeing
RNASeq.subset <- cbind(RNASeq.subset,rowMeans(RNASeq.subset[, -ncol(RNASeq.subset)]))
colnames(RNASeq.subset)[ncol(RNASeq.subset)] <- c("avg")
RNASeq.subset <- merge (RNASeq.subset,Consensus.class,by="row.names")
row.names(RNASeq.subset) <- RNASeq.subset$Row.names
RNASeq.subset$Row.names <- NULL
RNASeq.subset$PatientID <- NULL
RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$Group,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting withing cluster add ,RNASeq.subset$avg 
RNASeq.subset$avg <- NULL
RNASeq.subset$Group <- NULL
Consensus.class<-Consensus.class[rownames(RNASeq.subset),]

# Heatmap
patientcolors <- Consensus.class
levels (patientcolors$Group) <- c(levels (patientcolors$Group),c("#FF0000","#FFA500","#00FF00","#0000FF"))                            #Aply color scheme to patients
patientcolors$Group[patientcolors$Group=="ICR4"] <- "#FF0000"
patientcolors$Group[patientcolors$Group=="ICR3"] <- "#FFA500"
patientcolors$Group[patientcolors$Group=="ICR2"] <- "#00FF00"
patientcolors$Group[patientcolors$Group=="ICR1"] <- "#0000FF"
patientcolors$Group <- droplevels(patientcolors$Group)
patientcolors <- as.character(patientcolors$Group)
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors = c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100))
png("./4 FIGURES/Heatmaps/Heatmap.RNASeq.TCGA.LIHC.ISGS.png",res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(RNASeq.subset),
          main = paste0("Heatmap RNASeq - ",Geneset," sel., K=",K),
          col=my.palette,breaks=my.colors,                                 #set color sheme RED High, GREEN low
          ColSideColors=patientcolors,                                    #set goup colors
          key=TRUE,symm=FALSE,symkey=FALSE,symbreaks=TRUE,             
          scale="row", density.info="none", trace="none",
          labCol=FALSE,cexRow=1.3,cexCol=0.1,margins=c(2,7),
          Colv=FALSE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

