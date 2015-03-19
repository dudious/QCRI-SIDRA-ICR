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
Geneset <- "27G" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 4             # SET K here
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/BRCA.TCGA.MA.k4.",Geneset,".reps1000/BRCA.TCGA.MA.k4.",Geneset,".reps1000.k=",K,".consensusClass.csv"),header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
load (paste0("./2 DATA/SUBSETS/MA.subset.",Geneset,".RData"))
MA.subset <- as.matrix(MA.subset)

#ordeing
MA.subset <- cbind(MA.subset,rowMeans(MA.subset[, -ncol(MA.subset)]))
colnames(MA.subset)[ncol(MA.subset)] <- c("avg")
MA.subset <- merge (MA.subset,Consensus.class,by="row.names")
row.names(MA.subset) <- MA.subset$Row.names
MA.subset$Row.names <- NULL
MA.subset$PatientID <- NULL
MA.subset <- MA.subset[order(factor(MA.subset$Group,levels = c(3,2,1,4)),MA.subset$avg),] 
MA.subset$avg <- NULL
MA.subset$Group <- NULL
Consensus.class<-Consensus.class[rownames(MA.subset),]

# Heatmap
color.map <- function(Consensus.class) { if (Consensus.class=="4") "#FF0000" else "#0000FF" }   #Set color scheme
patientcolors <- unlist(lapply(Consensus.class$Group, color.map))                               #Aply color scheme to patients
my.palette <- colorRampPalette(c("blue", "white", "orange"))(n = 299)
my.colors = c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100))
png(paste0("./4 FIGURES/Heatmaps/Heatmap.MA.R1000.",Geneset,".k=",K,".png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.2(t(MA.subset),
          main = paste0("Heatmap MA - ",Geneset," sel., K=",K),
          col=my.palette,breaks=my.colors,                                 #set color sheme RED High, GREEN low
          ColSideColors=patientcolors,                                    #set goup colors
          key=TRUE,symm=FALSE,symkey=FALSE,symbreaks=TRUE,             
          scale="row", density.info="none", trace="none",
          labCol=FALSE,cexRow=1.3,cexCol=0.1,margins=c(2,7),
          Colv=FALSE)
par(lend = 1)        
legend("topright",legend = c("group 1", "other groups"),
       col = c("red", "blue"),lty= 1,lwd = 5)
dev.off()

