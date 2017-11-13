#################################################################
###
### This draws immunescore heatmaps by cancer 
###
### 
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")
#setwd("/mnt3/wouter/BREAST-QATAR/")
#Dependencies
source ("~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

# Set Parameters
Cancersets        = "BRCA"
Geneset           = "DBGS3.FLTR"   
DL.Method         = "BIOLINKS"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Selected"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer
Genedatabase      = "Gene_selection_v2.6.txt"

#Load data
dir.create("./3 ANALISYS/IMMUNOSCORE/xCell scoring/",showWarnings=FALSE)
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
gene.list <- read.csv (paste0("./2 DATA/SUBSETS/",Genedatabase))                                 # Select subset here !!!!! and change filename below !!!!
gene.list.ALL <- as.character(gene.list[which(gene.list[,"DBGS3"]==1),1])
gene.list.INH <- as.character(gene.list[which(gene.list[,"ImSuGS"]==1),1])
a<-which(gene.list[,"DBGS3"]==1)
b<-which(gene.list[,"ImSuGS"]==1)
gene.list.ACT <- as.character(gene.list[a[-which(a%in%b)],1])

#test
Cancerset = Cancersets

# DO ALL
#for (i in 1:N.sets) {
#  Cancerset = Cancersets[i]
#  if (Cancerset %in% c("LAML","FPPP")) {next}

## Load Data
xCell.data <- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/xCell scoring/xCell_FINAL_",Cancerset,".RNASeq.TCGA.",DL.Method,".Selected.NORMALIZED.TP_FILTERED_LOG2_xCell_0556040617.txt"),sep = "\t")
rownames(xCell.data) <- xCell.data$X
xCell.data$X <- NULL
colnames(xCell.data) <- gsub("\\.","-",colnames(xCell.data))
xCell.matrix <- as.matrix (xCell.data)


Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class$X <- NULL
rownames(Consensus.class) <- Consensus.class$PatientID
#filter ICR1 and 4
Consensus.class<-Consensus.class[Consensus.class$Group %in% c("ICR1","ICR4"),]
Consensus.class<-Consensus.class[order(Consensus.class$Group),]
xCell.matrix <- xCell.matrix[,row.names(Consensus.class)]
Consensus.class <- Consensus.class[colnames(xCell.matrix),]


patientcolors <- Consensus.class
patientcolors$PatientID <- NULL
levels (patientcolors$Group) <- c(levels (patientcolors$Group),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to Group
patientcolors$Group[patientcolors$Group=="ICR4"] <- "#FF0000"
patientcolors$Group[patientcolors$Group=="ICR3"] <- "#FFA500"
patientcolors$Group[patientcolors$Group=="ICR2"] <- "#00FF00"
patientcolors$Group[patientcolors$Group=="ICR1"] <- "#0000FF"

patientcolors <- as.matrix(patientcolors)

my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
my.colors = unique(c(seq(-5,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,5,length=100)))

#scale
xCell.matrix <- t(scale(t(xCell.matrix)))

dev.new()
heatmap.3(xCell.matrix,
          main = paste0("Heatmap immune enrichment in",Cancerset),
          col=my.palette,
          breaks = my.colors,
          #scale = "row",
          #Colv=FALSE,
          ColSideColors=patientcolors,
          labCol=FALSE,
          cexRow=1.3,cexCol=0.1,margins=c(2,7)
)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR1"),
       col = c("red","blue"),lty= 1,lwd = 5,cex = 0.7)

