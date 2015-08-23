#################################################################
###
### This Script PLots Heatmaps based on 
### Consensus Clustering clustering of RNASeq Data and mutation data
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
# Dependencies
required.packages <- c("gplots","plyr","beepr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library("plyr")
library("beepr")

## Parameters
Cancerset   = "BRCA.BSF"     # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset     = "DBGS3.FLTR"   # SET GENESET HERE 
matrix.type = "Any"     # Alterantives "Any" , "Missense"
plot.type   = "selected"         # Alterantives "low" , "high" , "373genes"  ,"auto"," selected"


# Load Data
#mutation data
load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))
#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
if (Cancerset %in% c("COAD","READ","UCEC")) {
  #GA data
  Cancerset <- paste0(Cancerset,"-GA")
  Consensus.class.GA <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class.GA <- Consensus.class.GA[,-1]
  colnames (Consensus.class.GA) <- c("Patient_ID","Cluster")
  rownames(Consensus.class.GA) <- Consensus.class.GA[,1]
  Cancerset <- substring(Cancerset,1,4)
  #hiseq data
  Cancerset <- paste0(Cancerset,"-hiseq")
  Consensus.class.hiseq <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class.hiseq <- Consensus.class.hiseq[,-1]
  colnames (Consensus.class.hiseq) <- c("Patient_ID","Cluster")
  rownames(Consensus.class.hiseq) <- Consensus.class.hiseq[,1]
  Cancerset <- substring(Cancerset,1,4)
  #merge GA-hiseq
  Consensus.class <- unique(rbind (Consensus.class.hiseq,Consensus.class.GA))
} else {
  Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class <- Consensus.class[,-1]
  colnames (Consensus.class) <- c("Patient_ID","Cluster")
  rownames(Consensus.class) <- Consensus.class[,1]
}
# select data to plot
if (plot.type == "low"){allmuts.mutatedgenes <- genes.mutations.low}
if (plot.type == "high"){allmuts.mutatedgenes <- genes.mutations.high}
if (plot.type == "auto"){allmuts.mutatedgenes <- genes.mutations.auto}
if (plot.type == "373genes"){allmuts.mutatedgenes <- genes.mutations.373genes}
if (plot.type == "selected"){allmuts.mutatedgenes <- genes.mutations.selected}
allmuts.mutatedgenes[is.na(allmuts.mutatedgenes)] = 0

#merge data
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,Consensus.class,by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,ClinicalData.subset[,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE],by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes$Patient_ID <- NULL

#ordeing rows (patients)
allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                   factor(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"))),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg

#calculate frequency (mean)
allmuts.mutatedgenes.mean <- ddply(allmuts.mutatedgenes[,-ncol(allmuts.mutatedgenes)],.(Cluster),colwise(mean))                             # Calculate frequency (mean)

#generate numeric mutation matrix
allmuts.mutatedgenes$Cluster <- NULL
allmuts.mutatedgenes.mean$Cluster <- NULL
allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq <- NULL
allmuts.mutatedgenes <- as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean <- as.matrix(allmuts.mutatedgenes.mean)
mode(allmuts.mutatedgenes) <- "numeric"
mode(allmuts.mutatedgenes.mean) <- "numeric"

#ordering columns (genes)
allmuts.mutatedgenes.sd <- as.data.frame(apply(allmuts.mutatedgenes.mean,2,sd))                             # Calculate frequency SD 
colnames (allmuts.mutatedgenes.sd) <- c("SD")
allmuts.mutatedgenes.sd <- allmuts.mutatedgenes.sd[order(allmuts.mutatedgenes.sd$SD),,drop = FALSE]

allmuts.mutatedgenes.mean <- as.data.frame(allmuts.mutatedgenes.mean[,rownames(allmuts.mutatedgenes.sd)])   # order mutaion frequency by SD
#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,rownames(allmuts.mutatedgenes.sd)])            # order mutation count by SD freq
allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,order(colnames(allmuts.mutatedgenes))])         # order mutation count alphabeticaly
Consensus.class<-as.data.frame(Consensus.class[rownames(allmuts.mutatedgenes),])                            # sort cluster asignments like mutation matrix

#enforce numeric mutation matrix
allmuts.mutatedgenes = as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean = as.matrix(allmuts.mutatedgenes.mean)
mode(allmuts.mutatedgenes)="numeric"
mode(allmuts.mutatedgenes.mean)="numeric"

#lookup subtype
subtype <- ClinicalData.subset[rownames(allmuts.mutatedgenes),"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#000099","#0066ff","#9933cc","#cc9933"))
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"] <- "#000099"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"] <- "#0066ff"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"] <- "#9933cc"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#cc9933"
subtypecolors <- as.character(subtype$TCGA.PAM50.RMethod.RNASeq)
# Binary Heatmap for selected gene mutations by patient
patientcolors <- Consensus.class
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
#patientcolors$Cluster <- droplevels(patientcolors$cluster)
patientcolors <- as.character(patientcolors$Cluster)

my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".reordered_alphabetic_IMS_A.png"),res=600,height=9,width=25,unit="in")     # set filename
heatmap.2(allmuts.mutatedgenes,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=patientcolors,                        # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2,
          margins=c(10,2),
          labRow=FALSE,
          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
          )
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.5)
dev.off()

my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".reordered_alphabetic_IMS_B.png"),res=600,height=9,width=25,unit="in")     # set filename
heatmap.2(allmuts.mutatedgenes,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=subtypecolors,                        # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2,
          margins=c(10,2),
          labRow=FALSE,
          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched"),
       col = c("#000099","#0066ff","#9933cc","#cc9933"),lty= 1,lwd = 5,cex = 1.5)
dev.off()



