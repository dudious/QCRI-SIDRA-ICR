#################################################################
###
### This Script generates a heatmap from the 
### cancer hallmark patways gene mutations
###
#################################################################

# Setup environment
rm(list=ls())

#Dependencies
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")

# Set Parameters
Cancerset <- "OV"           # do not use -GA or -hiseq (data is merged)
BRCA.Filter <- "BSF2"         # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
IMS.filter = "All"            # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
matrix.type    = "NonSilent"         # Alterantives "Any" , "Missense", "NonSilent"

#load mutation data
load (file =  paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Hallmark.cancer.pathway.Mutation.Matrix.",matrix.type,".Rdata"))

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

#merge data
allmuts.mutatedgenes <- merge (hallmark.pathways.Zscores,Consensus.class,by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes$Patient_ID <- NULL

#ordeing rows (patients)
allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg

#sidecolors
patientcolors <- allmuts.mutatedgenes[,ncol(allmuts.mutatedgenes),drop=FALSE]
allmuts.mutatedgenes$Cluster <- NULL
levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
#patientcolors$Cluster <- droplevels(patientcolors$cluster)
patientcolors <- as.character(patientcolors$Cluster)

#enforce numeric mutation matrix
allmuts.mutatedgenes = as.matrix(allmuts.mutatedgenes)
mode(allmuts.mutatedgenes)="numeric"


my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
dev.new()
#png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".reordered_alphabetic.png"),res=600,height=9,width=25,unit="in")     # set filename
heatmap.2(allmuts.mutatedgenes,
          main = "Cancer Hallmark Pathways",
          #col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=patientcolors,                        # set goup colors
          #key=FALSE,
          #symm=FALSE,
          #symkey=FALSE,
          #symbreaks=TRUE,             
          #scale="col", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2,
          margins=c(45,2),
          labRow=FALSE,
          Colv=TRUE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
legend("bottomleft",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.5)
dev.off()
