#################################################################
###
### This Script generates a heatmap from the 
### cancer hallmark patways gene mutations
###
#################################################################

# Setup environment
rm(list=ls())
#setwd("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")

#Dependencies
required.packages <- c("gplots","plyr","fmsb")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library("plyr")
library("fmsb")

# Set Parameters
Cancerset <- "OV"           # do not use -GA or -hiseq (data is merged)
BRCA.Filter <- "PCF"         # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
IMS.filter = "All"            # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
matrix.type    = "NonSilent"         # Alterantives "Any" , "Missense", "NonSilent"
statfilter = "TRUE"
Binary = "TRUE"

if (Cancerset =="BRCA") {
  Cancerset <- paste0(Cancerset,".",BRCA.Filter)
}

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
#load mutation data
load (file =  paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Hallmark.cancer.pathway.Mutation.Matrix.",matrix.type,".Rdata"))

#stats Mann-Whitney-Wilcoxon Test
test.matrix <- rbind (Consensus.class[Consensus.class$Cluster=="ICR1" & Consensus.class$Patient_ID %in% rownames(hallmark.pathways.Zscores),"Cluster",drop=FALSE],
                      Consensus.class[Consensus.class$Cluster=="ICR4" & Consensus.class$Patient_ID %in% rownames(hallmark.pathways.Zscores),"Cluster",drop=FALSE])
test.matrix <- merge(test.matrix,hallmark.pathways.Zscores,by = "row.names",all.x = TRUE,all.y = FALSE)
rownames(test.matrix) <- test.matrix$Row.names
test.matrix$Row.names <- NULL
test.matrix <- test.matrix[order(test.matrix$Cluster),]
stats <- apply(test.matrix[-1],2,function(x) wilcox.test(x ~ Cluster, data=test.matrix))
stats.df <- as.data.frame(t(as.data.frame(lapply(stats, function(x) {x$p.value}))))
colnames(stats.df) <- "p.value"
stats.df.sign <- stats.df[stats.df$p.value < 0.05,,drop=FALSE]
stats.df.sign <- stats.df.sign[order(stats.df.sign$p.value),,drop=FALSE]

#Pathway oddsratio
ICR1.count <- apply(test.matrix[test.matrix$Cluster=="ICR1",c(2:ncol(test.matrix))],2,function (x) sum(x > 0))
ICR4.count <- apply(test.matrix[test.matrix$Cluster=="ICR4",c(2:ncol(test.matrix))],2,function (x) sum(x > 0))
count.table <- cbind (ICR1.count,ICR4.count,nrow(test.matrix[test.matrix$Cluster=="ICR1",]),nrow(test.matrix[test.matrix$Cluster=="ICR4",]))
colnames(count.table)<- c("ICR1_MUT","ICR4_MUT","ICR1_N","ICR4_N")
count.table<- as.data.frame(count.table)
count.table$ICR1_WT <- count.table$ICR1_N-count.table$ICR1_MUT
count.table$ICR4_WT <- count.table$ICR4_N-count.table$ICR4_MUT
count.table$pathways <- rownames(count.table)
Pathway.oddsratios <- ddply(count.table, .(pathways), .fun=function(x) {v1=oddsratio(x$ICR1_MUT,x$ICR1_WT,x$ICR4_MUT,x$ICR4_WT);v2=c(v1$estimate, v1$p.value, v1$conf.int); names(v2) = c("ratio", "p_value","Lower", "Upper"); return(v2)})
Pathway.oddsratios <- Pathway.oddsratios[order(Pathway.oddsratios$p_value),]
Pathway.oddsratios$ratio_inv <- 1/Pathway.oddsratios$ratio

save(count.table,Pathway.oddsratios,file=paste0("./3 ANALISYS/Pathway ORs/Pathway.mutation.ORs.",Cancerset,".",Geneset,".Rdata"))
#merge data
allmuts.mutatedgenes <- merge (hallmark.pathways.Zscores,Consensus.class,by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes$Patient_ID <- NULL

#ordeing rows (patients)
allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg

#filter by stats
if (statfilter == "TRUE"){
  allmuts.mutatedgenes <- allmuts.mutatedgenes[allmuts.mutatedgenes$Cluster=="ICR4"|allmuts.mutatedgenes$Cluster=="ICR1" ,c("Cluster",rownames(stats.df.sign))]
}

#sidecolors
patientcolors <- allmuts.mutatedgenes[,"Cluster",drop=FALSE]
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

#Make Binary
if (Binary == "TRUE"){
  allmuts.mutatedgenes[allmuts.mutatedgenes > 1] <- 1
}



my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 20)
#dev.new()
png(paste0("./4 FIGURES/Heatmaps/selected_pathways/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.StatFilter_",statfilter,".Binary_",Binary,".png"),res=600,height=25,width=25,unit="in")     # set filename
heatmap.2(allmuts.mutatedgenes,
          main = paste0("Cancer Selected Pathways - ",Cancerset),
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=patientcolors,                        # set goup colors
          key=FALSE,
          #symm=FALSE,
          #symkey=FALSE,
          #symbreaks=TRUE,             
          #scale="col", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2,
          margins=c(46,2),
          labRow=FALSE,
          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
legend("bottomleft",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.5)
dev.off()
