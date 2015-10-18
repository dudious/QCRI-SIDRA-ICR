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
required.packages <- c("gplots","plyr","Hmisc")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library("plyr")
library("Hmisc")

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "chisqr"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict", "chisqr"
IMS.filter     = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "All"         # Alternatives "1vs4" , "All"
gene.filter    = "FALSE"        # Alternatives "TRUE" , "FALSE"

# Load Data
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))
#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]
# load genelist
freq.mut.list <- read.csv ("./3 ANALISYS/Mutations/Frequently mutated cancer genes.csv",header = FALSE)
colnames(freq.mut.list) <- "Gene"
mutsig.list.TCGA <- as.data.frame(read.csv ("./3 ANALISYS/MutSig/brca-tcga-mutsig-genes.csv",header = TRUE))
mutsig.list.GDAC <- as.data.frame(read.csv ("./3 ANALISYS/MutSig/brca-gdac-mutsig-genes.csv",header = TRUE))

#load patient / mutation frequencies
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",IMS.filter,".",Geneset,".Frequencies.RDATA"))

#mutation deciles
Mut.freq.decile <- Mutation.Frequency.Patient[,c("Patient_ID","Freq.NonSilent")]
Mut.freq.decile <- rbind(Mut.freq.decile,c("TCGA-E9-A54Y",0.0))
rownames(Mut.freq.decile) <- Mut.freq.decile$Patient_ID
class(Mut.freq.decile$Freq.NonSilent) <- "numeric"
Mut.freq.decile$decile <- with(Mut.freq.decile,cut(Freq.NonSilent,breaks=quantile(Freq.NonSilent, type=5 , probs=seq(0,1, by=0.1)), 
                                include.lowest=TRUE))
colfunc<-colorRampPalette(rev(c("red","yellow","springgreen","royalblue")))
decile.colors=(colfunc(10))
Mut.freq.decile$label <- Mut.freq.decile$decile
levels(Mut.freq.decile$label) <- decile.colors
Mut.freq.decile$order <- as.numeric(Mut.freq.decile$decile)
#cluster selection
if (cluster.select == "1vs4") {
  Consensus.class <- Consensus.class[Consensus.class$Cluster %in% c("ICR1","ICR4"),]
}

# select data to plot
if (plot.type == "low"){allmuts.mutatedgenes <- genes.mutations.low}
if (plot.type == "high"){allmuts.mutatedgenes <- genes.mutations.high}
if (plot.type == "auto"){allmuts.mutatedgenes <- genes.mutations.auto}
if (plot.type == "db.test"){allmuts.mutatedgenes <- genes.mutations.dbtest}
if (plot.type == "db.test.strict"){allmuts.mutatedgenes <- genes.mutations.dbtest.strict}
if (plot.type == "373genes"){allmuts.mutatedgenes <- genes.mutations.373genes}
if (plot.type == "selected"){allmuts.mutatedgenes <- genes.mutations.selected}
if (plot.type == "chisqr"){allmuts.mutatedgenes <- genes.mutations.chisqr}
allmuts.mutatedgenes[is.na(allmuts.mutatedgenes)] = 0

#filter genes
if (gene.filter == "TRUE") {
  allmuts.mutatedgenes <- allmuts.mutatedgenes[,colnames(allmuts.mutatedgenes) %in% freq.mut.list$Gene]
}

#merge data
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,Consensus.class,by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,ClinicalData.subset[,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE],by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,Mut.freq.decile[,"order",drop=FALSE],by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes$Patient_ID <- NULL

#ordeing rows (patients)
allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                   factor(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                                   allmuts.mutatedgenes$order),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg
#allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
#                                                   factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg

#calculate frequency (mean)
allmuts.mutatedgenes$order <- NULL
if (length(which(is.na(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq)))>0) {
allmuts.mutatedgenes <- allmuts.mutatedgenes[-which(is.na(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq)),]
}
allmuts.mutatedgenes.mean <- ddply(allmuts.mutatedgenes[,-ncol(allmuts.mutatedgenes)],.(Cluster),colwise(mean))                             # Calculate frequency (mean)
allmuts.mutatedgenes.mean.byIMS <- ddply(allmuts.mutatedgenes,.(Cluster,TCGA.PAM50.RMethod.RNASeq),colwise(mean))
Meansorder <- t(allmuts.mutatedgenes.mean[nrow(allmuts.mutatedgenes.mean),-1])-t(allmuts.mutatedgenes.mean[1,-1])
colnames(Meansorder) <- "order.value"
Meansorder[which(Meansorder[,"order.value"]<0),] <- Meansorder[which(Meansorder[,"order.value"]<0),]-1
Meansorder <- as.data.frame(Meansorder[order(-abs(Meansorder)),,drop=FALSE])
Meansorder <- Meansorder[Meansorder$order.value!=0,,drop=FALSE]

#ordering means.byIMS
allmuts.mutatedgenes.mean.byIMS <- allmuts.mutatedgenes.mean.byIMS[order(factor(allmuts.mutatedgenes.mean.byIMS$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                                         factor(allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"))),]
#generate numeric mutation matrix
allmuts.mutatedgenes$Cluster <- NULL
allmuts.mutatedgenes.mean$Cluster <- NULL
Meanby.IMS.colors <- data.frame (Cluster = allmuts.mutatedgenes.mean.byIMS$Cluster,subtype = allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq)
allmuts.mutatedgenes.mean.byIMS$Cluster <- NULL
allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq <- NULL
allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq <- NULL
allmuts.mutatedgenes <- as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean <- as.matrix(allmuts.mutatedgenes.mean)
allmuts.mutatedgenes.mean.byIMS <- as.matrix(allmuts.mutatedgenes.mean.byIMS)
mode(allmuts.mutatedgenes) <- "numeric"
mode(allmuts.mutatedgenes.mean) <- "numeric"
mode(allmuts.mutatedgenes.mean.byIMS) <- "numeric"

#ordering columns (genes)
allmuts.mutatedgenes.sd <- as.data.frame(apply(allmuts.mutatedgenes.mean,2,sd))                             # Calculate frequency SD 
colnames (allmuts.mutatedgenes.sd) <- c("SD")
allmuts.mutatedgenes.sd <- allmuts.mutatedgenes.sd[order(allmuts.mutatedgenes.sd$SD),,drop = FALSE]
allmuts.mutatedgenes.mean <- as.data.frame(allmuts.mutatedgenes.mean[,rownames(allmuts.mutatedgenes.sd)])   # order mutation frequency by SD

#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,rownames(allmuts.mutatedgenes.sd)])            # order mutation count by SD freq
#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,order(colnames(allmuts.mutatedgenes))])         # order mutation count alphabeticaly
allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,rownames(Meansorder)])
allmuts.mutatedgenes.mean.byIMS <- as.data.frame(allmuts.mutatedgenes.mean.byIMS[,rownames(Meansorder)])
Consensus.class<-as.data.frame(Consensus.class[rownames(allmuts.mutatedgenes),])                            # sort cluster asignments like mutation matrix

#enforce numeric mutation matrix
allmuts.mutatedgenes = as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean = as.matrix(allmuts.mutatedgenes.mean)
mode(allmuts.mutatedgenes)="numeric"
mode(allmuts.mutatedgenes.mean)="numeric"

#lookup subtype
subtype <- ClinicalData.subset[rownames(allmuts.mutatedgenes),"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"]     <- "#eaff00"    
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"]     <- "#00c0ff"     
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"]    <- "#da70d6"   
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#daa520"            
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Normal-like"]   <- "#d3d3d3"      
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
#lookup mutation frequency decile
decile <- as.character(Mut.freq.decile[rownames(allmuts.mutatedgenes),"label",drop=TRUE])

#Mutation lists
match.matrix <- cbind(colnames(allmuts.mutatedgenes) %in% freq.mut.list$Gene,
                      colnames(allmuts.mutatedgenes) %in% mutsig.list.TCGA$gene,
                      colnames(allmuts.mutatedgenes) %in% mutsig.list.GDAC$gene)
#colnames(match.matrix) <- c("373g","TCGA","GDAC")
rownames(match.matrix) <- colnames(allmuts.mutatedgenes)
match.matrix[match.matrix[,1]==TRUE,1] <- "#d8d8d8"
match.matrix[match.matrix[,2]==TRUE,2] <- "#b1b1b1"
match.matrix[match.matrix[,3]==TRUE,3] <- "#808080"
match.matrix[match.matrix == FALSE] <- "white"

#Plot size automation
plot.width <- log(ncol(allmuts.mutatedgenes))*4 + 10
plot.height<- 10 

# Heatmap for gene mutation 
color.matrix <- as.matrix(rbind (patientcolors,subtypecolors,decile))
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".ICR4-1mean_IMS_combo_DEC.png"),res=600,height=plot.height,width=plot.width,unit="in")     # set filename
heatmap.3(allmuts.mutatedgenes,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=color.matrix,                         # set goup colors
          ColSideColors=match.matrix,
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=1.8,
          margins=c(10,2),
          labRow=FALSE,
          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topleft",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","",levels(Mut.freq.decile$decile),"","ICR4","ICR3","ICR2","ICR1","","Freq.mut","Mutsig.TCGA","Mutsig.GDAC"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white",levels(Mut.freq.decile$label),"white","red","orange","green","blue","white","#d8d8d8","#b1b1b1","#808080"),lty= 1,lwd = 5,cex = 1)
dev.off()

# Heatmap for gene mutation frequency by cluster/subtype for the same selection of genes

levels(Meanby.IMS.colors$Cluster) <- rev(c("#FF0000","#FFA500","#00FF00","#0000FF"))
levels (Meanby.IMS.colors$subtype) <- c(levels (Meanby.IMS.colors$subtype),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
Meanby.IMS.colors$subtype[Meanby.IMS.colors$subtype=="Luminal A"]     <- "#eaff00"    
Meanby.IMS.colors$subtype[Meanby.IMS.colors$subtype=="Luminal B"]     <- "#00c0ff"     
Meanby.IMS.colors$subtype[Meanby.IMS.colors$subtype=="Basal-like"]    <- "#da70d6"   
Meanby.IMS.colors$subtype[Meanby.IMS.colors$subtype=="HER2-enriched"] <- "#daa520"            
Meanby.IMS.colors$subtype[Meanby.IMS.colors$subtype=="Normal-like"]   <- "#d3d3d3"
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors <- c(seq(0,0.01,length=100),seq(0.01,0.05,length=100),seq(0.05,1,length=100))
Meanby.IMS.colors <- as.matrix(t(Meanby.IMS.colors))
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".ICR4-1mean_IMS_MEAN.png"),res=600,height=plot.height,width=plot.width+8,unit="in")     # set filename
par(mar=c(1,1,1,1))
heatmap.3(allmuts.mutatedgenes.mean.byIMS,
          main = "HeatMap-MutatedGenes frequency",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          #breaks=my.colors,                                   # set manual color gradient of color scheme
          RowSideColors=Meanby.IMS.colors,                        # set goup colors
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="col", 
          #density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2.6,
          margins=c(15,15),
          labRow=FALSE,
          Colv=FALSE,Rowv=FALSE                               # reorder row/columns by dendogram
)
par(lend = 1)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue"),lty= 1,lwd = 5,cex = 1)
dev.off()

if (plot.type == "selected"){
# Heatmap showing mutual exlusion
allmuts.exclusion.matrix <- cbind (t(color.matrix),allmuts.mutatedgenes)
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[,-3]
allmuts.exclusion.matrix.reordered <- allmuts.exclusion.matrix[allmuts.exclusion.matrix[,"MAP3K1"]==1,]
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[-which(rownames(allmuts.exclusion.matrix) %in% rownames(allmuts.exclusion.matrix.reordered)),]
allmuts.exclusion.matrix.reordered <- rbind(allmuts.exclusion.matrix.reordered,allmuts.exclusion.matrix[allmuts.exclusion.matrix[,"MAP2K4"]==1,])
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[-which(rownames(allmuts.exclusion.matrix) %in% rownames(allmuts.exclusion.matrix.reordered)),]
allmuts.exclusion.matrix.reordered <- rbind(allmuts.exclusion.matrix.reordered,allmuts.exclusion.matrix[allmuts.exclusion.matrix[,"CTCF"]==1,])
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[-which(rownames(allmuts.exclusion.matrix) %in% rownames(allmuts.exclusion.matrix.reordered)),]
allmuts.exclusion.matrix.reordered <- rbind(allmuts.exclusion.matrix.reordered,allmuts.exclusion.matrix[allmuts.exclusion.matrix[,"FCGBP"]==1,])
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[-which(rownames(allmuts.exclusion.matrix) %in% rownames(allmuts.exclusion.matrix.reordered)),]
allmuts.exclusion.matrix.reordered <- rbind(allmuts.exclusion.matrix.reordered,allmuts.exclusion.matrix[allmuts.exclusion.matrix[,"TP53"]==1,])
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[-which(rownames(allmuts.exclusion.matrix) %in% rownames(allmuts.exclusion.matrix.reordered)),]
allmuts.exclusion.matrix.reordered <- rbind(allmuts.exclusion.matrix.reordered,allmuts.exclusion.matrix)
allmuts.exclusion.matrix <- allmuts.exclusion.matrix[-which(rownames(allmuts.exclusion.matrix) %in% rownames(allmuts.exclusion.matrix.reordered)),]

color.matrix.exlusion <- as.matrix(t(allmuts.exclusion.matrix.reordered[,c("patientcolors","subtypecolors")]))
allmuts.exclusion.matrix.reordered <- as.matrix(allmuts.exclusion.matrix.reordered[,-c(1,2)])
mode(allmuts.exclusion.matrix.reordered) <- "numeric"
colnames(color.matrix.exlusion) <- NULL

my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".mutual_exclusion.png"),res=600,height=plot.height,width=plot.width,unit="in")     # set filename
heatmap.3(allmuts.exclusion.matrix.reordered,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=color.matrix.exlusion,                         # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes),
          cexRow=1,cexCol=2.2,
          margins=c(10,2),
          labRow=FALSE,
          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
#legend("topleft",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","",levels(Mut.freq.decile$decile),"","ICR4","ICR3","ICR2","ICR1"),
#       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white",levels(Mut.freq.decile$label),"white","red","orange","green","blue"),lty= 1,lwd = 5,cex = 1)
dev.off()
}