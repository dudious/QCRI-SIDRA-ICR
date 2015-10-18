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

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "db.test"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict"
IMS.filter     = "Basal"     # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "1vs4"        # Alternatives "1vs4" , "All"
method         = "GO"

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
allmuts.mutatedgenes[is.na(allmuts.mutatedgenes)] = 0

#select all data
#allmuts.mutatedgenes <- genes.mutations

#merge data
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,Consensus.class,by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes <- merge (allmuts.mutatedgenes,ClinicalData.subset[,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE],by="row.names")
row.names(allmuts.mutatedgenes) <- allmuts.mutatedgenes$Row.names
allmuts.mutatedgenes$Row.names <- NULL
allmuts.mutatedgenes$Patient_ID <- NULL

#ordeing rows (patients)
allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                                   factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]    # if sorting within cluster add ,allmuts.mutatedgenes$avg
if (IMS.filter == "All") {
  allmuts.mutatedgenes <- allmuts.mutatedgenes[order(factor(allmuts.mutatedgenes$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                   factor(allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"))),] 
}
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
levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#000099","#0066ff","#9933cc","#cc9933","#000000"))
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"] <- "#000099"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"] <- "#0066ff"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"] <- "#9933cc"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#cc9933"
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Normal-like"] <- "#000000"
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

#my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
#png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".reordered_alphabetic_IMS_A.png"),res=600,height=9,width=25,unit="in")     # set filename
#heatmap.2(allmuts.mutatedgenes,
#          main = "HeatMap-MutatedGenes",
#         col=my.palette,                                     # set color scheme RED High, GREEN low
#          RowSideColors=patientcolors,                        # set goup colors
#          key=FALSE,
#          symm=FALSE,
#          symkey=FALSE,
#          symbreaks=TRUE,             
#          #scale="row", 
#          density.info="none",
#          trace="none",
#          labCol=colnames(allmuts.mutatedgenes),
#          cexRow=1,cexCol=2,
#          margins=c(10,2),
#         labRow=FALSE,
#          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
#          )
#par(lend = 1)
#legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.5)
#dev.off()
#
#my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 3)
#png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".reordered_alphabetic_IMS_B.png"),res=600,height=9,width=25,unit="in")     # set filename
#heatmap.2(allmuts.mutatedgenes,
#          main = "HeatMap-MutatedGenes",
#          col=my.palette,                                     # set color scheme RED High, GREEN low
#          RowSideColors=subtypecolors,                        # set goup colors
#          key=FALSE,
#          symm=FALSE,
#          symkey=FALSE,
#          symbreaks=TRUE,             
#          #scale="row", 
#          density.info="none",
#          trace="none",
#          labCol=colnames(allmuts.mutatedgenes),
#          cexRow=1,cexCol=2,
#          margins=c(10,2),
#          labRow=FALSE,
#          Colv=FALSE, Rowv=FALSE                              # reorder row/columns by dendogram
# )
#par(lend = 1)
#legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"),
#       col = c("#000099","#0066ff","#9933cc","#cc9933","#000000"),lty= 1,lwd = 5,cex = 1.3)
#dev.off()

#convertion of gene to GO
save (allmuts.mutatedgenes,file="./3 ANALISYS/Mutations/BRCA.BSF2/Binary_mutation_matrix.All.RDAta")
load ("./3 ANALISYS/Mutations/BRCA.BSF2/Binary_mutation_matrix.All.RDAta")

pathways <- read.csv("./3 ANALISYS/Mutations/BRCA.BSF2/Enriched.pathways.csv")
splitted.genes <- strsplit(as.character(pathways$genes), "/")
pathways.melted <- data.frame(ID = rep.int(pathways$ID, sapply(splitted.genes, length)), Gene = unlist(splitted.genes))
PW.matrix <- as.data.frame(t(allmuts.mutatedgenes))
PW.matrix.2 <- merge(pathways.melted,PW.matrix,by.x="Gene",by.y="row.names",all.x=TRUE,All.y=TRUE)
PW.matrix.2[is.na(PW.matrix.2)] <- 0
PW.matrix.2$Gene <-NULL
PW.matrix.summed <- aggregate(PW.matrix.2,by=list(PW.matrix.2$ID), FUN = sum)
PW.matrix.summed$ID <-NULL
PW.genecounts <- as.data.frame(table (pathways.melted$ID))
PW.matrix.summed$PW.norm <- PW.genecounts$Freq[match(PW.matrix.summed$Group.1,PW.genecounts$Var1)]
PW.matrix.summed.norm <- (PW.matrix.summed / PW.matrix.summed$PW.norm)*100
PW.matrix.summed.norm$Group.1 <- (PW.matrix.summed.norm$Group.1 /100)*PW.matrix.summed$PW.norm
PW.matrix.summed$PW.norm <- NULL
PW.matrix.summed.norm$PW.norm <- NULL
PW.matrix.summed$PW.name <- pathways$class[match(PW.matrix.summed$Group.1,pathways$ID)]
PW.matrix.summed.norm$PW.name <- pathways$class[match(PW.matrix.summed.norm$Group.1,pathways$ID)]
rownames(PW.matrix.summed) <- PW.matrix.summed$PW.name
rownames(PW.matrix.summed.norm) <- PW.matrix.summed$PW.name
PW.matrix.summed$Group.1 <-NULL
PW.matrix.summed$PW.name <-NULL
allmuts.mutated.PW <- t(PW.matrix.summed)
PW.matrix.summed.norm$Group.1 <-NULL
PW.matrix.summed.norm$PW.name <-NULL
allmuts.mutated.PW.norm <- t(PW.matrix.summed.norm)


load ("./3 ANALISYS/Mutations/BRCA.BSF2/diferentially.mutated.dbtest.loose.GO.Rdata")
GO.matrix <- as.data.frame(t(allmuts.mutatedgenes))
GO.names <- unique(GO.data.filtered[,c(2,3)])
#GO.matrix$GO <- GO.data.filtered$name_1006[match(rownames(GO.matrix),GO.data.filtered$hgnc_symbol)]
GO.matrix.2 <- merge(GO.data.filtered,GO.matrix,by.x="hgnc_symbol",by.y="row.names",all.x=TRUE,All.y=TRUE)
GO.matrix.2[is.na(GO.matrix.2)] <- 0
GO.matrix.2$hgnc_symbol <-NULL
GO.matrix.2$name_1006 <-NULL
GO.matrix.summed <- aggregate(data=GO.matrix.2,.~go_id, FUN = sum)
GO.genecounts <- as.data.frame(table (GO.matrix.2$go_id))
GO.matrix.summed$GO.norm <- GO.genecounts$Freq[match(GO.matrix.summed$go_id,GO.genecounts$Var1)]
GO.matrix.summed.norm <- (GO.matrix.summed[,-1] / GO.matrix.summed$GO.norm)*100
GO.matrix.summed$GO.name <- GO.names$name_1006[match(GO.matrix.summed$go_id,GO.names$go_id)]
GO.matrix.summed.norm$GO.norm <- NULL
GO.matrix.summed$GO.norm <- NULL
GO.matrix.summed$go_id <- NULL
rownames(GO.matrix.summed) <- GO.matrix.summed$GO.name
rownames(GO.matrix.summed.norm) <- GO.matrix.summed$GO.name
GO.matrix.summed$GO.name <- NULL

allmuts.mutated.GO <- t(GO.matrix.summed)
allmuts.mutated.GO.norm <- t(GO.matrix.summed.norm)

allmuts.mutated.plot <- allmuts.mutated.PW
if (method == "GO"){ allmuts.mutated.plot <- allmuts.mutated.GO }

#pathway stats & filter
allmuts.mutated.stats <- as.data.frame(allmuts.mutated.plot)
allmuts.mutated.stats$cluster <- Consensus.class$Cluster[match(rownames(allmuts.mutated.stats),Consensus.class$Patient_ID)]
Stats <- as.data.frame(t(sapply(allmuts.mutated.stats[-ncol(allmuts.mutated.stats)], function(x) 
  unlist(t.test(x~allmuts.mutated.stats$cluster)[c("estimate","p.value","statistic","conf.int")]))))
nrow(Stats)
significant.subset <- rownames(Stats[Stats$p.value<0.01,])
length(significant.subset)
allmuts.mutated.plot <- allmuts.mutated.plot[,colnames(allmuts.mutated.plot)%in%significant.subset]


color.matrix <- as.matrix(rbind (patientcolors,subtypecolors))
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(0,2,length=100),seq(2,15,length=100),seq(15,100,length=100)))
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".IMS_clustered_combo_",method,".ttested.png"),
    res=300,height=25,width=25,unit="in")     # set filename
heatmap.3(allmuts.mutated.plot,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          breaks=my.colors,
          RowSideColors=color.matrix,                         # set goup colors
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="col", 
          density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutated.plot),
          cexRow=1,cexCol=2,
          margins=c(65,2),
          labRow=FALSE,
          Colv=TRUE, Rowv=FALSE                              # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topleft",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1"),
       col = c("#000099","#0066ff","#9933cc","#cc9933","#000000","white","red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
dev.off()

