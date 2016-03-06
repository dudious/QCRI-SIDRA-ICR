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
plot.type      = "selected"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict", "chisqr"
IMS.filter     = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "All"         # Alternatives "1vs4" , "All"
gene.filter    = "FALSE"        # Alternatives "TRUE" , "FALSE"
selected.genes = c("TP53","MAP3K1","MAP2K4","CTCF","FCGBP","CXCL9")

# Load Data
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
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
Mut.freq.decile$HiLo <- Mut.freq.decile$order
Mut.freq.decile$HiLo[Mut.freq.decile$HiLo>7] <- "High"
Mut.freq.decile$HiLo[Mut.freq.decile$HiLo<8 & Mut.freq.decile$HiLo>3] <- "Medium"
Mut.freq.decile$HiLo[Mut.freq.decile$HiLo<4] <- "Low"

#cluster selection
if (cluster.select == "1vs4") {
  Consensus.class <- Consensus.class[Consensus.class$Cluster %in% c("ICR1","ICR4"),]
}

## Read the gene list (373 genes)
genes.list = read.table("./2 DATA/Frequently mutated cancer genes.csv")

## Load the mutation variation data
load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".",matrix.type,".VariationTables.RData"))

# select data to plot
## Pick the genes from the provided (373 genes)list or significant variation File
genes.mutations.373genes = genes.mutations[,colnames(genes.mutations) %in% as.character(genes.list$V1)]
genes.mutations.low = genes.mutations[,colnames(genes.mutations) %in% as.character(low.significant.variation.table$Gene)]
genes.mutations.high = genes.mutations[,colnames(genes.mutations) %in% as.character(high.significant.variation.table$Gene)]
genes.mutations.auto = genes.mutations[,colnames(genes.mutations) %in% as.character(auto.significant.variation.table$Gene)]
genes.mutations.selected = genes.mutations[,colnames(genes.mutations) %in% selected.genes]
genes.mutations.dbtest = genes.mutations[,colnames(genes.mutations) %in% db.test.significant.variation.table$Gene]
genes.mutations.dbtest.strict = genes.mutations[,colnames(genes.mutations) %in% db.test.strict.significant.variation.table$Gene]
genes.mutations.chisqr = genes.mutations[,colnames(genes.mutations) %in% chisq.significant.variation.table$Gene]

if (plot.type == "low"){allmuts.mutatedgenes <- genes.mutations.low}
if (plot.type == "high"){allmuts.mutatedgenes <- genes.mutations.high}
if (plot.type == "auto"){allmuts.mutatedgenes <- genes.mutations.auto}
if (plot.type == "db.test"){allmuts.mutatedgenes <- genes.mutations.dbtest}
if (plot.type == "db.test.strict"){allmuts.mutatedgenes <- genes.mutations.dbtest.strict}
if (plot.type == "373genes"){allmuts.mutatedgenes <- genes.mutations.373genes}
if (plot.type == "selected"){allmuts.mutatedgenes <- genes.mutations.selected}
if (plot.type == "chisqr"){allmuts.mutatedgenes <- genes.mutations.chisqr}
allmuts.mutatedgenes[allmuts.mutatedgenes=="MUT"] = 1
allmuts.mutatedgenes[is.na(allmuts.mutatedgenes)] = 0
allmuts.mutatedgenes <- as.matrix(allmuts.mutatedgenes)
mode(allmuts.mutatedgenes) <- "numeric"

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
allmuts.mutatedgenes$HiLo <- Mut.freq.decile$HiLo[match(rownames(allmuts.mutatedgenes),rownames(Mut.freq.decile))]
allmuts.mutatedgenes.mean.byIMS.byHilo <- ddply(allmuts.mutatedgenes,.(Cluster,TCGA.PAM50.RMethod.RNASeq,HiLo),colwise(mean))
allmuts.mutatedgenes.meta <- allmuts.mutatedgenes[,(ncol(allmuts.mutatedgenes)-2):ncol(allmuts.mutatedgenes)]
allmuts.mutatedgenes$HiLo <- NULL
Meansorder <- t(allmuts.mutatedgenes.mean[nrow(allmuts.mutatedgenes.mean),-1])-t(allmuts.mutatedgenes.mean[1,-1])
colnames(Meansorder) <- "order.value"
Meansorder[which(Meansorder[,"order.value"]<0),] <- Meansorder[which(Meansorder[,"order.value"]<0),]-1
Meansorder <- as.data.frame(Meansorder[order(-abs(Meansorder)),,drop=FALSE])
#Meansorder <- Meansorder[Meansorder$order.value!=0,,drop=FALSE]

#ordering means.byIMS
allmuts.mutatedgenes.mean.byIMS <- allmuts.mutatedgenes.mean.byIMS[order(factor(allmuts.mutatedgenes.mean.byIMS$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                                         factor(allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"))),]
allmuts.mutatedgenes.mean.byIMS.byHilo <- allmuts.mutatedgenes.mean.byIMS.byHilo[order(factor(allmuts.mutatedgenes.mean.byIMS.byHilo$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                                                       factor(allmuts.mutatedgenes.mean.byIMS.byHilo$TCGA.PAM50.RMethod.RNASeq,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                                                                       factor(allmuts.mutatedgenes.mean.byIMS.byHilo$HiLo,levels = c("Low","Medium","High"))),]

#generate numeric mutation matrix
allmuts.mutatedgenes$Cluster <- NULL
allmuts.mutatedgenes.mean$Cluster <- NULL
allmuts.mutatedgenes.mean.byIMS <- allmuts.mutatedgenes.mean.byIMS[allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq != "Normal-like",] #remove normal like for meansplot
allmuts.mutatedgenes.mean.byIMS.byHilo <- allmuts.mutatedgenes.mean.byIMS.byHilo[allmuts.mutatedgenes.mean.byIMS.byHilo$TCGA.PAM50.RMethod.RNASeq != "Normal-like",] #remove normal like for meansplot
Meanby.IMS.colors <- data.frame (Cluster = allmuts.mutatedgenes.mean.byIMS$Cluster,subtype = allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq)
Meanby.IMS.Hilo.colors <- data.frame (Cluster = allmuts.mutatedgenes.mean.byIMS.byHilo$Cluster,
                                      subtype = allmuts.mutatedgenes.mean.byIMS.byHilo$TCGA.PAM50.RMethod.RNASeq,
                                      Hilo = allmuts.mutatedgenes.mean.byIMS.byHilo$HiLo)
allmuts.mutatedgenes.mean.byIMS$Cluster <- NULL
allmuts.mutatedgenes.mean.byIMS$TCGA.PAM50.RMethod.RNASeq <- NULL
allmuts.mutatedgenes.mean.byIMS.byHilo$Cluster <- NULL
allmuts.mutatedgenes.mean.byIMS.byHilo$TCGA.PAM50.RMethod.RNASeq <- NULL
allmuts.mutatedgenes.mean.byIMS.byHilo$HiLo <- NULL
allmuts.mutatedgenes$TCGA.PAM50.RMethod.RNASeq <- NULL
allmuts.mutatedgenes <- as.matrix(allmuts.mutatedgenes)
allmuts.mutatedgenes.mean <- as.matrix(allmuts.mutatedgenes.mean)
allmuts.mutatedgenes.mean.byIMS <- as.matrix(allmuts.mutatedgenes.mean.byIMS)
allmuts.mutatedgenes.mean.byIMS.byHilo <- as.matrix(allmuts.mutatedgenes.mean.byIMS.byHilo)
mode(allmuts.mutatedgenes) <- "numeric"
mode(allmuts.mutatedgenes.mean) <- "numeric"
mode(allmuts.mutatedgenes.mean.byIMS) <- "numeric"
mode(allmuts.mutatedgenes.mean.byIMS.byHilo) <- "numeric"

#ordering columns (genes)
allmuts.mutatedgenes.sd <- as.data.frame(apply(allmuts.mutatedgenes.mean,2,sd))                             # Calculate frequency SD 
colnames (allmuts.mutatedgenes.sd) <- c("SD")
allmuts.mutatedgenes.sd <- allmuts.mutatedgenes.sd[order(allmuts.mutatedgenes.sd$SD),,drop = FALSE]
allmuts.mutatedgenes.mean <- as.data.frame(allmuts.mutatedgenes.mean[,rownames(allmuts.mutatedgenes.sd)])   # order mutation frequency by SD

#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,rownames(allmuts.mutatedgenes.sd)])            # order mutation count by SD freq
#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,order(colnames(allmuts.mutatedgenes))])         # order mutation count alphabeticaly
allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,rownames(Meansorder)])
#allmuts.mutatedgenes <- as.data.frame(allmuts.mutatedgenes[,c("TP53","MAP3K1","CXCL9")])                    # manual overwrite
#order depending tables
allmuts.mutatedgenes.mean.byIMS <- as.data.frame(allmuts.mutatedgenes.mean.byIMS[,rownames(Meansorder)])
allmuts.mutatedgenes.mean.byIMS.byHilo <- as.data.frame(allmuts.mutatedgenes.mean.byIMS.byHilo[,rownames(Meansorder)])
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
my.palette <- colorRampPalette(c("blue", "#c4f5ff", "#8b0000"))(n = 3)
#dev.new()
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".ICR4-1mean_IMS_combo_DEC.v2.png"),res=600,height=plot.height,width=(plot.width/2)+16,unit="in")     # set filename
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
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".ICR4-1mean_IMS_MEAN.png"),res=600,height=plot.height,width=(plot.width/2)+30,unit="in")     # set filename
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

# Heatmap for gene mutation frequency by cluster/subtype for the same selection of genes HiLo edition

levels(Meanby.IMS.Hilo.colors$Cluster) <- rev(c("#FF0000","#FFA500","#00FF00","#0000FF"))
levels (Meanby.IMS.Hilo.colors$subtype) <- c(levels (Meanby.IMS.Hilo.colors$subtype),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
levels (Meanby.IMS.Hilo.colors$Hilo) <- c(levels (Meanby.IMS.Hilo.colors$Hilo),c("green","yellow","red"))
Meanby.IMS.Hilo.colors$Hilo[Meanby.IMS.Hilo.colors$Hilo=="Low"]     <- "green"    
Meanby.IMS.Hilo.colors$Hilo[Meanby.IMS.Hilo.colors$Hilo=="Medium"]  <- "yellow"     
Meanby.IMS.Hilo.colors$Hilo[Meanby.IMS.Hilo.colors$Hilo=="High"]    <- "red"
Meanby.IMS.Hilo.colors$subtype[Meanby.IMS.Hilo.colors$subtype=="Luminal A"]     <- "#eaff00"    
Meanby.IMS.Hilo.colors$subtype[Meanby.IMS.Hilo.colors$subtype=="Luminal B"]     <- "#00c0ff"     
Meanby.IMS.Hilo.colors$subtype[Meanby.IMS.Hilo.colors$subtype=="Basal-like"]    <- "#da70d6"   
Meanby.IMS.Hilo.colors$subtype[Meanby.IMS.Hilo.colors$subtype=="HER2-enriched"] <- "#daa520"            
Meanby.IMS.Hilo.colors$subtype[Meanby.IMS.Hilo.colors$subtype=="Normal-like"]   <- "#d3d3d3"
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
my.colors <- c(seq(0,0.01,length=100),seq(0.01,0.05,length=100),seq(0.05,1,length=100))
Meanby.IMS.Hilo.colors <- as.matrix(t(Meanby.IMS.Hilo.colors))
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".ICR4-1mean_IMS_MEAN_HiLo.png"),res=600,height=plot.height,width=(plot.width/2)+30,unit="in")     # set filename
par(mar=c(1,1,1,1))
heatmap.3(allmuts.mutatedgenes.mean.byIMS.byHilo,
          main = "HeatMap-MutatedGenes frequency",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          #breaks=my.colors,                                   # set manual color gradient of color scheme
          RowSideColors=Meanby.IMS.Hilo.colors,                        # set goup colors
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="col", 
          #density.info="none",
          trace="none",
          labCol=colnames(allmuts.mutatedgenes.mean.byIMS.byHilo),
          cexRow=1,cexCol=2.6,
          margins=c(15,15),
          labRow=FALSE,
          Colv=FALSE,Rowv=FALSE                               # reorder row/columns by dendogram
        )
par(lend = 1)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1","","Low (1-3)","Medium (4-7)","High (8-10)"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue","white","green","yellow","red"),lty= 1,lwd = 5,cex = 0.9)
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
  png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".mutual_exclusion.png"),res=600,height=plot.height,width=plot.width/2,unit="in")     # set filename
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

# Heatmap for matching GISTC data
gistic.genes <- c(colnames(allmuts.exclusion.matrix.reordered))
selected.data <- merged.matrix[rownames(allmuts.mutatedgenes),gistic.genes]
selected.data[selected.data=="MUT"] <- NA
selected.data <- gsub("MUT;","",as.matrix(selected.data))
selected.data[selected.data=="HOMDEL"] <- -1
selected.data[selected.data=="AMP"] <- 1
selected.data <- as.matrix(selected.data)
selected.data[is.na(selected.data)] <- 0
mode(selected.data) <-"numeric"

#HiLo Gisic
selected.data <- as.data.frame(selected.data)
selected.data$subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(rownames(selected.data),rownames(ClinicalData.subset))]
selected.data$HiLo <- Mut.freq.decile$HiLo[match(rownames(selected.data),rownames(Mut.freq.decile))]
selected.data$Decile <- Mut.freq.decile$decile[match(rownames(selected.data),rownames(Mut.freq.decile))]
selected.data$label<- Mut.freq.decile$label[match(rownames(selected.data),rownames(Mut.freq.decile))]
selected.data$cluster <- Consensus.class$Cluster[match(rownames(selected.data),rownames(Consensus.class))]
selected.data <- selected.data[order(factor(selected.data$cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                     factor(selected.data$subtype,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                     selected.data$Decile
                                     ),]
selected.data.mean.byIMS.byHilo <- ddply(selected.data[,-c(9,10)],.(cluster,subtype,HiLo),colwise(mean))
selected.data.mean.byIMS.byHilo <- selected.data.mean.byIMS.byHilo[order(factor(selected.data.mean.byIMS.byHilo$cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                                         factor(selected.data.mean.byIMS.byHilo$subtype,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                                                         factor(selected.data.mean.byIMS.byHilo$HiLo   ,levels = c("Low","Medium","High"))
                                                                         ),]
selected.data.mean.byIMS.byHilo <- selected.data.mean.byIMS.byHilo[selected.data.mean.byIMS.byHilo$subtype != "Normal-like",] #remove normal like for meansplot

gistic.Meanby.IMS.Hilo <- as.matrix(selected.data.mean.byIMS.byHilo[,4:ncol(selected.data.mean.byIMS.byHilo)])
gistic.selected.data <- as.matrix(selected.data[,1:6])
mode(gistic.Meanby.IMS.Hilo) <- "numeric"
mode(gistic.selected.data) <- "numeric"
write.csv(selected.data, "./3 ANALISYS/selected.csv")
write.csv(selected.data.mean.byIMS.byHilo,"./3 ANALISYS/selected.mean.csv")

gistic.Meanby.IMS.Hilo.colors <- selected.data.mean.byIMS.byHilo[,1:3]
levels(gistic.Meanby.IMS.Hilo.colors$cluster) <- c(levels(gistic.Meanby.IMS.Hilo.colors$cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))
levels (gistic.Meanby.IMS.Hilo.colors$subtype) <- c(levels (gistic.Meanby.IMS.Hilo.colors$subtype),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
levels (gistic.Meanby.IMS.Hilo.colors$HiLo) <- c(levels (gistic.Meanby.IMS.Hilo.colors$HiLo),c("green","yellow","red"))
gistic.Meanby.IMS.Hilo.colors$HiLo[gistic.Meanby.IMS.Hilo.colors$HiLo=="Low"]     <- "green"    
gistic.Meanby.IMS.Hilo.colors$HiLo[gistic.Meanby.IMS.Hilo.colors$HiLo=="Medium"]  <- "yellow"     
gistic.Meanby.IMS.Hilo.colors$HiLo[gistic.Meanby.IMS.Hilo.colors$HiLo=="High"]    <- "red"
gistic.Meanby.IMS.Hilo.colors$subtype[gistic.Meanby.IMS.Hilo.colors$subtype=="Luminal A"]     <- "#eaff00"    
gistic.Meanby.IMS.Hilo.colors$subtype[gistic.Meanby.IMS.Hilo.colors$subtype=="Luminal B"]     <- "#00c0ff"     
gistic.Meanby.IMS.Hilo.colors$subtype[gistic.Meanby.IMS.Hilo.colors$subtype=="Basal-like"]    <- "#da70d6"   
gistic.Meanby.IMS.Hilo.colors$subtype[gistic.Meanby.IMS.Hilo.colors$subtype=="HER2-enriched"] <- "#daa520"            
gistic.Meanby.IMS.Hilo.colors$subtype[gistic.Meanby.IMS.Hilo.colors$subtype=="Normal-like"]   <- "#d3d3d3"
gistic.Meanby.IMS.Hilo.colors$cluster[gistic.Meanby.IMS.Hilo.colors$cluster=="ICR4"]    <- "#FF0000"    
gistic.Meanby.IMS.Hilo.colors$cluster[gistic.Meanby.IMS.Hilo.colors$cluster=="ICR3"]    <- "#FFA500"     
gistic.Meanby.IMS.Hilo.colors$cluster[gistic.Meanby.IMS.Hilo.colors$cluster=="ICR2"]    <- "#00FF00"   
gistic.Meanby.IMS.Hilo.colors$cluster[gistic.Meanby.IMS.Hilo.colors$cluster=="ICR1"]    <- "#0000FF"

gistic.colors <- selected.data[,c(11,7,10)]
levels (gistic.colors$cluster) <- c(levels(gistic.colors$cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))
levels (gistic.colors$subtype) <- c(levels (gistic.colors$subtype),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
gistic.colors$subtype[gistic.colors$subtype=="Luminal A"]     <- "#eaff00"    
gistic.colors$subtype[gistic.colors$subtype=="Luminal B"]     <- "#00c0ff"     
gistic.colors$subtype[gistic.colors$subtype=="Basal-like"]    <- "#da70d6"   
gistic.colors$subtype[gistic.colors$subtype=="HER2-enriched"] <- "#daa520"            
gistic.colors$subtype[gistic.colors$subtype=="Normal-like"]   <- "#d3d3d3"
gistic.colors$cluster[gistic.colors$cluster=="ICR4"]    <- "#FF0000"    
gistic.colors$cluster[gistic.colors$cluster=="ICR3"]    <- "#FFA500"     
gistic.colors$cluster[gistic.colors$cluster=="ICR2"]    <- "#00FF00"   
gistic.colors$cluster[gistic.colors$cluster=="ICR1"]    <- "#0000FF"

my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 3)
gistic.colors <- as.matrix(t(gistic.colors))
colnames(gistic.colors)<-NULL
png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".gistic.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".mutual_exclusion.png"),res=600,height=plot.height,width=plot.width/2,unit="in")     # set filename
heatmap.3(gistic.selected.data,
          main = "HeatMap-MutatedGenes",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=gistic.colors,                         # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="row", 
          density.info="none",
          trace="none",
          labCol=gistic.genes,
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


png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".gistic.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".Mean.HiLo.png"),res=600,height=plot.height,width=plot.width/2,unit="in")     # set filename
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
gistic.Meanby.IMS.Hilo.colors <- as.matrix(t(gistic.Meanby.IMS.Hilo.colors))
colnames(gistic.Meanby.IMS.Hilo.colors) <-NULL
par(mar=c(1,1,1,1))
heatmap.3(gistic.Meanby.IMS.Hilo,
          main = "HeatMap-MutatedGenes frequency",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          #breaks=my.colors,                                   # set manual color gradient of color scheme
          RowSideColors=gistic.Meanby.IMS.Hilo.colors,                        # set goup colors
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          #scale="col", 
          #density.info="none",
          trace="none",
          labCol=colnames(gistic.Meanby.IMS.Hilo),
          cexRow=1,cexCol=2.6,
          margins=c(15,15),
          labRow=FALSE,
          Colv=FALSE,Rowv=FALSE                               # reorder row/columns by dendogram
)
par(lend = 1)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1","","Low (1-3)","Medium (4-7)","High (8-10)"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue","white","green","yellow","red"),lty= 1,lwd = 5,cex = 0.9)
dev.off()


# Heatmap for matching RNAseq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
seleceted.genes.RNAseq <- t(RNASeq.NORM_Log2[selected.genes,])
seleceted.genes.RNAseq <- seleceted.genes.RNAseq[rownames(allmuts.mutatedgenes),]

my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))

png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".RNAseq.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".png"),res=600,height=plot.height,plot.width/2,unit="in")     # set filename
heatmap.3(seleceted.genes.RNAseq[,c(2,3,4,5,1,6)],
          main = paste0("Heatmap RNASeq - ",Geneset),
          col=my.palette,                   #set color sheme RED High, GREEN low
          #breaks=my.colors,                                 
          RowSideColors=color.matrix,      #set goup colors                 
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol = colnames(seleceted.genes.RNAseq[,c(2,3,4,5,1,6)]),labRow = FALSE,
          cexRow=1.3,cexCol=2,margins=c(8,7),
          Colv=FALSE,Rowv=FALSE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

#RNASEQ HiLo edition
RNAseq.table <- as.data.frame(seleceted.genes.RNAseq)
RNAseq.table$subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(rownames(RNAseq.table),rownames(ClinicalData.subset))]
RNAseq.table$HiLo <- Mut.freq.decile$HiLo[match(rownames(RNAseq.table),rownames(Mut.freq.decile))]
RNAseq.table$Decile <- Mut.freq.decile$decile[match(rownames(RNAseq.table),rownames(Mut.freq.decile))]
RNAseq.table$label<- Mut.freq.decile$label[match(rownames(RNAseq.table),rownames(Mut.freq.decile))]
RNAseq.table$cluster <- Consensus.class$Cluster[match(rownames(RNAseq.table),rownames(Consensus.class))]
RNAseq.table <- RNAseq.table[order(factor(RNAseq.table$cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                     factor(RNAseq.table$subtype,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                     RNAseq.table$Decile
),]
RNAseq.table.mean.byIMS.byHilo <- ddply(RNAseq.table[,-c(9,10)],.(cluster,subtype,HiLo),colwise(mean))
RNAseq.table.mean.byIMS.byHilo <- RNAseq.table.mean.byIMS.byHilo[order(factor(RNAseq.table.mean.byIMS.byHilo$cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                                                         factor(RNAseq.table.mean.byIMS.byHilo$subtype,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like")),
                                                                         factor(RNAseq.table.mean.byIMS.byHilo$HiLo   ,levels = c("Low","Medium","High"))
),]
RNAseq.table.mean.byIMS.byHilo <- RNAseq.table.mean.byIMS.byHilo[RNAseq.table.mean.byIMS.byHilo$subtype != "Normal-like",] #remove normal like for meansplot

RNAseq.Meanby.IMS.Hilo <- as.matrix(RNAseq.table.mean.byIMS.byHilo[,4:ncol(RNAseq.table.mean.byIMS.byHilo)])
#RNAseq.table <- as.matrix(RNAseq.table[,1:6])
mode(RNAseq.Meanby.IMS.Hilo) <- "numeric"
#mode(RNAseq.table) <- "numeric"

png(paste0("./4 FIGURES/Heatmaps/mutations/",Cancerset,".",IMS.filter,".",Geneset,".RNASEQ.HeatMap.",matrix.type,".",plot.type,".",cluster.select,".gene.filter_",gene.filter,".png"),res=600,height=plot.height,width=plot.width/2,unit="in")     # set filename
heatmap.3(RNAseq.Meanby.IMS.Hilo[,c(2,3,4,5,1,6)],
          main = "HeatMap-RNASeq-HiLo",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          RowSideColors=gistic.Meanby.IMS.Hilo.colors,                         # set goup colors
          key=FALSE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row", 
          density.info="none",
          trace="none",
          labCol=colnames(RNAseq.Meanby.IMS.Hilo[,c(2,3,4,5,1,6)]),
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

### CNV MAPKs stacked bar charts

##filter luminal
#selected.data <- selected.data[selected.data$subtype %in% c("Luminal A","Luminal B"),]

ICR.count <- count(selected.data,vars = c("cluster"))

MAP3K1.count <- count(selected.data,vars = c("cluster","MAP3K1"))
MAP3K1.count$N <- ICR.count$freq[match(MAP3K1.count$cluster,ICR.count$cluster)]
MAP3K1.count$pct <- MAP3K1.count$freq/MAP3K1.count$N*100
MAP3K1.count[MAP3K1.count$MAP3K1 == 1,][,"MAP3K1"] <- "Amplified"
MAP3K1.count[MAP3K1.count$MAP3K1 == -1,][,"MAP3K1"] <- "Deleted"

MAP2K4.count <- count(selected.data,vars = c("cluster","MAP2K4"))
MAP2K4.count$N <- ICR.count$freq[match(MAP2K4.count$cluster,ICR.count$cluster)]
MAP2K4.count$pct <- MAP2K4.count$freq/MAP2K4.count$N*100
MAP2K4.count[MAP2K4.count$MAP2K4 == 1,][,"MAP2K4"] <- "Amplified"
MAP2K4.count[MAP2K4.count$MAP2K4 == -1,][,"MAP2K4"] <- "Deleted"

ICR1.DEL.count <- MAP3K1.count[MAP3K1.count$cluster=="ICR1" & MAP3K1.count$MAP3K1 =="Deleted","freq"]
ICR1.NOTDEL.count <- sum(MAP3K1.count[MAP3K1.count$cluster=="ICR1" & MAP3K1.count$MAP3K1 !="Deleted","freq"])
ICR4.DEL.count <- MAP3K1.count[MAP3K1.count$cluster=="ICR4" & MAP3K1.count$MAP3K1 =="Deleted","freq"]
ICR4.NOTDEL.count <- sum(MAP3K1.count[MAP3K1.count$cluster=="ICR4" & MAP3K1.count$MAP3K1 !="Deleted","freq"])
chi.matrix.ICR.MAP3K1.DEL <- matrix(c(ICR1.DEL.count,ICR1.NOTDEL.count,ICR4.DEL.count,ICR4.NOTDEL.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix.ICR.MAP3K1.DEL)
chi.matrix.ICR.MAP3K1.DEL.p <- test$p.value

ICR1.DEL.count <- MAP2K4.count[MAP2K4.count$cluster=="ICR1" & MAP2K4.count$MAP2K4 =="Deleted","freq"]
ICR1.NOTDEL.count <- sum(MAP2K4.count[MAP2K4.count$cluster=="ICR1" & MAP2K4.count$MAP2K4 !="Deleted","freq"])
ICR4.DEL.count <- MAP2K4.count[MAP2K4.count$cluster=="ICR4" & MAP2K4.count$MAP2K4 =="Deleted","freq"]
ICR4.NOTDEL.count <- sum(MAP2K4.count[MAP2K4.count$cluster=="ICR4" & MAP2K4.count$MAP2K4 !="Deleted","freq"])
chi.matrix.ICR.MAP2K4.DEL <- matrix(c(ICR1.DEL.count,ICR1.NOTDEL.count,ICR4.DEL.count,ICR4.NOTDEL.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix.ICR.MAP2K4.DEL)
chi.matrix.ICR.MAP2K4.DEL.p <- test$p.value

ICR1.AMP.count <- MAP3K1.count[MAP3K1.count$cluster=="ICR1" & MAP3K1.count$MAP3K1 =="Amplified","freq"]
ICR1.NOTAMP.count <- sum(MAP3K1.count[MAP3K1.count$cluster=="ICR1" & MAP3K1.count$MAP3K1 !="Amplified","freq"])
ICR4.AMP.count <- MAP3K1.count[MAP3K1.count$cluster=="ICR4" & MAP3K1.count$MAP3K1 =="Amplified","freq"]
ICR4.NOTAMP.count <- sum(MAP3K1.count[MAP3K1.count$cluster=="ICR4" & MAP3K1.count$MAP3K1 !="Amplified","freq"])
chi.matrix.ICR.MAP3K1.AMP <- matrix(c(ICR1.AMP.count,ICR1.NOTAMP.count,ICR4.AMP.count,ICR4.NOTAMP.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix.ICR.MAP3K1.AMP)
chi.matrix.ICR.MAP3K1.AMP.p <- test$p.value

ICR1.AMP.count <- MAP2K4.count[MAP2K4.count$cluster=="ICR1" & MAP2K4.count$MAP2K4 =="Amplified","freq"]
ICR1.NOTAMP.count <- sum(MAP2K4.count[MAP2K4.count$cluster=="ICR1" & MAP2K4.count$MAP2K4 !="Amplified","freq"])
ICR4.AMP.count <- MAP2K4.count[MAP2K4.count$cluster=="ICR4" & MAP2K4.count$MAP2K4 =="Amplified","freq"]
ICR4.NOTAMP.count <- sum(MAP2K4.count[MAP2K4.count$cluster=="ICR4" & MAP2K4.count$MAP2K4 !="Amplified","freq"])
chi.matrix.ICR.MAP2K4.AMP <- matrix(c(ICR1.AMP.count,ICR1.NOTAMP.count,ICR4.AMP.count,ICR4.NOTAMP.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix.ICR.MAP2K4.AMP)
chi.matrix.ICR.MAP2K4.AMP.p <- test$p.value

gg = ggplot(MAP3K1.count[MAP3K1.count$MAP3K1!=0,],aes(x=cluster,y=pct,fill=MAP3K1))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("red","blue")) +
  ggtitle("MAP3K1 CNV") +
  ylab("Copy Number Variation in percent") +
  theme_bw()
print(gg)

gg = ggplot(MAP2K4.count[MAP2K4.count$MAP2K4!=0,],aes(x=cluster,y=pct,fill=MAP2K4))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("red","blue")) +
  ggtitle("MAP2K4 CNV") +
  ylab("Copy Number Variation in percent") +
  theme_bw()
print(gg)

#CNV in MUTvsWT
mutstat.data <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
mutstat.data$X  <-NULL
selected.data$MUTSTAT <- mutstat.data$MAP2K4.MAP3K1[match(rownames(selected.data),mutstat.data$Sample)]
MUT.count <- count(selected.data,vars = c("MUTSTAT"))
MUT.count.MAP3K1 <- count(selected.data,vars = c("MUTSTAT","MAP3K1"))
MUT.count.MAP2K4 <- count(selected.data,vars = c("MUTSTAT","MAP2K4"))
MUT.count.MAP3K1$N <- MUT.count$freq[match(MUT.count.MAP3K1$MUTSTAT,MUT.count$MUTSTAT)]
MUT.count.MAP2K4$N <- MUT.count$freq[match(MUT.count.MAP2K4$MUTSTAT,MUT.count$MUTSTAT)]
MUT.count.MAP3K1$pct <- MUT.count.MAP3K1$freq/MUT.count.MAP3K1$N*100
MUT.count.MAP2K4$pct <- MUT.count.MAP2K4$freq/MUT.count.MAP2K4$N*100
MUT.count.MAP3K1[MUT.count.MAP3K1$MAP3K1 == 1,][,"MAP3K1"] <- "Amplified"
MUT.count.MAP3K1[MUT.count.MAP3K1$MAP3K1 == -1,][,"MAP3K1"] <- "Deleted"
MUT.count.MAP2K4[MUT.count.MAP2K4$MAP2K4 == 1,][,"MAP2K4"] <- "Amplified"
MUT.count.MAP2K4[MUT.count.MAP2K4$MAP2K4 == -1,][,"MAP2K4"] <- "Deleted"

MUT.DEL.count <- MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="MUT" & MUT.count.MAP3K1$MAP3K1 =="Deleted","freq"]
MUT.NOTDEL.count <- sum(MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="MUT" & MUT.count.MAP3K1$MAP3K1 !="Deleted","freq"])
WT.DEL.count <- MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="WT" & MUT.count.MAP3K1$MAP3K1 =="Deleted","freq"]
WT.NOTDEL.count <- sum(MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="WT" & MUT.count.MAP3K1$MAP3K1 !="Deleted","freq"])
chi.matrix <- matrix(c(MUT.DEL.count,MUT.NOTDEL.count,WT.DEL.count,WT.NOTDEL.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix)
test$p.value

MUT.DEL.count <- MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="MUT" & MUT.count.MAP2K4$MAP2K4 =="Deleted","freq"]
MUT.NOTDEL.count <- sum(MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="MUT" & MUT.count.MAP2K4$MAP2K4 !="Deleted","freq"])
WT.DEL.count <- MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="WT" & MUT.count.MAP2K4$MAP2K4 =="Deleted","freq"]
WT.NOTDEL.count <- sum(MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="WT" & MUT.count.MAP2K4$MAP2K4 !="Deleted","freq"])
chi.matrix <- matrix(c(MUT.DEL.count,MUT.NOTDEL.count,WT.DEL.count,WT.NOTDEL.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix)
test$p.value

MUT.AMP.count <- MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="MUT" & MUT.count.MAP3K1$MAP3K1 =="Amplified","freq"]
MUT.NOTAMP.count <- sum(MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="MUT" & MUT.count.MAP3K1$MAP3K1 !="AMPeted","freq"])
WT.AMP.count <- MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="WT" & MUT.count.MAP3K1$MAP3K1 =="Amplified","freq"]
WT.NOTAMP.count <- sum(MUT.count.MAP3K1[MUT.count.MAP3K1$MUTSTAT=="WT" & MUT.count.MAP3K1$MAP3K1 !="AMPeted","freq"])
chi.matrix <- matrix(c(MUT.AMP.count,MUT.NOTAMP.count,WT.AMP.count,WT.NOTAMP.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix)
test$p.value

MUT.AMP.count <- MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="MUT" & MUT.count.MAP2K4$MAP2K4 =="Amplified","freq"]
MUT.NOTAMP.count <- sum(MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="MUT" & MUT.count.MAP2K4$MAP2K4 !="Amplified","freq"])
WT.AMP.count <- MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="WT" & MUT.count.MAP2K4$MAP2K4 =="Amplified","freq"]
WT.NOTAMP.count <- sum(MUT.count.MAP2K4[MUT.count.MAP2K4$MUTSTAT=="WT" & MUT.count.MAP2K4$MAP2K4 !="Amplified","freq"])
chi.matrix <- matrix(c(MUT.AMP.count,MUT.NOTAMP.count,WT.AMP.count,WT.NOTAMP.count),nrow = 2, ncol = 2)
test<-chisq.test(chi.matrix)
test$p.value

gg = ggplot(MUT.count.MAP3K1[MUT.count.MAP3K1$MAP3K1!=0,],aes(x=MUTSTAT,y=pct,fill=MAP3K1))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("red","blue")) +
  ggtitle("MAP3K1 CNV") +
  ylab("Copy Number Variation in percent") +
  theme_bw()
print(gg)

gg = ggplot(MUT.count.MAP2K4[MUT.count.MAP2K4$MAP2K4!=0,],aes(x=MUTSTAT,y=pct,fill=MAP2K4))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("red","blue")) +
  ggtitle("MAP2K4 CNV") +
  ylab("Copy Number Variation in percent") +
  theme_bw()
print(gg)

#CNV in Subtype
IMS.count <-count(selected.data,vars=c("subtype"))
IMS.count.MAP3K1 <- count(selected.data,vars = c("subtype","MAP3K1"))
IMS.count.MAP2K4 <- count(selected.data,vars = c("subtype","MAP2K4"))
IMS.count.MAP3K1$N <- IMS.count$freq[match(IMS.count.MAP3K1$subtype,IMS.count$subtype)]
IMS.count.MAP2K4$N <- IMS.count$freq[match(IMS.count.MAP2K4$subtype,IMS.count$subtype)]
IMS.count.MAP3K1$pct <- IMS.count.MAP3K1$freq/IMS.count.MAP3K1$N*100
IMS.count.MAP2K4$pct <- IMS.count.MAP2K4$freq/IMS.count.MAP2K4$N*100
IMS.count.MAP3K1[IMS.count.MAP3K1$MAP3K1 == 1,][,"MAP3K1"] <- "Amplified"
IMS.count.MAP3K1[IMS.count.MAP3K1$MAP3K1 == -1,][,"MAP3K1"] <- "Deleted"
IMS.count.MAP2K4[IMS.count.MAP2K4$MAP2K4 == 1,][,"MAP2K4"] <- "Amplified"
IMS.count.MAP2K4[IMS.count.MAP2K4$MAP2K4 == -1,][,"MAP2K4"] <- "Deleted"
IMS.count.MAP3K1 <- IMS.count.MAP3K1[order(factor(IMS.count.MAP3K1$subtype,levels = c("Luminal A","Luminal B","HER2-enriched","Basal-like"))),]
IMS.count.MAP2K4 <- IMS.count.MAP2K4[order(factor(IMS.count.MAP2K4$subtype,levels = c("Luminal A","Luminal B","HER2-enriched","Basal-like"))),]
IMS.count.MAP3K1$subtype <- factor(IMS.count.MAP3K1$subtype, levels=unique(as.character(IMS.count.MAP3K1$subtype)))
IMS.count.MAP2K4$subtype <- factor(IMS.count.MAP2K4$subtype, levels=unique(as.character(IMS.count.MAP2K4$subtype)))
IMS.count.MAP3K1 <- IMS.count.MAP3K1[IMS.count.MAP3K1$subtype != "Normal-like",]
IMS.count.MAP2K4 <- IMS.count.MAP2K4[IMS.count.MAP2K4$subtype != "Normal-like",]

gg = ggplot(IMS.count.MAP3K1[IMS.count.MAP3K1$MAP3K1!=0,],aes(x=subtype,y=pct,fill=MAP3K1))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("red","blue")) +
  ggtitle("MAP3K1 CNV") +
  ylab("Copy Number Variation in percent") +
  theme_bw()
print(gg)

gg = ggplot(IMS.count.MAP2K4[IMS.count.MAP2K4$MAP2K4!=0,],aes(x=subtype,y=pct,fill=MAP2K4))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("red","blue")) +
  ggtitle("MAP2K4 CNV") +
  ylab("Copy Number Variation in percent") +
  theme_bw()
print(gg)

#expression vs CNV
selected.data$MAP3K1.exp <- RNAseq.table$MAP3K1[match(rownames(selected.data),rownames(RNAseq.table))]
selected.data$MAP2K4.exp <- RNAseq.table$MAP2K4[match(rownames(selected.data),rownames(RNAseq.table))]

ICR.MAP3K1.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","MAP3K1.exp")],mean)
ICR.MAP2K4.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","MAP2K4.exp")],mean)
MUTSTAT.MAP3K1.exp <- aggregate(.~MUTSTAT,data=selected.data[,c("MUTSTAT","MAP3K1.exp")],mean)
MUTSTAT.MAP2K4.exp <- aggregate(.~MUTSTAT,data=selected.data[,c("MUTSTAT","MAP2K4.exp")],mean)
CNV.MAP3K1.exp <- aggregate(.~MAP3K1,data=selected.data[,c("MAP3K1","MAP3K1.exp")],mean)
CNV.MAP2K4.exp <- aggregate(.~MAP2K4,data=selected.data[,c("MAP2K4","MAP2K4.exp")],mean)

CNV.MAP3K1.exp[CNV.MAP3K1.exp$MAP3K1 == 1,][,"MAP3K1"] <- "Amplified"
CNV.MAP3K1.exp[CNV.MAP3K1.exp$MAP3K1 == -1,][,"MAP3K1"] <- "Deleted"
CNV.MAP3K1.exp[CNV.MAP3K1.exp$MAP3K1 == 0,][,"MAP3K1"] <- "Normal"
CNV.MAP2K4.exp[CNV.MAP2K4.exp$MAP2K4 == 1,][,"MAP2K4"] <- "Amplified"
CNV.MAP2K4.exp[CNV.MAP2K4.exp$MAP2K4 == -1,][,"MAP2K4"] <- "Deleted"
CNV.MAP2K4.exp[CNV.MAP2K4.exp$MAP2K4 == 0,][,"MAP2K4"] <- "Normal"

T1<-cbind("MAP3K1.exp",ICR.MAP3K1.exp)
T2<-cbind("MAP2K4.exp",ICR.MAP2K4.exp)
T3<-cbind("MAP3K1.exp",MUTSTAT.MAP3K1.exp)
T4<-cbind("MAP2K4.exp",MUTSTAT.MAP2K4.exp)
T5<-cbind("MAP3K1.exp",CNV.MAP3K1.exp)
T6<-cbind("MAP2K4.exp",CNV.MAP2K4.exp)

colnames(T1) = c("MAPK","Group","Expression")
colnames(T2) = c("MAPK","Group","Expression")
colnames(T3) = c("MAPK","Group","Expression")
colnames(T4) = c("MAPK","Group","Expression")
colnames(T5) = c("MAPK","Group","Expression")
colnames(T6) = c("MAPK","Group","Expression")

expression.table<-rbind (T1,T2,T3,T4,T5,T6)

gg<-ggplot(expression.table,aes(x=Group,y=Expression,fill=MAPK))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("MAPK expression") +
  ylab("Expression (Log2)") +
  theme_bw()
print(gg)
  