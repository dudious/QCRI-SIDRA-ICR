# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

#load data
MAPK.PW.genes <- read.csv ("./3 ANALISYS/MAPK/MAPK.pathway.genes.csv",header = FALSE)

load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGinICR4vs1.RDATA")
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGsLM_ICR4vs1.rdata")
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.rdata")

#select significant MAPK.PW genes in TCGA ICR DEG  or MAPK DEG

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

###### UPREGULATED Mutated in MAP3K1 or MAP2K4
DEGinMUTvsWT.lum.MAPK.PW<-DEGs[rownames(DEGs)%in%MAPK.PW.genes$V2,]
DEGinMUTvsWT.lum.MAPK.PW<-DEGinMUTvsWT.lum.MAPK.PW[DEGinMUTvsWT.lum.MAPK.PW$logFC>0,]
DEGinMUTvsWT.lum.MAPK.PW<-DEGinMUTvsWT.lum.MAPK.PW[order(DEGinMUTvsWT.lum.MAPK.PW$FDR),]
DEGinMUTvsWT.lum.MAPK.PW<-DEGinMUTvsWT.lum.MAPK.PW[DEGinMUTvsWT.lum.MAPK.PW$FDR<0.05,]

#load extra data
#RNASeq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
RNASeq.subset <- RNASeq.NORM_Log2[rownames(DEGinMUTvsWT.lum.MAPK.PW),]

#MAPKMUT data
MAPKMUT <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
rownames(MAPKMUT) <- MAPKMUT$Sample
MAPKMUT$X <- NULL
MAPKMUT$Sample <- NULL

#clinical data
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL

#clustering data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#lookup MUTstat
mutstat <- MAPKMUT[colnames(RNASeq.subset),"MAP2K4.MAP3K1",drop=FALSE]
levels (mutstat$MAP2K4.MAP3K1) <- c(levels (mutstat$MAP2K4.MAP3K1),c("#8b0000","grey"))
mutstat<-mutstat[!is.na(mutstat$MAP2K4.MAP3K1),,drop=F]
mutstat$MAP2K4.MAP3K1[mutstat$MAP2K4.MAP3K1=="WT"] <- "grey"
mutstat$MAP2K4.MAP3K1[mutstat$MAP2K4.MAP3K1=="MUT"] <- "#8b0000"
#drop data without mutation data and sort expression data
RNASeq.subset<-RNASeq.subset[,rownames(mutstat)]
means <- colMeans(RNASeq.subset)
RNASeq.subset <- RNASeq.subset[,order(means)]
#reorder mutstat
mutstat<- mutstat[colnames(RNASeq.subset),,drop=FALSE]

#lookup cluster
cluster <- Consensus.class[colnames(RNASeq.subset),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#lookup subtype
RNASeq.order <- colnames(RNASeq.subset)[which(colnames(RNASeq.subset) %in% rownames(ClinicalData.subset))]
subtype <- ClinicalData.subset[RNASeq.order,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
subtype <-subtype [!is.na(subtype$TCGA.PAM50.RMethod.RNASeq),,drop=F]
levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"]     <- "#eaff00"    
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"]     <- "#00c0ff"     
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"]    <- "#da70d6"   
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#daa520"            
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Normal-like"]   <- "#d3d3d3"

#cluster filter
cluster <- cluster[cluster$Cluster=="#FF0000"|cluster$Cluster=="#0000FF",,drop=FALSE]
subtype <- subtype[rownames(cluster),,drop=FALSE]
mutstat <- mutstat[rownames(cluster),,drop=FALSE]
RNASeq.subset <- RNASeq.subset[,rownames(cluster)]

#subtype filter
subtype <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="#eaff00"|subtype$TCGA.PAM50.RMethod.RNASeq=="#00c0ff",,drop=FALSE]
cluster <- cluster[rownames(subtype),,drop=FALSE]
mutstat <- mutstat[rownames(subtype),,drop=FALSE]
RNASeq.subset <- RNASeq.subset[,rownames(subtype)]

#color matrix
mutcolors <- as.character(mutstat$MAP2K4.MAP3K1)
clustercolors <- as.character(cluster$Cluster)
subtypecolors <- as.character(subtype$TCGA.PAM50.RMethod.RNASeq)
color.matrix <- as.matrix(rbind (clustercolors,subtypecolors,mutcolors))
rownames(color.matrix) <- c("Cluster","Subtype","MAPK-mutation")

#sample order
sample.order <- colnames(RNASeq.subset)

#heatmap
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))


png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.TCGA.SIGN.UP.MUTvsWT.luminal.png"),res=600,height=6,width=6,unit="in")     # set filename
heatmap.3(RNASeq.subset,
          main = "SIGN.UP.MUTvsWT.luminal",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          ColSideColors=t(color.matrix),                         # set goup colors
          breaks = my.colors,
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1.2,cexCol=0.1,margins=c(9,9),
          Colv=FALSE,Rowv=TRUE                           # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1","","Mutated","Wild Type"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue","white","#8b0000","grey"),lty= 1,lwd = 5,cex = 0.3)
dev.off()
dim(RNASeq.subset)

#######DOWNREGULATED Mutated in MAP3K1 or MAP2K4
DEGinMUTvsWT.lum.MAPK.PW<-DEGs[rownames(DEGs)%in%MAPK.PW.genes$V2,]
DEGinMUTvsWT.lum.MAPK.PW<-DEGinMUTvsWT.lum.MAPK.PW[DEGinMUTvsWT.lum.MAPK.PW$logFC<0,]
DEGinMUTvsWT.lum.MAPK.PW<-DEGinMUTvsWT.lum.MAPK.PW[order(DEGinMUTvsWT.lum.MAPK.PW$FDR),]
DEGinMUTvsWT.lum.MAPK.PW<-DEGinMUTvsWT.lum.MAPK.PW[DEGinMUTvsWT.lum.MAPK.PW$FDR<0.05,]

#load extra data
#RNASeq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
RNASeq.subset <- RNASeq.NORM_Log2[rownames(DEGinMUTvsWT.lum.MAPK.PW),]

#reverse expression
RNASeq.subset <- -(RNASeq.subset)

#MAPKMUT data
MAPKMUT <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
rownames(MAPKMUT) <- MAPKMUT$Sample
MAPKMUT$X <- NULL
MAPKMUT$Sample <- NULL

#clinical data
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL

#clustering data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#lookup MUTstat
mutstat <- MAPKMUT[colnames(RNASeq.subset),"MAP2K4.MAP3K1",drop=FALSE]
levels (mutstat$MAP2K4.MAP3K1) <- c(levels (mutstat$MAP2K4.MAP3K1),c("#8b0000","grey"))
mutstat<-mutstat[!is.na(mutstat$MAP2K4.MAP3K1),,drop=F]
mutstat$MAP2K4.MAP3K1[mutstat$MAP2K4.MAP3K1=="WT"] <- "grey"
mutstat$MAP2K4.MAP3K1[mutstat$MAP2K4.MAP3K1=="MUT"] <- "#8b0000"
#drop data without mutation data and sort expression data
RNASeq.subset<-RNASeq.subset[,rownames(mutstat)]
means <- colMeans(RNASeq.subset)
RNASeq.subset <- RNASeq.subset[,order(means)]
#reorder mutstat
mutstat<- mutstat[colnames(RNASeq.subset),,drop=FALSE]

#lookup cluster
cluster <- Consensus.class[colnames(RNASeq.subset),"Cluster",drop=FALSE]
levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"

#lookup subtype
RNASeq.order <- colnames(RNASeq.subset)[which(colnames(RNASeq.subset) %in% rownames(ClinicalData.subset))]
subtype <- ClinicalData.subset[RNASeq.order,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
subtype <-subtype [!is.na(subtype$TCGA.PAM50.RMethod.RNASeq),,drop=F]
levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"]     <- "#eaff00"    
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"]     <- "#00c0ff"     
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"]    <- "#da70d6"   
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#daa520"            
subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Normal-like"]   <- "#d3d3d3"

#cluster filter
cluster <- cluster[cluster$Cluster=="#FF0000"|cluster$Cluster=="#0000FF",,drop=FALSE]
subtype <- subtype[rownames(cluster),,drop=FALSE]
mutstat <- mutstat[rownames(cluster),,drop=FALSE]
RNASeq.subset <- RNASeq.subset[,rownames(cluster)]

#subtype filter
subtype <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="#eaff00"|subtype$TCGA.PAM50.RMethod.RNASeq=="#00c0ff",,drop=FALSE]
cluster <- cluster[rownames(subtype),,drop=FALSE]
mutstat <- mutstat[rownames(subtype),,drop=FALSE]
RNASeq.subset <- RNASeq.subset[,rownames(subtype)]

#reorder to fit Downregulated
RNASeq.subset <- RNASeq.subset[,sample.order]
subtype <- subtype[sample.order,,drop=FALSE]
cluster <- cluster[sample.order,,drop=FALSE]
mutstat <- mutstat[sample.order,,drop=FALSE]

#color matrix
mutcolors <- as.character(mutstat$MAP2K4.MAP3K1)
clustercolors <- as.character(cluster$Cluster)
subtypecolors <- as.character(subtype$TCGA.PAM50.RMethod.RNASeq)
color.matrix <- as.matrix(rbind (clustercolors,subtypecolors,mutcolors))
rownames(color.matrix) <- c("Cluster","Subtype","MAPK-mutation")

#heatmap
my.palette <- colorRampPalette(c("red", "yellow", "blue"))(n = 297)
my.colors = unique(c(seq(-6,-1,length=100),seq(-1,1,length=100),seq(1,6,length=100)))


png(paste0("./4 FIGURES/Heatmaps/DEG_heatmap.TCGA.SIGN.DOWN.rev.uporder.MUTvsWT.luminal.png"),res=600,height=9,width=6,unit="in")     # set filename
heatmap.3(RNASeq.subset,
          main = "SIGN.DOWN.rev.uporder.MUTvsWT.luminal",
          col=my.palette,                                     # set color scheme RED High, GREEN low
          ColSideColors=t(color.matrix),                         # set goup colors
          breaks = my.colors,
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          cexRow=1.2,cexCol=0.1,margins=c(9,9),
          Colv=FALSE,Rowv=TRUE                           # reorder row/columns by dendogram
)
par(lend = 1)
#legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
#       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1","","Mutated","Wild Type"),
       col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue","white","#8b0000","grey"),lty= 1,lwd = 5,cex = 0.3)
dev.off()
dim(RNASeq.subset)
