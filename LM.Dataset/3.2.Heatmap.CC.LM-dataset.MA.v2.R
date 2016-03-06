#################################################################
###
### This Script PLots Heatmaps based on 
### Consensus Clustering grouping of LM.BCRA MA Data
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/...
### Data is saved :
### ./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/...
### Figures are saved :
### ./4 FIGURES/Heatmaps
###
### Parameters to modified : Geneset, K 
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("gplots","GMD")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library (corrplot)

# Set Parameters
Cancerset <- "LM.Dataset"
Geneset <- "DBGS3" # SET GENESET HERE !!!!!!!!!!!!!!
K <- 4             # SET K here

# Load Data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/",Cancerset,".MA.k7.",Geneset,".reps5000/",Cancerset,".MA.k7.",Geneset,".reps5000.k=4.consensusClass.csv"),header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/",Cancerset,".MA.subset.",Geneset,".RData"))
MA.subset <- as.matrix(MA.subset)
load (paste0("./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/",Cancerset,".MA.k7.",Geneset,".reps5000/ConsensusClusterObject.Rdata"))
load ("./2 DATA/LM.BRCA/LM.Dataset.split.Rdata") 

## Code to reorder within cluster, expression data not used
consensusClusters <- as.factor(ConsensusClusterObject[[K]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
hhc <- ConsensusClusterObject[[K]]$consensusTree
sampleOrder <- consensusClusters[hhc$order]                                            ## get the order of cluser assignments based on the consensus tree
ConsensusClusterObject.oGE <- t(MA.subset[names(sampleOrder),])

#Add cluster assignment to data
MA.subset <- merge (MA.subset,Consensus.class,by="row.names")
row.names(MA.subset) <- MA.subset$Row.names
MA.subset$Row.names <- NULL
MA.subset$PatientID <- NULL

#Rename ICR clusters
MA.subset <- MA.subset[-1783,] # remove ALL-NA sample
Cluster.order <- data.frame(Group=MA.subset[,ncol(MA.subset)], avg=rowMeans (MA.subset[,1:(ncol(MA.subset)-1)]))
Cluster.order <- aggregate(Cluster.order,by=list(Cluster.order$Group),FUN=mean)
Cluster.order <- cbind(Cluster.order[order(Cluster.order$avg),c(2,3)],ICR.name=c("ICR1","ICR2","ICR3","ICR4"))
Consensus.class$Group[Consensus.class$Group==Cluster.order[1,1]] <- as.character(Cluster.order[1,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[2,1]] <- as.character(Cluster.order[2,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[3,1]] <- as.character(Cluster.order[3,3])
Consensus.class$Group[Consensus.class$Group==Cluster.order[4,1]] <- as.character(Cluster.order[4,3])
write.csv (Consensus.class,paste0(file="./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/",Cancerset,".MA.k7.",Geneset,".reps5000/",Cancerset,".MA.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"))       

#Update Cluster names
MA.subset$Group <- NULL
MA.subset <- merge (MA.subset,Consensus.class,by="row.names")
row.names(MA.subset) <- MA.subset$Row.names
MA.subset$Row.names <- NULL
MA.subset$PatientID <- NULL

#ordeing within clusters (comment out if no reordering within cluster is required)
# MA.subset   <- MA.subset[colnames(ConsensusClusterObject.oGE),]

#ordering of the clusters
MA.subset <- MA.subset[order(factor(MA.subset$Group,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
MA.subset$Group <- NULL

#re-order the labels
Consensus.class <- Consensus.class[rownames(MA.subset),]

MA.bygene.matrix <- data.frame(t(MA.subset))
MA.bygene.matrix$Gene <- Gene.Meta.data$Symbol[match(rownames(MA.bygene.matrix),Gene.Meta.data$Affy_Probe_ID)]
#rownames(MA.bygene.matrix) <- MA.bygene.matrix$Gene  # Genbe names are not unique
Genes.names <- paste0(MA.bygene.matrix$Gene,"(",sapply(strsplit(rownames(MA.bygene.matrix),"AFFX-HUMISGF3A/"),tail,1),")")
MA.bygene.matrix.ag <- aggregate(MA.bygene.matrix,list(MA.bygene.matrix$Gene),FUN=mean)
rownames(MA.bygene.matrix.ag) <- MA.bygene.matrix.ag$Group.1
MA.bygene.matrix.ag$Group.1 <- NULL
MA.bygene.matrix.ag$Gene <-NULL
MA.bygene.matrix.ag <- t(MA.bygene.matrix.ag)
MA.bygene.matrix$Gene <-NULL
MA.bygene.matrix <- as.matrix(MA.bygene.matrix)
                              
# Heatmap 2 (simple no extra annotations)
patientcolors <- Consensus.class
levels (patientcolors$Group) <- c(levels (patientcolors$Group),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
patientcolors$Group[patientcolors$Group=="ICR4"] <- "#FF0000"
patientcolors$Group[patientcolors$Group=="ICR3"] <- "#FFA500"
patientcolors$Group[patientcolors$Group=="ICR2"] <- "#00FF00"
patientcolors$Group[patientcolors$Group=="ICR1"] <- "#0000FF"
#patientcolors$Group <- droplevels(patientcolors$Group)
patientcolors <- patientcolors$Group
my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
my.colors = unique(c(seq(-4,-0.5,length=100),seq(-0.5,1,length=100),seq(1,4,length=100)))
png(paste0("./4 FIGURES/Heatmaps/Heatmap.MA.",Cancerset,".",Geneset,".png"),res=600,height=6,width=7,unit="in")     # set filename
heatmap.2(MA.bygene.matrix,
          main = paste0("Heatmap MA - ",Geneset," sel., K=",K),
          col=my.palette,                   #set color sheme red High, Yellow low
          breaks=my.colors,                                 
          ColSideColors=patientcolors,      #set goup colors                 
          key=TRUE,
          symm=FALSE,
          symkey=FALSE,
          symbreaks=TRUE,             
          scale="row",
          density.info="none",
          trace="none",
          labCol=FALSE,
          labRow=Genes.names,
          cexRow=1,cexCol=0.1,margins=c(2,12),
          Colv=FALSE)
par(lend = 1)
legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
dev.off()

immunescore<-data.frame(unscaled.score=colMeans(MA.bygene.matrix))
write.csv(immunescore,"./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BRCA.LMDATA.csv")

# Corelation matrix
MA.subset.cor <- cor (MA.bygene.matrix.ag,method="spearman")

# cor significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
MA.subset.cor.sign <- cor.mtest(MA.subset.cor, 0.95)




# Correlation plot
png(paste0("./4 FIGURES/CORRELATION/",Geneset,"/correlation.",Cancerset,".",Geneset,".png"),res=600,height=6,width=6,unit="in")  #adjust output file names here !!!!!
cex.before <- par("cex")
par(cex = 0.45)
col1 = colorRampPalette(c("blue", "white", "#009900"))
lims=c(-1,1)
if (length(MA.subset.cor[MA.subset.cor<0]) == 0) {lims=c(0,1)} 
corrplot.mixed (MA.subset.cor,
                #type="lower",
                #p.mat = RNASeq.subset.cor.sign[[1]],    # add significance to correlations
                col = col1(100),
                lower = "square",
                upper = "number",
                order="FPC",
                cl.lim=lims,                      # only positive correlations
                tl.pos ="lt",
                tl.col = "#c00000",
                #insig="pch",                          # remove insignificant correlations
                tl.cex = 1/par("cex"),
                cl.cex = 1/par("cex"),
                title = paste0("Spearman LM.MA (",Cancerset,".",Geneset," avg:",round(mean(MA.subset.cor),2),")"),
                cex.main = 1.4/par("cex"),
                mar=c(5.1,4.1,4.1,2.1))
par(cex = cex.before)
dev.off()
