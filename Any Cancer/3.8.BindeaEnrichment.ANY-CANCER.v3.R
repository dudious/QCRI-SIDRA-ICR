# Setup environment
rm(list=ls())
#setwd("~/Dropbox/BREAST_QATAR")
#setwd("/mnt3/wouter/BREAST-QATAR/")
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR")
## dependencies
required.packages.BioC <- c("GSVA")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
if(length(missing.packages)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(missing.packages)
}
required.packages <- c("heatmap3")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)

library(GSVA)
require(gplots)
require(heatmap3)

# Set Parameters
Cancersets        = "ALL"
Geneset           = "DBGS3"   
DL.Method         = "BIOLINKS"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Selected"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer
K                 = 4

TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)

#test
Cancerset = "BRCA"

#load data
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".",sample.types,".NORMALIZED.TP_FILTERED_LOG2.RData"))
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000/",
                                   Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".FLTR.reps5000.k=",K,".consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
rownames(Consensus.class) <- Consensus.class$PatientID

load ("./2 DATA/marker_list_bindea2.RData")

#fit variables
normCounts <- RNASeq.NORM.TP_Log2
DatiCli <- Consensus.class

#upload the variability in the attachment

#mettere nello stesso ordine matrice e dati clinici
DatiCli <- DatiCli[colnames(normCounts),]

ES <- gsva(normCounts,marker_list,method="ssgsea")
dim(ES)#
ESz <- ES #come se rinormalizzasse per tutti i campioni
for(i in 1: nrow(ES))
{
  ESz[i,]<- (ES[i,]-mean(ES[i,]))/sd(ES[i,])
}

sHc <- hclust(ddist <- dist(t(ESz)), method = "ward.D2")
plot(sHc,labels=FALSE)
rawCluster<- cutree(sHc,k = 2)

pd2<-DatiCli[names(rawCluster),]

table(pd2[names(rawCluster),"Group"],rawCluster)
ICR.col<-as.character(pd2[,"Group"])
names(ICR.col)<-rownames(pd2)
table(ICR.col)
ICR.col[ICR.col=="ICR1"]<-"blue"
ICR.col[ICR.col=="ICR2"]<-"green"
ICR.col[ICR.col=="ICR3"]<-"orange"
ICR.col[ICR.col=="ICR4"]<-"red"

annot<-as.vector(ICR.col)

#for(cc in unique(subtype.col)) rug(which(subtype.col[sHc$order] == cc), col = cc, lwd = 3, ticksize=-.06)
for(cc in unique(ICR.col)) rug(which(ICR.col[sHc$order] == cc), col = cc, lwd = 3, ticksize=-.06)

#log o nolog= stessi risultati
pdf("./3 ANALISYS/IMMUNOSCORE/BINDEA/HM_ssGSEA_bindea_TCGA_ICR4-1-noLOG.pdf",width = 16,height = 8 )

sHc2<-reorder(as.dendrogram(sHc),c(length(sHc$labels):1))
plot(sHc2)
dev.off()
pdf("./3 ANALISYS/IMMUNOSCORE/BINDEA/bindea.pdf",width = 16,height = 8 )
heatmap3((as.matrix(ESz)),
         main="ssGSEA/bindea signatures",
         ColSideColors=annot,
         Colv= as.dendrogram(sHc),col=bluered(75) , labCol=NA,
         scale='col',margins = c(12, 7))
dev.off()
pdf("./3 ANALISYS/IMMUNOSCORE/BINDEA/bindea_reverse.pdf",width = 16,height = 8 )
heatmap3((as.matrix(ESz)),
         main="ssGSEA/bindea signatures",
         ColSideColors=annot,
         Colv= as.dendrogram(sHc2),col=bluered(75) , labCol=NA,
         scale='col',margins = c(12, 7))
dev.off()
