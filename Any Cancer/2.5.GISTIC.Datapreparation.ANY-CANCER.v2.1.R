################################################################
###
### This scripts prepares the TCGA CNV data for GISTIC analysis
###
###############################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

# Set Parameters
Cancerset   <- "BRCA"
Geneset     <- "DBGS1"

# Load Data

SegTable <- read.table(paste0("./2 DATA/TCGA CNV/TCGA CNV_",Cancerset,"_ASSEMBLER/",Cancerset,".CN.TCGA.ASSEMBLER.DATA.nocnv.hg19.txt"),as.is=TRUE,header = TRUE)
load(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/ConsensusClusterObject.Rdata"))
consensusClusters <- as.factor(ConsensusClusterObject[[4]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
clusters <- levels(consensusClusters)
dir.create(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset))
dir.create(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC"))
setwd(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/"))
for ( cl in clusters){
  
  samples_cl <- names(consensusClusters[consensusClusters==cl])
  seg_cl <- SegTable[substr(SegTable$Sample,1,12) %in% samples_cl,]
  SegFileName <- paste0(Geneset,".Cluster_",cl,"_seg.txt")
  write.table(seg_cl,file=SegFileName,row.names=FALSE,sep="\t",quote=FALSE)
}

setwd("~/Dropbox/BREAST_QATAR")
ICRCLA<-read.delim(file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",
                               Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",
                               Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),sep=",")
ICRCLA<-ICRCLA[order(ICRCLA[,3]),]
rownames(ICRCLA)<-ICRCLA$X
ICRCLA$PatientID<-NULL
ICRCLA$X<-NULL

Cluster1F78B4<-read.delim(file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_#1F78B4_seg.txt"))
tmp<-which(rownames(ICRCLA) %in% substr(Cluster1F78B4$Sample,1,12))#ICR3

Cluster33A02C<-read.delim(file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_#33A02C_seg.txt"))
tmp<-which(rownames(ICRCLA) %in% substr(Cluster33A02C$Sample,1,12))#ICR2
# 
ClusterA6CEE3<-read.delim(file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_#A6CEE3_seg.txt"))
tmp<-which(rownames(ICRCLA) %in% substr(ClusterA6CEE3$Sample,1,12))#ICR1
#
ClusterB2DF8A<-read.delim(file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_#B2DF8A_seg.txt"))
tmp<-which(rownames(ICRCLA) %in% substr(ClusterB2DF8A$Sample,1,12))#ICR4
#RENAME THE FILE IN THE FOLDER WITH THE ICRCLUSTER FLAG
