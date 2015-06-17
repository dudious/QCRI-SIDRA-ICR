#################################################################
###
### This script downloads the Breast cancer copy number Data 
### from the TCGA database
### It will process the data into a  single table. 
### Data is saved :
### ../2 DATA/TCGA TCGA CNV/TCGA CNV_Brca_ASSEMBLER/...
### File to use :
### Brca.CN.TCGA.ASSEMBLER.DATA.hg19.txt
### Brca.GeneLevel.CNA.hg19.txt
###
#################################################################

#Download the CNA?CNV data using TCGA Assembler
# Setup environment
rm(list=ls())
getwd()
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("HGNChelper","RCurl","httr","stringr","digest","bitops")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")


#Download data
start.time <- Sys.time ()
CNVRawData <- DownloadCNAData(traverseResultFile="./2 DATA/DirectoryTraverseResult_May-06-2015.rda",
                              saveFolderName="./2 DATA/TCGA CNV/TCGA CNV_Brca_ASSEMBLER/",
                              cancerType= "Brca",
                              assayPlatform="genome_wide_snp_6",
                              outputFileName="Brca.CN.TCGA.ASSEMBLER.DATA")
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

### manually rename file (remove _Brca_enc.edu.....)
cat ("Copy Number Data downloaded,...
     Please rename the file to : Brca.CN.TCGA.ASSEMBLER.DATA.txt - remove _Brca_enc.edu -,...
     Press ENTER when ready... ")


SegTable <- read.table("./2 DATA/TCGA CNV/TCGA CNV_BRCA_ASSEMBLER/BRCA.CN.TCGA.ASSEMBLER.DATA.nocnv.hg19.txt",as.is=T,header = T)
# dim(SegTable) #181619      6


################################################################################################################


#FOR DAVIDE BEDOGNETTI SIGNATURE
load("./3 ANALISYS/CLUSTERING/RNAseq/BRCA/BRCA.TCGA.EDASeq.k7.DBGS1.reps5000/ConsensusClusterObject.Rdata")
consensusClusters <- as.factor(ConsensusClusterObject[[4]]$clrs[[1]])
names(consensusClusters) <- attr(ddist, "Labels")
consensusClusters
clusters <- levels(consensusClusters)
clusters
dir.create("./3 ANALISYS/GISTIC/GISTIC.BRCA")
dir.create("./3 ANALISYS/GISTIC/GISTIC.BRCA/FILE.FOR.GISTIC")
setwd("./3 ANALISYS/GISTIC/GISTIC.BRCA/FILE.FOR.GISTIC/")
for ( cl in clusters){
  
  samples_cl <- names(consensusClusters[consensusClusters==cl])
  seg_cl <- SegTable[substr(SegTable$Sample,1,12) %in% samples_cl,]
  SegFileName <- paste("Cluster_",cl,"_seg.txt",sep="")
  write.table(seg_cl,file=SegFileName,row.names=FALSE,sep="\t",quote=FALSE)
}

setwd("~/Dropbox/BREAST_QATAR")
ICRCLA<-read.delim(file="./3 ANALISYS/CLUSTERING/RNAseq/BRCA/BRCA.TCGA.EDASeq.k7.DBGS1.reps5000/BRCA.TCGA.EDASeq.k7.DBGS1.reps5000.k=4.consensusClass.ICR.csv",sep=",")
ICRCLA<-ICRCLA[order(ICRCLA[,3]),]
rownames(ICRCLA)<-ICRCLA$X
ICRCLA$PatientID<-NULL
ICRCLA$X<-NULL

Cluster1F78B4<-read.delim(file="./3 ANALISYS/GISTIC/GISTIC.BRCA/FILE.FOR.GISTIC//Cluster_#1F78B4_seg.txt")
tmp<-which(rownames(ICRCLA) %in% substr(Cluster1F78B4$Sample,1,12))#ICR3

Cluster33A02C<-read.delim(file="./3 ANALISYS/GISTIC/GISTIC.BRCA/FILE.FOR.GISTIC//Cluster_#33A02C_seg.txt")
tmp<-which(rownames(ICRCLA) %in% substr(Cluster33A02C$Sample,1,12))#ICR2
# 
ClusterA6CEE3<-read.delim(file="./3 ANALISYS/GISTIC/GISTIC.BRCA/FILE.FOR.GISTIC//Cluster_#A6CEE3_seg.txt")
tmp<-which(rownames(ICRCLA) %in% substr(ClusterA6CEE3$Sample,1,12))#ICR1
# 
ClusterB2DF8A<-read.delim(file="./3 ANALISYS/GISTIC/GISTIC.BRCA/FILE.FOR.GISTIC//Cluster_#B2DF8A_seg.txt")
tmp<-which(rownames(ICRCLA) %in% substr(ClusterB2DF8A$Sample,1,12))#ICR4
#RENAME THE FILE IN THE FOLDER WITH THE ICRCLUSTER FLAG
