#input CircosPlot
# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

#load BRCA BSF Master File
#ClusterAss<-read.csv("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF/BRCA.BSF.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv", sep=",")
Master<-read.delim("./3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF.RNASeq_subset_DBGS3.FLTR.Master.csv",sep=",")
colnames(Master)# to select the col of interst
subMaster<- as.matrix(Master[,c(1,15,30,34)])
#write only Stage I II III IV
unique(subMaster[,2])
tmp<-which(subMaster[,2]=="Stage IA")#84
subMaster[tmp,2]<-"Stage I"
tmp<-which(subMaster[,2]=="Stage IB")#9
subMaster[tmp,2]<-"Stage I"
tmp<-which(subMaster[,2]=="Stage IIA")#360
subMaster[tmp,2]<-"Stage II"
tmp<-which(subMaster[,2]=="Stage IIB")#251
subMaster[tmp,2]<-"Stage II"
tmp<-which(subMaster[,2]=="Stage IIIA")#153
subMaster[tmp,2]<-"Stage III"
tmp<-which(subMaster[,2]=="Stage IIIB")#29
subMaster[tmp,2]<-"Stage III"
tmp<-which(subMaster[,2]=="Stage IIIC")#63
subMaster[tmp,2]<-"Stage III"
unique(subMaster[,2])
rownames(subMaster)<-subMaster[,1]
subMaster<-subMaster[,-1]
###########################################
#absolute values
subMaster<-as.data.frame(subMaster)
Stage<-table(subMaster$Cluster.DBGS3.FLTR.RNSeq,subMaster$ajcc_pathologic_tumor_stage)
save(Stage, file="StageForBRCA-BSF-23July2015.RDATA")
Type<-table(subMaster$Cluster.DBGS3.FLTR.RNSeq,subMaster$TCGA.PAM50.RMethod.RNASeq)
save(Type,file="SubtypeForBRCA-BSF-23July2015.RDATA")
###########################
