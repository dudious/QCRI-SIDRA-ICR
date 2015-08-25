#Wilcoxcon Test BC subtype
#set directory
setwd("~/Dropbox/BREAST_QATAR")

# Load data
#1)DBGS3
#2)Master File
MasterFile<-read.csv("./3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF.RNASeq_subset_DBGS3.FLTR.Master.csv") 
MasterFile<-as.matrix(MasterFile)
subMasterFile<-MasterFile[,c(1,15,30,34)]
rownames(subMasterFile)<-subMasterFile[,1]
subMasterFile<-subMasterFile[,-1]
#3)No silent/silent/missense table
#load("......../BRCA.BSF.genes.mutation.nonsilent.table.RData")
#4)Log2 expression matrix RNASeq
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
RNASeq.subset<-RNASeq.NORM_Log2[DBGS3,]
#If i want results using z-score
#must have genes x samples
RNASeq.subsetZ <- (RNASeq.subset - rowMeans(RNASeq.subset))/apply(RNASeq.subset, 1, sd)
RNASeq.subset<-RNASeq.subsetZ
######################################################################
#Pool the  selected subtype (Basal-like, LuminalA-B, HER2-enriched, Normal-like)
who1<-which(subMasterFile[,2]=="Basal-like")
BL<-subMasterFile[who1,]
who2<-which(subMasterFile[,2]=="Luminal A")
LA<-subMasterFile[who2,]
who3<-which(subMasterFile[,2]=="HER2-enriched")
Her<-subMasterFile[who3,]
who4<-which(subMasterFile[,2]=="Luminal B")
LB<-subMasterFile[who4,]
who5<-which(subMasterFile[,2]=="Normal-like")
NL<-subMasterFile[who5,]
###########################################################################
#select the file of interest
#MutStatus<-BRCA.genes.mutations.nonsilent
#MutStatus<-BRCA.genes.mutations.silent
#MutStatus<-BRCA.genes.mutations.missense
MutStatus<-as.matrix(MutStatus)
Mutmatrix<- matrix(0,nrow=length(rownames(MutStatus)),ncol=length(colnames(MutStatus)))
rownames(Mutmatrix)<- rownames(MutStatus)
colnames(Mutmatrix)<- colnames(MutStatus)
Mutmatrix[MutStatus == 1] <- 1
#1=Mut; 0= Wt

#Compute ICR and CYT scores
ICRscore<-apply(RNASeq.subset,2,mean)
ICRscore<-as.matrix(ICRscore)
colnames(ICRscore)<-"DBGS3mean"
CYTScore<-apply(RNASeq.subset[c("PRF1","GZMA"),],2,mean)
CYTScore<-as.matrix(CYTScore)
colnames(CYTScore)<-"PRF1andGZMAmean"
#####################################################

#All samples
Mutmatrix<-t(Mutmatrix)
WilcoxonAll <- matrix(0, nrow= nrow(Mutmatrix), ncol=2)
rownames(WilcoxonAll) <- rownames(Mutmatrix)
colnames(WilcoxonAll) <- c("W","pValue")

for(i in 1:nrow(Mutmatrix)){
  Mut<-names(which(Mutmatrix[i,]==1))
  Wt<-names(which(Mutmatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(ICRscore[Mut,], ICRscore[Wt,])
    WilcoxonAll[i, 1] <- tmp$statistic
    WilcoxonAll[i, 2] <- tmp$p.value
  } else {
    WilcoxonAll[i, 1] <- "not calculated"
    WilcoxonAll[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonAll[,2],method="fdr")
WilcoxonAll <- cbind(WilcoxonAll, FDR=pvfdr)
#Order
Wilcoxon_O<-WilcoxonAll[order(WilcoxonAll[,3]),]
wwho<-which(Wilcoxon_O [,3] <= 0.05)#0
# Save
save(Wilcoxon_O,file="Wilcoxon-test_BRCA-ALL_DBGS3_missense.RDATA")
#NO CYT Analysis

#BASAL-LIKE
Mutmatrix<-t(Mutmatrix)
who<-which(rownames(Mutmatrix)%in% rownames(BL))#N=168
BLMutMatrix<-Mutmatrix[who,]
BLMutMatrix<-t(BLMutMatrix)
#output file ICR
WilcoxonBL <- matrix(0, nrow= nrow(BLMutMatrix), ncol=2)
rownames(WilcoxonBL) <- rownames(BLMutMatrix)
colnames(WilcoxonBL) <- c("W","pValue")

#test DBGS3
for(i in 1:nrow(BLMutMatrix)){
  Mut<-names(which(BLMutMatrix[i,]==1))
  Wt<-names(which(BLMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(ICRscore[Mut,], ICRscore[Wt,])
    WilcoxonBL[i, 1] <- tmp$statistic
    WilcoxonBL[i, 2] <- tmp$p.value
  } else {
    WilcoxonBL[i, 1] <- "not calculated"
    WilcoxonBL[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonBL[,2],method="fdr")
WilcoxonBL <- cbind(WilcoxonBL, FDR=pvfdr)
#Order
Wilcoxon_O<-WilcoxonBL[order(WilcoxonBL[,3]),]
wwho<-which(Wilcoxon_O [,3] <= 0.05)#0
# Save
save(Wilcoxon_O,file="Wilcoxon-test_BRCA-BASAL-LIKE_DBGS3_missense.RDATA")

#output file for CYT
WilcoxonBL_CYT <- matrix(0, nrow= nrow(BLMutMatrix), ncol=2)
rownames(WilcoxonBL_CYT) <- rownames(BLMutMatrix)
colnames(WilcoxonBL_CYT) <- c("W","pValue")

#test CYT
for(i in 1:nrow(BLMutMatrix)){
  Mut<-names(which(BLMutMatrix[i,]==1))
  Wt<-names(which(BLMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(CYTScore[Mut,], CYTScore[Wt,])
    WilcoxonBL_CYT[i, 1] <- tmp$statistic
    WilcoxonBL_CYT[i, 2] <- tmp$p.value
  } else {
    WilcoxonBL_CYT[i, 1] <- "not calculated"
    WilcoxonBL_CYT[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonBL_CYT[,2],method="fdr")
WilcoxonBL_CYT <- cbind(WilcoxonBL_CYT, FDR=pvfdr)
#Order
Wilcoxon_O_CYT<-WilcoxonBL_CYT[order(WilcoxonBL_CYT[,3]),]
wwho<-which(Wilcoxon_O_CYT [,3] <= 0.05)#0
# Save
save(Wilcoxon_O_CYT,file="Wilcoxon-test_BRCA-BASAL-LIKE_CYT-SILENT.RDATA")

################################################
#Her2
who<-which(rownames(Mutmatrix)%in% rownames(Her))#N=79
HerMutMatrix<-Mutmatrix[who,]
HerMutMatrix<-t(HerMutMatrix)

#output file
WilcoxonHer <- matrix(0, nrow= nrow(HerMutMatrix), ncol=2)
rownames(WilcoxonHer) <- rownames(HerMutMatrix)
colnames(WilcoxonHer) <- c("W","pValue")

#test DBGS3
for(i in 1:nrow(HerMutMatrix)){
  Mut<-names(which(HerMutMatrix[i,]==1))
  Wt<-names(which(HerMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(ICRscore[Mut,], ICRscore[Wt,])
    WilcoxonHer[i, 1] <- tmp$statistic
    WilcoxonHer[i, 2] <- tmp$p.value
  } else {
    WilcoxonHer[i, 1] <- "not calculated"
    WilcoxonHer[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonHer[,2],method="fdr")
WilcoxonHer <- cbind(WilcoxonHer, FDR=pvfdr)
#Order
Wilcoxon_O<-WilcoxonHer[order(WilcoxonHer[,3]),]
wwho<-which(Wilcoxon_O [,3] <= 0.05)#0
# Save
save(Wilcoxon_O,file="Wilcoxon-test_BRCA-Her2Enriched_DBGS3-missense.RDATA")


#output file for CYT
WilcoxonHer_CYT <- matrix(0, nrow= nrow(HerMutMatrix), ncol=2)
rownames(WilcoxonHer_CYT) <- rownames(HerMutMatrix)
colnames(WilcoxonHer_CYT) <- c("W","pValue")

#test CYT
for(i in 1:nrow(HerMutMatrix)){
  Mut<-names(which(HerMutMatrix[i,]==1))
  Wt<-names(which(HerMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(CYTScore[Mut,], CYTScore[Wt,])
    WilcoxonHer_CYT[i, 1] <- tmp$statistic
    WilcoxonHer_CYT[i, 2] <- tmp$p.value
  } else {
    WilcoxonHer_CYT[i, 1] <- "not calculated"
    WilcoxonHer_CYT[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonHer_CYT[,2],method="fdr")
WilcoxonHer_CYT <- cbind(WilcoxonHer_CYT, FDR=pvfdr)
#Order
Wilcoxon_O_CYT<-WilcoxonHer_CYT[order(WilcoxonHer_CYT[,3]),]
wwho<-which(Wilcoxon_O_CYT [,3] <= 0.05)#0
# Save
save(Wilcoxon_O_CYT,file="Wilcoxon-test_BRCA-Her2Enriched_CYT-SILENT.RDATA")
###########################################
#Luminal A
who<-which(rownames(Mutmatrix)%in% rownames(LA))#N=503
LAMutMatrix<-Mutmatrix[who,]
LAMutMatrix<-t(LAMutMatrix)

#output file
WilcoxonLA <- matrix(0, nrow= nrow(LAMutMatrix), ncol=2)
rownames(WilcoxonLA) <- rownames(LAMutMatrix)
colnames(WilcoxonLA) <- c("W","pValue")

#test DBGS3
for(i in 1:nrow(LAMutMatrix)){
  Mut<-names(which(LAMutMatrix[i,]==1))
  Wt<-names(which(LAMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(ICRscore[Mut,], ICRscore[Wt,])
    WilcoxonLA[i, 1] <- tmp$statistic
    WilcoxonLA[i, 2] <- tmp$p.value
  } else {
    WilcoxonLA[i, 1] <- "not calculated"
    WilcoxonLA[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonLA[,2],method="fdr")
WilcoxonLA <- cbind(WilcoxonLA, FDR=pvfdr)
#Order
Wilcoxon_O<-WilcoxonLA[order(WilcoxonLA[,3]),]
wwho<-which(Wilcoxon_O [,3] <= 0.05)#0
# Save
save(Wilcoxon_O,file="Wilcoxon-test_BRCA-LuminalA_DBGS3-missense.RDATA")

#output file for CYT
WilcoxonLA_CYT <- matrix(0, nrow= nrow(LAMutMatrix), ncol=2)
rownames(WilcoxonLA_CYT) <- rownames(LAMutMatrix)
colnames(WilcoxonLA_CYT) <- c("W","pValue")

#test CYT
for(i in 1:nrow(LAMutMatrix)){
  Mut<-names(which(LAMutMatrix[i,]==1))
  Wt<-names(which(LAMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(CYTScore[Mut,], CYTScore[Wt,])
    WilcoxonLA_CYT[i, 1] <- tmp$statistic
    WilcoxonLA_CYT[i, 2] <- tmp$p.value
  } else {
    WilcoxonLA_CYT[i, 1] <- "not calculated"
    WilcoxonLA_CYT[i, 2] <- "not calculated"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonLA_CYT[,2],method="fdr")
WilcoxonLA_CYT <- cbind(WilcoxonLA_CYT, FDR=pvfdr)
#Order
Wilcoxon_O_CYT<-WilcoxonLA_CYT[order(WilcoxonLA_CYT[,3]),]
wwho<-which(Wilcoxon_O_CYT [,3] <= 0.05)#0
# Save
save(Wilcoxon_O_CYT,file="Wilcoxon-test_BRCA-LuminalA_CYT-SILENT.RDATA")
########################################
#Luminal B

who<-which(rownames(Mutmatrix)%in% rownames(LB))#N=168
LBMutMatrix<-Mutmatrix[who,]
LBMutMatrix<-t(LBMutMatrix)

#output file
WilcoxonLB <- matrix(0, nrow= nrow(LBMutMatrix), ncol=2)
rownames(WilcoxonLB) <- rownames(LBMutMatrix)
colnames(WilcoxonLB) <- c("W","pValue")

#test DBGS3
for(i in 1:nrow(LBMutMatrix)){
  Mut<-names(which(LBMutMatrix[i,]==1))
  Wt<-names(which(LBMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(ICRscore[Mut,], ICRscore[Wt,])
    WilcoxonLB[i, 1] <- tmp$statistic
    WilcoxonLB[i, 2] <- tmp$p.value
  } else {
    WilcoxonLB[i, 1] <- "not calcuLBted"
    WilcoxonLB[i, 2] <- "not calcuLBted"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonLB[,2],method="fdr")
WilcoxonLB <- cbind(WilcoxonLB, FDR=pvfdr)
#Order
Wilcoxon_O<-WilcoxonLB[order(WilcoxonLB[,3]),]
wwho<-which(Wilcoxon_O [,3] <= 0.05)#0
# Save
save(Wilcoxon_O,file="Wilcoxon-test_BRCA-LuminaLB_DBGS3-missense.RDATA")


#output file for CYT
WilcoxonLB_CYT <- matrix(0, nrow= nrow(LBMutMatrix), ncol=2)
rownames(WilcoxonLB_CYT) <- rownames(LBMutMatrix)
colnames(WilcoxonLB_CYT) <- c("W","pValue")

#test CYT
for(i in 1:nrow(LBMutMatrix)){
  Mut<-names(which(LBMutMatrix[i,]==1))
  Wt<-names(which(LBMutMatrix[i,]==0))
  if(length(Mut)>0 & length(Wt) > 0){
    tmp <- wilcox.test(CYTScore[Mut,], CYTScore[Wt,])
    WilcoxonLB_CYT[i, 1] <- tmp$statistic
    WilcoxonLB_CYT[i, 2] <- tmp$p.value
  } else {
    WilcoxonLB_CYT[i, 1] <- "not calcuLBted"
    WilcoxonLB_CYT[i, 2] <- "not calcuLBted"
  }
  print(i)
}
pvfdr<-p.adjust(WilcoxonLB_CYT[,2],method="fdr")
WilcoxonLB_CYT <- cbind(WilcoxonLB_CYT, FDR=pvfdr)
#Order
Wilcoxon_O_CYT<-WilcoxonLB_CYT[order(WilcoxonLB_CYT[,3]),]
wwho<-which(Wilcoxon_O_CYT [,3] <= 0.05)#0
# Save
save(Wilcoxon_O_CYT,file="Wilcoxon-test_BRCA-LuminaLB_CYT-SILENT.RDATA")
###########################################
#Normal-like
who<-which(rownames(Mutmatrix)%in% rownames(NL))#N=0
##############################################################

#plotP53
tp53<-which(colnames(Mutmatrix)=="TP53")
P53<-as.matrix(Mutmatrix[,tp53])
colnames(P53)<-"MutStatusTP53ALL"
mut<-which(P53[,1]==1)
P53[mut,]<-"MUT"
wt<-which(P53[,1]==0)
P53[wt,]<-"WT"
P53<-merge(P53,ICRscore,by="row.names",all.x=TRUE)
rownames(P53)<-P53[,1]
P53$Row.names<-NULL
P53<-merge(P53,CYTScore,by="row.names",all.x=TRUE)
rownames(P53)<-P53[,1]
P53$Row.names<-NULL
################################
tp53BL<-which(rownames(BLMutMatrix)=="TP53")
Basal<-as.matrix(BLMutMatrix[tp53BL,])
colnames(Basal)<-"MutStatusTP53BL"
mut<-which(Basal[,1]==1)
Basal[mut,]<-"MUT"
wt<-which(Basal[,1]==0)
Basal[wt,]<-"WT"
Basal<-merge(Basal,ICRscore,by="row.names",all.x=TRUE)
rownames(Basal)<-Basal[,1]
Basal$Row.names<-NULL
Basal<-merge(Basal,CYTScore,by="row.names",all.x=TRUE)
rownames(Basal)<-Basal[,1]
Basal$Row.names<-NULL
####
tp53LA<-which(rownames(LAMutMatrix)=="TP53")
LumA<-as.matrix(LAMutMatrix[tp53LA,])
colnames(LumA)<-"MutStatusTP53LUMA"
mut<-which(LumA[,1]==1)
LumA[mut,]<-"MUT"
wt<-which(LumA[,1]==0)
LumA[wt,]<-"WT"
LumA<-merge(LumA,ICRscore,by="row.names",all.x=TRUE)
rownames(LumA)<-LumA[,1]
LumA$Row.names<-NULL
LumA<-merge(LumA,CYTScore,by="row.names",all.x=TRUE)
rownames(LumA)<-LumA[,1]
LumA$Row.names<-NULL
####
tp53LB<-which(rownames(LBMutMatrix)=="TP53")
LumB<-as.matrix(LBMutMatrix[tp53LB,])
colnames(LumB)<-"MutStatusTP53LUMB"
mut<-which(LumB[,1]==1)
LumB[mut,]<-"MUT"
wt<-which(LumB[,1]==0)
LumB[wt,]<-"WT"
LumB<-merge(LumB,ICRscore,by="row.names",all.x=TRUE)
rownames(LumB)<-LumB[,1]
LumB$Row.names<-NULL
LumB<-merge(LumB,CYTScore,by="row.names",all.x=TRUE)
rownames(LumB)<-LumB[,1]
LumB$Row.names<-NULL

####
tp53H2<-which(rownames(HerMutMatrix)=="TP53")
HER2E<-as.matrix(HerMutMatrix[tp53H2,])
colnames(HER2E)<-"MutStatusTP53HER2"
mut<-which(HER2E[,1]==1)
HER2E[mut,]<-"MUT"
wt<-which(HER2E[,1]==0)
HER2E[wt,]<-"WT"
HER2E<-merge(HER2E,ICRscore,by="row.names",all.x=TRUE)
rownames(HER2E)<-HER2E[,1]
HER2E$Row.names<-NULL
HER2E<-merge(HER2E,CYTScore,by="row.names",all.x=TRUE)
rownames(HER2E)<-HER2E[,1]
HER2E$Row.names<-NULL
###################################################
#Plotting:make the plot; change the data
pdf("ICRvsTP53Mutstatus-missenseinAllsubtype.pdf")
boxplot(DBGS3mean~ MutStatusTP53ALL, data= P53, col= c("gray","White"),
        main="Association between TP53 mutational status and ICR score",
        ylab="ICR score",
        xlab="All BC subtypes")
text(x=1,y=2.2,labels="q= 0.003", col="black",adj=1, pos=3,cex=1.2, font=4)
dev.off()
#######################
#subtype.order = c("Basal-like", "HER2-enriched","Luminal A", "Luminal B")
#colors = c( "darkorchid","darkgoldenrod3","darkblue","dodgerblue2")
#Change the data and the information
pdf("ICRvsTP53Mutstatus-missenseinHER2Enriched.pdf")
boxplot(DBGS3mean~ MutStatusTP53HER2, data= HER2E, col= c("darkgoldenrod3","White"),
        main="Association between TP53 mutational status and ICR score",
        ylab="ICR score",
        xlab="HER2-enriched")
text(x=1.5,labels="q= 0.964", col="black",adj=1, pos=3,cex=1.2, font=4)
dev.off()
