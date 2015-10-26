# 6/10/2015
# DEGs in ICR4 vs ICR1
# DEGs in ICR4 vs ICR123
#EDGEG Method


library(edgeR)

#load the functions

TimeUse <-function(func){ 
  time1 <- proc.time() 
  result <-func
  time2 <- proc.time()
  show(timeelapsed <- time2-time1)
}

edgeRDEA4GSEAOff <- function(geneCounts, groups, baseLine){
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  DEA <- DGEList(geneCounts, group = groups)
  DEA <- calcNormFactors(DEA)
  keep <- rowSums(cpm(DEA)>1)>=10
  DEA <- DEA[keep,]
  DEA$samples$lib.size <- colSums(DEA$counts)
  message("Estimating GLM Common Dispersion")
  TimeUse(DEA <- estimateGLMCommonDisp(DEA, design = design))
  message("Estimating GLM Tagwise Dispersion")
  TimeUse(DEA <- estimateGLMTagwiseDisp(DEA, design = design))
  message("Fitting the GLM Model")
  TimeUse(fit <- glmFit(DEA, design = design))
  contrast <- rep(1, 2)
  blId <- which(colnames(design) == baseLine)
  contrast[blId] <- -1
  message("Maximum Likelihood Estimate")
  TimeUse(test <- glmLRT(fit, contrast = contrast))
  ans <- topTags(test, n = nrow(DEA$counts))
  ans <- ans$table
  diffGenes <- abs(ans[, "logFC"]) >= 1 & ans[, "FDR"] <= 0.01
  detags <- rownames(ans)[diffGenes]
  pdf("smearPlot.pdf")
  plotSmear(test,de.tags=detags)
  abline(h=c(-1,1),col="blue")
  dev.off()
  return(ans)
}

##############################################################################

setwd("~/Dropbox/BREAST_QATAR")
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")
normCounts<-RNASeq.NORM
selectedGenes <- which(rowMeans(normCounts) >= quantile(rowMeans(normCounts), 0.25))
normCounts <- normCounts[selectedGenes, ]
normCounts <- normCounts[, order(colnames(normCounts))]
#View(head(normCounts))# genes x sample
dim(normCounts)#15241  1097
################################################################
FinalClassification <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/BRCA.BSF2/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000/BRCA.BSF2.TCGA.EDASeq.k7.DBGS3.FLTR.reps5000.k=4.consensusClass.ICR.csv",stringsAsFactors=F)
rownames(FinalClassification)<-FinalClassification[,1]
FinalClassification$X<-NULL
colnames(FinalClassification)<-c("Sample","Cluster")
FinalClassification<-FinalClassification[order(FinalClassification[,2]),]
#only for 4vs1
#tmp1<-which(FinalClassification$Cluster=="ICR1" )
#tmp2<-which(FinalClassification$Cluster=="ICR4" )
#FinalClassification<-FinalClassification[c(tmp1,tmp2),]

groups <- rep("",nrow(FinalClassification))
names(groups)<- rownames(FinalClassification)
#View(groups)

#ICR4 vs ICR123
G1<- FinalClassification$Sample[FinalClassification[,2]=="ICR4"]
groups[G1]<-"G1"
G2<- FinalClassification$Sample[FinalClassification[,2]!="ICR4"]
groups[G2]<-"G2"
groups<-as.factor(groups)

system.time(DEGinICR4vs123 <- edgeRDEA4GSEAOff(normCounts[,names(groups)], groups = groups,  
                                        baseLine = "G2"))

DEGinICR4vs123<-DEGinICR4vs123[order(DEGinICR4vs123$logFC,decreasing = T),]
write.csv(DEGinICR4vs123,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/DEGinICR4vs123BRCA_BSF2_Ines.csv")
save(DEGinICR4vs123,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/DEGinICR4vs123.RDATA")

#FILTER
ICR4vs123_UPs <- DEGinICR4vs123[DEGinICR4vs123$FDR < 0.05 & DEGinICR4vs123$logFC >= 1.5, ]
ICR4vs123_UPs <- ICR4vs123_UPs[order(-ICR4vs123_UPs$logFC), ]
write.csv(ICR4vs123_UPs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs123_UPs_BRCA_BSF2_Ines.csv")
save(ICR4vs123_UPs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs123_UPs.RDATA")


ICR4vs123_DOWNs <- DEGinICR4vs123[DEGinICR4vs123$FDR < 0.05 & DEGinICR4vs123$logFC <= -1.5, ]
ICR4vs123_DOWNs <- ICR4vs123_DOWNs[order(-ICR4vs123_DOWNs$logFC), ]
write.csv(ICR4vs123_DOWNs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs123_DOWNs_BRCA_BSF2_Ines.csv")
save(ICR4vs123_DOWNs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs123_DOWNs.RDATA")

#############################################################################################################

#ICR4 vs ICR1
G1<- FinalClassification$Sample[FinalClassification[,2]=="ICR4"]
groups[G1]<-"G1"
G2<- FinalClassification$Sample[FinalClassification[,2]=="ICR1"]
groups[G2]<-"G2"
groups<-as.factor(groups)

system.time(DEGinICR4vs1 <- edgeRDEA4GSEAOff(normCounts[,names(groups)], groups = groups,  
                                               baseLine = "G2"))

DEGinICR4vs1<-DEGinICR4vs1[order(DEGinICR4vs1$logFC,decreasing = T),]
write.csv(DEGinICR4vs1,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/DEGinICR4vs1BRCA_BSF2_Ines.csv")
save(DEGinICR4vs1,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/DEGinICR4vs1.RDATA")

#FILTER
ICR4vs1_UPs <- DEGinICR4vs1[DEGinICR4vs1$FDR < 0.05 & DEGinICR4vs1$logFC >= 1.5, ]
ICR4vs1_UPs <- ICR4vs1_UPs[order(-ICR4vs1_UPs$logFC), ]
write.csv(ICR4vs1_UPs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs1_UPs_BRCA_BSF2_Ines.csv")
save(ICR4vs1_UPs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs1_UPs.RDATA")


ICR4vs1_DOWNs <- DEGinICR4vs1[DEGinICR4vs1$FDR < 0.05 & DEGinICR4vs1$logFC <= -1.5, ]
ICR4vs1_DOWNs <- ICR4vs1_DOWNs[order(-ICR4vs1_DOWNs$logFC), ]
write.csv(ICR4vs1_DOWNs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs1_DOWNs_BRCA_BSF2_Ines.csv")
save(ICR4vs1_DOWNs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs1_DOWNs.RDATA")

