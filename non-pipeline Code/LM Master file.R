# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

load("./2 DATA/LM.BRCA/LM.Dataset.split.Rdata")


Expression.Data$Gene <- Gene.Meta.data$Symbol[match(rownames(Expression.Data),Gene.Meta.data$Affy_Probe_ID)]
Expression.Data <- Expression.Data[-which(Expression.Data$Gene == 0),]
class (Expression.Data$Gene) <- "factor"
Expression.Data.colapsed.toGenes <- aggregate( . ~ Gene,
                                              data = Expression.Data,
                                              FUN = mean)
Expression.Data.colapsed.toGenes <- t(Expression.Data.colapsed.toGenes)
colnames(Expression.Data.colapsed.toGenes) <- Expression.Data.colapsed.toGenes[1,]
Expression.Data.colapsed.toGenes <- Expression.Data.colapsed.toGenes[-1,]
Expression.Data.colapsed.toGenes <- as.matrix(Expression.Data.colapsed.toGenes)
mode(Expression.Data.colapsed.toGenes) <- "numeric"
cluster.assignment <- read.csv("./3 ANALISYS/CLUSTERING/MA/LM.Dataset/LM.Dataset.MA.k7.DBGS3.reps5000/LM.Dataset.MA.k7.DBGS3.reps5000.k=4.consensusClass.ICR.csv")
Sample.Meta.data$Cluster <- cluster.assignment$Group[match(rownames(Sample.Meta.data),cluster.assignment$PatientID)]
LM.MASTER <- merge(Sample.Meta.data,Expression.Data.colapsed.toGenes,by="row.names")

write.csv (LM.MASTER,file = "./3 ANALISYS/MASTER FILES/LM_DATA_bygene.csv")
