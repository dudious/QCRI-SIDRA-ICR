#################################################################
###
### This script analyses the duplicates in the PANCANCER RNASeq Matrix 
###
### 
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")
library(ConsensusClusterPlus)
library(clue)
install.packages("dendextend")
library(dendextend)

# Set Parameters
Cancerset         = "GBM"
Geneset           = "DBGS3"   
DL.Method         = "PANCANCER"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Split"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer

## Load Data
load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)

cor.test(RNASeq.subset["TCGA-06-0156-01A-02R-1849-01",],RNASeq.subset["TCGA-06-0156-01A-03R-1849-01",])
cor.test(RNASeq.subset["TCGA-06-0211-01A-01R-1849-01",],RNASeq.subset["TCGA-06-0211-01B-01R-1849-01",])
cor.test(RNASeq.subset["TCGA-06-0156-01A-02R-1849-01",],RNASeq.subset["TCGA-06-0211-01B-01R-1849-01",])

# Clustering 
setwd(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/"))
sHc <- hclust(ddist <- dist(RNASeq.subset), method = "ward.D2")  
ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                               maxK = 7,                                                                          # set K
                                               pItem = 0.8,
                                               reps=5000,                                                                         # set repeats
                                               title=paste0("",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000"),   # Output filename (no path)
                                               clusterAlg = "hc",                                                   # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                            # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                        # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE);
dev.new()
dend2 <- color_labels(as.dendrogram(sHc), labels = c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01") , col = 2) 
plot(as.dendrogram(sHc))
plot(dend2)

dend2 <- color_labels(as.dendrogram(sHc), labels = c("TCGA-06-0211-01A-01R-1849-01","TCGA-06-0211-01B-01R-1849-01") , col = 2) 
plot(dend2)

Cancerset         = "OV"
load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)
sHc <- hclust(ddist <- dist(RNASeq.subset), method = "ward.D2")
dev.new()
dend2 <- color_labels(as.dendrogram(sHc), labels = c("TCGA-23-1023-01A-02R-1564-13" , "TCGA-23-1023-01R-01R-1564-13" ) , col = 2) 
plot(dend2)

Cancerset         = "LUSC"
load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)
sHc <- hclust(ddist <- dist(RNASeq.subset), method = "ward.D2")
dev.new()
dend2 <- color_labels(as.dendrogram(sHc), labels = c("TCGA-21-1076-01A-01R-0692-07" , "TCGA-21-1076-01A-02R-0692-07" ) , col = 2) 
plot(dend2)
