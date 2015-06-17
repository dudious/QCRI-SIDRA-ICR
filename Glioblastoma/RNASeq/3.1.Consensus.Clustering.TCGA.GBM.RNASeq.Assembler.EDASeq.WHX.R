#################################################################
###
### This Script clusters the RNASeq date retrieved with TCGA assembler
### and normalized using EDASeq
### Input data :
### ./2 DATA/SUBSETS/...
### Data is saved :
### ./3 ANALYSIS/CLUSTERING/RNASeq
### Figures are saved :
### ./4 FIGURES/CLUSTERING
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR/")
  #Dependencies
required.packages <- c("clue")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages("./1 CODE/R tools/clue_0.3-49.zip", repos=NULL)

  required.packages.BioC <- c("ConsensusClusterPlus")
  missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
  source("http://bioconductor.org/biocLite.R")
  if(length(missing.packages)) biocLite(missing.packages)

library(ConsensusClusterPlus)
library(clue)
source("./1 CODE/R tools/stefanofunctions.R")
  
# Load data
Geneset <- "ISGS" # SET GENESET HERE !!!!!!!!!!!!!!
# load (paste0("./2 DATA/SUBSETS/GBM/TCGA.GBM.RNASeq.subset.",Geneset,".RData"))                                               # select subset to cluster
load("./2 DATA/SUBSETS/GBM/TCGA.GBM.RNASeq.subset.ISGS.RData")
# Hierarchical Clustering 
start.time <- Sys.time ()
  sHc <- hclust(ddist <- dist(RNASeq.subset), method = "ward.D2")                              # hierachical clustering time :1.94 hours
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

setwd("~/Dropbox/BREAST_QATAR/3 ANALISYS/CLUSTERING/RNAseq/GBM")                                # temporary change in WD as path in CCP doesnt work
aCalinsky <- calinsky(sHc,gMax=10)
pdf(file = ("GBM.TCGA.aCalinsky.pdf"), width = 16, height = 6);
plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
dev.off()
#optimal k=2
# Consensus Clustering
start.time <- Sys.time ()

ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                               
                             maxK = 7,                                                            # set K
                             pItem = 0.8,
                             reps=5000,                                                            # set repeats
                             title="GBM.TCGA.EDASeq.k7.ISGS.reps5000",                                 # Save location 
							               clusterAlg = "hc",                                                   # clustering Algorithm : Hierarchiocal clustering
							               innerLinkage = "ward.D2",                                            # for color coding the clusters use tmyPal = ...
							               finalLinkage = "complete",
							               plot = 'pdf',                                                        # write resut to pdf (Alt:png)
							               writeTable = TRUE,
                             verbose = TRUE); 
setwd("~/Dropbox/BREAST_QATAR/")
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

save (ConsensusClusterObject,ddist,file="./3 ANALISYS/CLUSTERING/RNAseq/GBM/GBM.TCGA.EDASeq.k7.ISGS.reps5000/ConsensusClusterObject.Rdata")


