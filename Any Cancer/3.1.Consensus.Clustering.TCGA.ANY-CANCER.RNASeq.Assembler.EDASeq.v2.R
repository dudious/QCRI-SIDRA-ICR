#################################################################
###
### This Script clusters the ",Cancerset," RNASeq date retrieved with 
### TCGA assembler and normalized using EDASeq
### Input data :
### ./2 DATA/SUBSETS/...
### Data is saved :
### ./3 ANALYSIS/CLUSTERING/RNASeq/",Cancerset,"
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

# Set Parameters
Cancerset <- "COAD-hiseq"
Geneset <- "ISGS1"                                                                           # Select the genset to use
  
# Load data
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".RData"))             # load subset data to cluster
dir.create(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/"), showWarnings = FALSE)

# Hierarchical Clustering 
sHc <- hclust(ddist <- dist(RNASeq.subset), method = "ward.D2")                              # hierachical clustering

# Change WD for Consensusclustering
setwd(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/"))                            # temporary change in WD as path in CCP doesnt work

# Perform Calinsky for optimal K
aCalinsky <- calinsky(sHc,gMax=10)
pdf(file = paste0("",Cancerset,".TCGA.",Geneset,".Calinsky.pdf"), width = 16, height = 6);
  plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
  text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
dev.off()

# Consensus Clustering
start.time <- Sys.time ()
ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                               
                             maxK = 7,                                                            # set K
                             pItem = 0.8,
                             reps=5000,                                                           # set repeats
                             title=paste0("",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000"),            # Output filename (no path)
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

save (ConsensusClusterObject,ddist,file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/ConsensusClusterObject.Rdata"))


