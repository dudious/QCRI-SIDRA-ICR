#################################################################
###
### This Script clusters the Micro arrray date retrieved with TCGA assembler
### Input data :
### ../2 DATA/SUBSETS/..."
### Data is saved :
### ./3 ANALYSIS/CLUSTERING/MA
### Figures are saved :
### ./4 FIGURES/CLUSTERING
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR/")
  #Dependencies
  required.packages <- c("heatmap.plus")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  required.packages.BioC <- c("ConsensusClusterPlus")
  missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
  #source("http://bioconductor.org/biocLite.R")
  if(length(missing.packages)) biocLite(missing.packages)
  library(heatmap.plus)
  library(ConsensusClusterPlus)
  
# Load data
Geneset <- "16G.I.ID2" # SET GENESET HERE !!!!!!!!!!!!!!
cat ("Free memory at start : ", (memory.limit()-memory.size())/1000,"GB","\n")
load (paste0("./2 DATA/SUBSETS/MA.subset.",Geneset,".RData"))                                               # select subset to cluster
cat ("Free memory after loading data : ", (memory.limit()-memory.size())/1000,"GB","\n") 
# Hierarchical Clustering 

start.time <- Sys.time ()
  sHc <- hclust(ddist <- dist(MA.subset), method = "ward.D2")                              # hierachical clustering time :1.94 hours
end.time <- Sys.time ()
time <- end.time - start.time
print (time)

png(paste0("./4 FIGURES/CLUSTERING/BRCA.TCGA.MA.hclust.k4.",Geneset,".png"))
plot(	sHc,                                                                                        # PLot Clustering tree
		hang = -1,
		labels = FALSE,
		xlab = "euclidean/ward",
		sub = "",
		main = "raw data clustering");
(wwhere <- 4)                                                                                     # set K
hh <- mean(c(sHc$height[length(sHc$height)-wwhere+2], sHc$height[length(sHc$height)-wwhere+1]))   # calculate height of K
abline(h = hh, col = "red")                                                                       # add red line to plot
rawClusters <- as.factor(cutree(sHc, k = wwhere))                                                 # cut clustering at K
names(rawClusters) <- sHc$labels                                                                  # exreact genelist from sHc (hclust element)
levels(rawClusters) <- c("red","Blue","yellow","green")                                           # assing 4 colors (defined) instread of numbers to clusters 
#ALTERNATIVE : levels(rawClusters) <- (#colors()[(1:wwhere)*floor(length(colors())/wwhere)] )     # assing K colors (automatic) instread of numbers to clusters 
for(cc in unique(rawClusters)) rug(which(rawClusters[sHc$order] == cc), col = cc, lwd = 3)        # add color code to plot
dev.off()

cat ("Free memory after hclust : ", (memory.limit()-memory.size())/1000,"GB","\n") 
gc ()
cat ("Free memory after gc () : ", (memory.limit()-memory.size())/1000,"GB","\n") 

# Consensus Clustering
start.time <- Sys.time ()
setwd("~/Dropbox/BREAST_QATAR/3 ANALISYS/CLUSTERING/MA/")                                # temporary change in WD as path in CCP doesnt work
ans1 <- ConsensusClusterPlus(ddist,                                                               
                             maxK = 4,                                                            # set K
                             pItem = 0.8,
                             reps=1000,                                                            # set repeats
                             title=paste0("BRCA.TCGA.MA.k4.",Geneset,".reps1000"),                                 # Save location 
							               clusterAlg = "hc",                                                   # clustering Algorithm : Hierarchiocal clustering
							               innerLinkage = "ward.D2",
							               finalLinkage = "complete",
							               plot = 'pdf',                                                        # write resut to pdf (Alt:png)
							               writeTable = TRUE,
                             verbose = TRUE); 
setwd("~/Dropbox/BREAST_QATAR/")
end.time <- Sys.time ()
time <- end.time - start.time
print (time)
cat ("Free memory after Consensusclusterplus : ", (memory.limit()-memory.size())/1000,"GB","\n") 


