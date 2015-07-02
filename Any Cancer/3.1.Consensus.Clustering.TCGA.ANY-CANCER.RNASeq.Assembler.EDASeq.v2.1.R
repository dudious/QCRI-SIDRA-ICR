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
setwd("~/Dropbox/BREAST_QATAR/") #using symbolic link on R-server, no actual dropbox (ln -s /cancer_data/Cancer\ Immunesignature\ QCRI-SIDRA/ /home/whendrickx/Dropbox/BREAST_QATAR)
#Dependencies
  required.packages <- c("clue")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages( "~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/clue_0.3-49.zip", repos=NULL,type= "win.binary") #windows
  if (Sys.info()['sysname']=="Linux"){if(length(missing.packages)) install.packages(missing.packages)} #Linux
  required.packages.BioC <- c("ConsensusClusterPlus")
  missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
  source("http://bioconductor.org/biocLite.R")
  if(length(missing.packages)) biocLite(missing.packages)

library(ConsensusClusterPlus)
library(clue)
source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/stefanofunctions.R")

# Set Parameters
Cancerset <- "BRCA"           # Select the Cancer to use
Geneset   <- "DBGS1"          # Select the genset to use
Filter    <- "TRUE"           # Use Pre-Clustering Filter "TRUE" OR "FALSE"  (setup filter in "2.3.Exclude.Clinical" script)

# Load data
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".RData"))                  # load subset data to cluster
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv")) 

# Filter for sample exlusion parameter
Selected.samples <- Clinical.data[Clinical.data$exclude.pre=="No",1]
if (Filter == "TRUE"){
  Geneset <- paste0(Geneset,".FLTR")
  RNASeq.subset <- RNASeq.subset[rownames(RNASeq.subset) %in% Selected.samples,]
}

# Hierarchical Clustering 
sHc <- hclust(ddist <- dist(RNASeq.subset), method = "ward.D2")                                                     # hierachical clustering

# Change WD for Consensusclustering
dir.create(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/"), showWarnings = FALSE)
setwd(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/"))                                                      # temporary change in WD as path in CCP doesnt work

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
                             title=paste0("",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000"),   # Output filename (no path)
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


