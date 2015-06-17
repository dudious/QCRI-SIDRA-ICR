#################################################################
###
### This Script creates a subset of the LM data 
### for a selected set  genes.(Gene_selection_XXX.txt)
### source data :
### "./2 DATA/SUBSETS/Gene_selection_v2.2.txt"
### ~/Dropbox/Data_LM/LM.Dataset.split.Rdata"
### 
### Results are saved in
### ./2 DATA/SUBSETS/
### File to use :
### 
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

required.packages <- c("clue")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages( "./1 CODE/R tools/clue_0.3-49.zip", repos=NULL,type= "win.binary")

required.packages.BioC <- c("ConsensusClusterPlus")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

library(ConsensusClusterPlus)
library(clue)
source("./1 CODE/R tools/stefanofunctions.R")

# Parameters
Cancerset <- "LM.Dataset"
Geneset <- "DBGS1"

# load data
load ("~/Dropbox/Data_LM/LM.Dataset.split.Rdata")                                
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_v2.2.txt")                                 # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,Geneset]==1),1])

# check availabilety of the genes in the dataset
available.genes.MA <- gene.list.selected[which(gene.list.selected %in% Gene.Meta.data$Symbol)]
unavailable.genes.MA <- gene.list.selected[-which(gene.list.selected %in% Gene.Meta.data$Symbol)]
c("B7-H","B7-H1","B7H1","PD-L1","PDL1") %in% Gene.Meta.data$Symbol                                # check for CD274 alternative names

## Subset data
MA.probes.subset <- Gene.Meta.data[Gene.Meta.data$Symbol %in% available.genes.MA,]
MA.subset <- t(Expression.Data[rownames(Expression.Data) %in% MA.probes.subset$Affy_Probe_ID,])
mode(MA.subset) <- "numeric"

## report
print (paste0("Geneset Selected : ",Geneset))
print ("Genes selected for MA : ")
print (available.genes.MA)
print ("Genes missing for MA :")
print (unavailable.genes.MA)

# save subsetted data
dir.create(paste0("./2 DATA/SUBSETS/",Cancerset,"/"), showWarnings = FALSE)
save (MA.subset,file=paste0("./2 DATA/SUBSETS/",Cancerset,"/",Cancerset,".MA.subset.",Geneset,".RData"))                #adjust output file names here !!!!!

# Hierarchical Clustering 
MA.subset[which(is.na(MA.subset))] <- 0 # set NA to 0
sHc <- hclust(ddist <- dist(MA.subset), method = "ward.D2")                              # hierachical clustering

# Change WD for Consensusclustering
dir.create(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/"), showWarnings = FALSE)
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
                                               title=paste0("",Cancerset,".MA.k7.",Geneset,".reps5000"),            # Output filename (no path)
                                               clusterAlg = "hc",                                                   # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                            # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                        # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE); 
setwd("~/Dropbox/BREAST_QATAR/")
end.time <- Sys.time ()
time <- end.time - start.time
print (time)  #Time difference of 2.179532 hours

save (ConsensusClusterObject,ddist,file=paste0("./3 ANALISYS/CLUSTERING/MA/",Cancerset,"/",Cancerset,".MA.k7.",Geneset,".reps5000/ConsensusClusterObject.Rdata"))

