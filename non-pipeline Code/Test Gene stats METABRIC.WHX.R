#################################################################
###
### This script creates a ..... for the selcted genes
### Data used:
### ./2 DATA/SUBSETS/...
### Output :
### ./4 FIGURES/...
###
#################################################################


# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
## dependencies
#required.packages <- c("ggplot")
#missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
#if(length(missing.packages)) install.packages(missing.packages)

required.packages.BioC <- c("ggplot2")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

library (ggplot2)

# load data
load ("./2 DATA/SUBSETS/METABRIC.RNASEQ.DATA.2.genesubset.16G.I.RData") # Select subset here !!!!!
stats <- rbind(sapply(RNASEQ.DATA.2.subset,mean),sapply(RNASEQ.DATA.2.subset,sd),sapply(RNASEQ.DATA.2.subset,min),sapply(RNASEQ.DATA.2.subset,max))
rownames (stats) <- c("Mean","sd","min","max")             
      

ggplot(RNASEQ.DATA.2.subset, aes(x=cond, y=rating)) + geom_boxplot()
boxplot (RNASEQ.DATA.2.subset)
