#################################################################
###
### This script creates a correlation matrix for the selcted genes
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
required.packages <- c("corrplot")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (corrplot)

# Set Parameters
Cancerset <- "SKCM"           # Select the Cancer to use
Geneset   <- "DBGS3"          # Select the genset to use
Filter    <- "TRUE"           # Use Pre-Clustering Filter "TRUE" OR "FALSE"  (setup filter in "2.3.Exclude.Clinical" script)
BRCA.Filter <- "PCF"          # "PCF" or "BSF" Pancer or Breast specific

# load data
#load ("./2 DATA/SUBSETS/METABRIC/METABRIC.RNASEQ.DATA.1.genesubset.25G.DB.RData") # Select subset here !!!!!
#load ("./2 DATA/SUBSETS/METABRIC/METABRIC.RNASEQ.DATA.2.genesubset.25G.DB.RData") # Select subset here !!!!!
#RNASEQ.DATA.ALL.subset <- rbind(RNASEQ.DATA.1.subset,RNASEQ.DATA.2.subset)
load (paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".RData")) # Select subset here !!!!!

# Corelation matrix
RNASeq.subset.cor <- cor (RNASeq.subset,method="spearman")

# cor significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
RNASeq.subset.cor.sign <- cor.mtest(RNASeq.subset.cor, 0.95)

# Correlation plot
png(paste0("./4 FIGURES/CORRELATION/",Geneset,"/correlation.",Cancerset,".",Geneset,".png"),res=600,height=6,width=6,unit="in")  #adjust output file names here !!!!!
cex.before <- par("cex")
par(cex = 0.45)
col1 = colorRampPalette(c("green", "white", "orange"))
lims=c(-1,1)
if (length(RNASeq.subset.cor[RNASeq.subset.cor<0]) == 0) {lims=c(0,1)} 
corrplot.mixed (RNASeq.subset.cor,
               #type="lower",
               #p.mat = RNASeq.subset.cor.sign[[1]],    # add significance to correlations
               col = col1(100),
               lower = "pie",
               upper ="circle",
               order="FPC",
               cl.lim=lims,                      # only positive correlations
               tl.pos ="lt",
               tl.col = "#c00000",
               #insig="pch",                          # remove insignificant correlations
               tl.cex = 1/par("cex"),
               cl.cex = 1/par("cex"),
               title = paste0("Spearman TCGA-RNASeq (",Cancerset,".",Geneset," avg:",round(mean(RNASeq.subset.cor),2),")"),
               cex.main = 1.4/par("cex"),
               mar=c(5.1,4.1,4.1,2.1))
par(cex = cex.before)
dev.off()




