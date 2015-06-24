#################################################################
###
### This script creates a correlation matrix for the selcted genes
### Data is saved :
### 
### File to use :
### 
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

# load data
load ("./2 DATA/SUBSETS/RNASeq.subset.366G.RData")

# Corelation matrix
RNASeq.subset.cor <- cor (RNASeq.subset,method="pearson")

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
png("./4 FIGURES/correlation.366G.RNASeq.png",res=1200,height=12,width=12,unit="in")
corrplot(RNASeq.subset.cor,
         order="hclust",
         hclust.method="average",
         addrect=2,
         addgrid.col=NULL,
         method="square",
         tl.cex=0.23
         # addshade = NULL
        )
dev.off()




