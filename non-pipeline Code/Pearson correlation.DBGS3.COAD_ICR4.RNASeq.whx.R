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

# load data
RNASeq.data <- read.csv ("./3 ANALISYS/TCGA.COAD-ICR4_DBGS3_RNASeq.Data.csv") # Select subset here !!!!!
#RNASeq.subset <- as.matrix(RNASeq.data[,1:(ncol(RNASeq.data)-1)])
#RNASeq.subset <- as.matrix(RNASeq.data[RNASeq.data$X=="LOW",1:(ncol(RNASeq.data)-1)])
inhib.genes <- c("PDCD1","CTLA4","CD274","FOXP3","IDO1")
which(colnames(RNASeq.data) %in% inhib.genes)
RNASeq.subset <- as.matrix(RNASeq.data[,which(colnames(RNASeq.data) %in% inhib.genes)])
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
png("./4 FIGURES/CORRELATION/correlation.COAD_ICR4.all.IMSUG.RNASeq.png",res=600,height=6,width=6,unit="in")  #adjust output file names here !!!!!
cex.before <- par("cex")
par(cex = 0.45)
corrplot.mixed(RNASeq.subset.cor,
               #p.mat = RNASeq.subset.cor.sign[[1]],    # add significance to correlations
               lower = "circle",
               upper ="number",
               order="FPC",
               #cl.lim=c(0,1),                          # only positive correlations
               tl.pos ="lt",
               #insig="blank",                          # remove insignificant correlations
               tl.cex = 1/par("cex"),
               cl.cex = 1/par("cex"),
               title = "Pearson COAD DBGS3",
               cex.main = 1.4/par("cex"),
               mar=c(5.1,4.1,4.1,2.1)
               )
par(cex = cex.before)
dev.off()




