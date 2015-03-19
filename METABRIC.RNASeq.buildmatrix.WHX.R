#############################################################################
### Script for Probe to gene colapse of RNASeq Probe level data
### 
### Input data :
### ./2 DATA/METABRIC/FROM GABRIELE ZOPPOLI/GeneExpression_METABRIC2012_a.RData
### Data is saved :
### ./2 DATA/METABRIC/RNASEQ.DATA1.RData
###
#############################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
# Dependencies
required.packages <- c("plyr","data.table")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("plyr")
library("data.table")

# Load Data
load ("./2 DATA/METABRIC/FROM GABRIELE ZOPPOLI/GeneExpression_METABRIC2012_b.RData")           # load file a/b 

# append gene name
dim (GeneExpression)                                                                           # 997/995 x 48803
RNASEQ.T <- t(GeneExpression)                                                                  # transpose data
dim (RNASEQ.T)                                                                                 # 48803 x 997/995

RNASEQ.T.GN <- merge (RNASEQ.T,probeset_annotation[,2,drop=FALSE],by="row.names",all=TRUE,allow=TRUE)   # Append Gene name to RNASeq Data
rownames (RNASEQ.T.GN) <- RNASEQ.T.GN$Row.names
RNASEQ.T <- RNASEQ.T.GN[,-1]

#RNASEQ.T <- cbind(RNASEQ.T,rownames(RNASEQ.T))
#probeset_annotation$rn <- rownames(probeset_annotation)
#RNASEQ.T.GN2 <- join (RNASEQ.T,probeset_annotation[,2,drop=FALSE],by="rn",match="all",type="left")

RNASEQ.T <- RNASEQ.T [,c(996,1:995)]                                                           # reorder genename column
dim (RNASEQ.T)                                                                                 # 48803 x 998
RNASEQ.T <- RNASEQ.T [with(RNASEQ.T,order(GeneName)),]                                         # order by genename
RNASEQ.T <- RNASEQ.T[!(RNASEQ.T$GeneName == ""), ]                                             # drop rows with no gene name
dim (RNASEQ.T)                                                                                 # 36157 x 998

# collapse to gene level by average of probes

#RNASEQ.T.AVG <- aggregate(x = RNASEQ.T[, 2:ncol(RNASEQ.T)], by = list(GeneName = RNASEQ.T$GeneName), FUN = "mean", na.rm = TRUE) # average by gene
#dim(RNASEQ.AVG)                                                                                # 25159 x 998

RNASEQ.T.dt <- as.data.table(RNASEQ.T)                                                         # Much faste way to average by gene using DATA.TABLE
setkey(RNASEQ.T.dt,GeneName)
RNASEQ.T.dt.AVG <- RNASEQ.T.dt[, lapply(.SD,mean), by=GeneName]
dim(RNASEQ.T.dt.AVG)                                                                           # 25159 x 998

# transpose and save

RNASEQ.DATA <- t(RNASEQ.T.dt.AVG)                                                              # transpose data
dim(RNASEQ.DATA)                                                                               # 998 x 25159
colnames(RNASEQ.DATA) <- RNASEQ.DATA[1,]                                                       # set column names to gene names
RNASEQ.DATA <- RNASEQ.DATA[-1,]                                                                # remove gene name row
mode(RNASEQ.DATA) <- "numeric"                                                                 # convert matrix to numeric
save (RNASEQ.DATA,file="./2 DATA/METABRIC/RNASEQ.DATA2.RData")                                 # save
