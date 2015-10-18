rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
required.packages <- c("devtools","stringi")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(devtools)
install_github("dakl/oncoprint")
library(oncoprint)
load("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.All.DBGS3.FLTR.Mutation.Matrixes.NonSilent.oncoplot.Rdata")

selected.data <- t(merged.matrix[,colnames(genes.mutations.dbtest)])
dim(selected.data)
# vertical x-labels
vert_x <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
png ("./4 FIGURES/oncoprint.v1.png",height=1000, width = 1500)
oncoprint(selected.data) + coord_fixed() + vert_x
dev.off()
