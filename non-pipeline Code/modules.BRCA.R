# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

# methode chausebell
library(tmod)

load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/ICR4vs1_UPs.RDATA")
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")

fg = rownames(ICR4vs1_UPs)
bg = unique (rownames(RNASeq.NORM_Log2))
tmod.result = tmodHGtest(fg=fg, bg=bg)
tmodCERNOtest(fg)

hgEnrichmentPlot (fg,bg,"LI.M7.0" )
evidencePlot(fg, "LI.M7.0")

res.1 <- tmodUtest( fg)

tmodPanelPlot(res.1, text.cex=0.8)

#write.csv(ICR4vs1_UPs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/")