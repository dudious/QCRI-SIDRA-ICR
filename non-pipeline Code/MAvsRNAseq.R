# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")
library(ggplot2)
library(reshape)

#load data
load ("./2 DATA/SUBSETS/BRCA/TCGA.BRCA.MA.subset.DBGS3.RData")
load ("./2 DATA/SUBSETS/BRCA/TCGA.BRCA.RNASeq.subset.DBGS3.RData")

#drop non overlapping data and re-order alike
RNASeq.subset.MAset <- RNASeq.subset[rownames(RNASeq.subset) %in% rownames(MA.subset),]
MA.subset.RNAseqset <- MA.subset[rownames(MA.subset) %in% rownames(RNASeq.subset),]
RNASeq.subset.MAset <- RNASeq.subset.MAset[,colnames(MA.subset.RNAseqset)]
RNASeq.subset.MAset <- RNASeq.subset.MAset[rownames(MA.subset.RNAseqset),]

#ranked matrixes
RNASeq.subset.MAset.ranked <- as.matrix(t(apply(RNASeq.subset.MAset, 1, rank, ties.method='min')))
MA.subset.RNAseqset.ranked <- as.matrix(t(apply(MA.subset.RNAseqset, 1, rank, ties.method='min')))

#imunoscore (rowmeans)
RNASeq.subset.MAset <- cbind(RNASeq.subset.MAset,rowMeans(RNASeq.subset.MAset))
colnames(RNASeq.subset.MAset)[ncol(RNASeq.subset.MAset)] <- c("IGS")
MA.subset.RNAseqset <- cbind(MA.subset.RNAseqset,rowMeans(MA.subset.RNAseqset))
colnames(MA.subset.RNAseqset)[ncol(MA.subset.RNAseqset)] <- c("IGS")

#plot
for (gene in colnames(RNASeq.subset.MAset)){
  print(gene)
  
  df<-data.frame(RNASeq.subset.MAset[,gene],MA.subset.RNAseqset[,gene])
  colnames(df) <- c("RNASEQ","MA")
  mypath <- file.path(paste0("./4 FIGURES/CORRELATION/RNASeq_vs_MA/scatterplot.",gene,".png"))
  png(file=mypath)
    plot <- ggplot(df, aes(x=RNASEQ, y=MA)) + geom_point() + stat_smooth(method=lm)
    print(plot)
  dev.off()
  
  melted.df<-melt(as.matrix(df))
  colnames(melted.df)<-c("ID","tech","value")
  mypath <- file.path(paste0("./4 FIGURES/CORRELATION/RNASeq_vs_MA/histogram.",gene,".png"))
  png(file=mypath)
   plot <- ggplot(melted.df, aes(x=value)) + geom_histogram(binwidth=0.25, fill="white", colour="black") + facet_grid(tech ~ .)
   print(plot)
  dev.off()
  
  ranked.df <- data.frame(RNASeq.subset.MAset.ranked[,gene],MA.subset.RNAseqset.ranked[,gene])
  colnames(ranked.df) <- c("RNASEQ","MA")
  mypath <- file.path(paste0("./4 FIGURES/CORRELATION/RNASeq_vs_MA/rankplot.",gene,".png"))
  png(file=mypath)
   plot <- ggplot(ranked.df, aes(x=RNASEQ, y=MA)) +
    geom_point(position="jitter") +
    stat_smooth(method=lm) +
    scale_y_continuous(limits=c(0,18),breaks=seq(0,18,3)) +
    scale_x_continuous(limits=c(0,18),breaks=seq(0,18,3))
   print(plot)
  dev.off()
  
  assign(gene,df)
}




