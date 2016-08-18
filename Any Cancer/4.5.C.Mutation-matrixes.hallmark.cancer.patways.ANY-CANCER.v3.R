# Setup environment
rm(list=ls())
#setwd("~/Dropbox (Research team)/BREAST_QATAR/")
setwd("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
# Dependencies
required.packages <- c("beepr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("beepr")

## Parameters
Cancerset      = "BRCA"
Geneset        = "DBGS3.FLTR"
BRCA.Filter    = "PCF"
matrix.type    = "NonSilent"         # Alterantives "Any" , "Missense", "NonSilent"
IMS.filter     = "All"           # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"

## Load Data
if (Cancerset =="BRCA") {
  Cancerset <- paste0(Cancerset,".",BRCA.Filter)
}
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))
#load hallmark pathways data
load("./2 Data/Hallmark Cancer Pathways/cancer.halmark.pathways.R")
#add barier genes
hallmark.pathways$BARRIER_GENES <- c("FLG","TACSTD2","DSC3","DST","DSP","PPL","PKP3","JUP")

hallmark.pathways.Zscores <- data.frame(matrix(ncol=length(hallmark.pathways), nrow=nrow(genes.mutations)))
rownames(hallmark.pathways.Zscores) = rownames(genes.mutations)
colnames(hallmark.pathways.Zscores) = names(hallmark.pathways)

for(i in 1:length(hallmark.pathways)){
pathway.genes <- unlist(hallmark.pathways[i],use.names = FALSE)
present.genes <- pathway.genes[which(pathway.genes %in% colnames(genes.mutations))]
assign(paste(names(hallmark.pathways[i]),".Mutation.matrix", sep=""),genes.mutations[,present.genes])
hallmark.pathways.Zscores[,names(hallmark.pathways[i])] <- rowSums(genes.mutations[,present.genes],na.rm = TRUE)
}

## Save Data
save (hallmark.pathways.Zscores ,file =  paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Hallmark.cancer.pathway.Mutation.Matrix.",matrix.type,".Rdata"))

