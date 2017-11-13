# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
# Dependencies
required.packages <- c("beepr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("beepr")

## Parameters
Cancerset      = "OV"
Geneset        = "DBGS3.FLTR"
BRCA.Filter    = "BSF2"
matrix.type    = "NonSilent"         # Alterantives "Any" , "Missense", "NonSilent"
IMS.filter     = "All"           # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
selected.genes = c("TP53","MAP2K4","MAP3K1","CTCF","FCGBP")

#### creat submatrixes

load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))

## Read the gene list (373 genes)
genes.list = read.table("./2 DATA/Frequently mutated cancer genes.csv")

## Load the mutation variation data
load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".",matrix.type,".VariationTables.RData"))

## Pick the genes from the provided (373 genes)list or significant variation File
genes.mutations.373genes = genes.mutations[,colnames(genes.mutations) %in% as.character(genes.list$V1)]
genes.mutations.low = genes.mutations[,colnames(genes.mutations) %in% as.character(low.significant.variation.table$Gene)]
genes.mutations.high = genes.mutations[,colnames(genes.mutations) %in% as.character(high.significant.variation.table$Gene)]
genes.mutations.auto = genes.mutations[,colnames(genes.mutations) %in% as.character(auto.significant.variation.table$Gene)]
genes.mutations.selected = genes.mutations[,colnames(genes.mutations) %in% selected.genes]
genes.mutations.dbtest = genes.mutations[,colnames(genes.mutations) %in% db.test.significant.variation.table$Gene]
genes.mutations.dbtest.strict = genes.mutations[,colnames(genes.mutations) %in% db.test.strict.significant.variation.table$Gene]
genes.mutations.chisqr = genes.mutations[,colnames(genes.mutations) %in% chisq.significant.variation.table$Gene]
save (genes.mutations,
      genes.mutations.373genes,
      genes.mutations.low,
      genes.mutations.high,
      genes.mutations.auto,
      genes.mutations.dbtest,
      genes.mutations.dbtest.strict,
      genes.mutations.chisqr,
      genes.mutations.selected,
      file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.sub.Matrixes.",matrix.type,".Rdata"))

#### Get the number of genes that are mutated in more than 2,3, 5% of the samples
## for each gene, compute the frequency and percentage of samples that have the gene mutated
#geneCounts = colSums(genes.mutations, na.rm=T)
#num.samples = length(all.samples)
#geneCounts.percent = (geneCounts/num.samples)*100

#percentages = c(2,3,5,10)
#num.genes = rep(0, length(percentages))
#for(i in 1:length(percentages)){
#num.genes[i] = length(which(geneCounts.percent>percentages[i]))
#}

## Table with the percentage and number of genes 
#stats = data.frame(percentages, num.genes)

## Look at that genes and their percentages
#View(skcm.geneCounts.percent[which(skcm.geneCounts.percent>2)])