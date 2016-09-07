#################################################################
###
### This Script creates a subset of the MA and RNAseq data 
### for a selected set  genes.(Gene_selection_XXX.txt)
### source data :
### "./2 DATA/TCGA ",Cancerset," MA/",Cancerset,".MA.TCGA.ASSEMBLER.CLEANED.RData" 
### "./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"
### "./2 DATA/SUBSETS/Gene_selection_xxx.txt" (SELECTED GENES)
### Results are saved in
### ./2 DATA/SUBSETS/
### File to use :
### "TCGA.",Cancerset,".MA.subset.16G.RData"
### "TCGA.",Cancerset,".RNASeq.subset.16G.RData"
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR")

# Parameters
Cancerset <- "BRCA"
Geneset <- "CHKPNT_IL"
Genedatabase <- "Gene_selection_v2.7.txt"
DL.method <- "ASSEMBLER"

# Load data
gene.list <- read.csv (paste0("./2 DATA/SUBSETS/",Genedatabase))                                 # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,Geneset]==1),1])

# RNAseq
## load data
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.method,".NORMALIZED.LOG2.RData"))

# check availabilety of the genes in the dataset
available.genes.RNAseq <- gene.list.selected[which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]
unavailable.genes.RNAseq <- gene.list.selected[-which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]

## Subset data
RNASeq.subset <- t(RNASeq.NORM_Log2[available.genes.RNAseq,])

## report
print (paste0("Geneset Selected : ",Geneset))
print ("Genes selected for RNASeq : ")
print (available.genes.RNAseq)
print ("Genes missing for RNASeq :")
print (unavailable.genes.RNAseq)

# save subsetted data
dir.create(paste0("./2 DATA/SUBSETS/",Cancerset,"/"), showWarnings = FALSE)
save (RNASeq.subset,file=paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".for_CHKPNT_heatmap.RData"))    #adjust output file names here !!!!!

