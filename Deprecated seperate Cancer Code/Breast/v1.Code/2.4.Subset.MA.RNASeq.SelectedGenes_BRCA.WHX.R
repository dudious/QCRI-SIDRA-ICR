#################################################################
###
### This Script creates a subset of the MA and RNAseq data 
### for a selected set  genes.(Gene_selection_XXX.txt)
### source data :
### "./2 DATA/TCGA BRCA MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.RData" 
### "./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"
### "./2 DATA/SUBSETS/Gene_selection_xxx.txt" (SELECTED GENES)
### Results are saved in
### ./2 DATA/SUBSETS/
### File to use :
### "TCGA.BRCA.MA.subset.16G.RData"
### "TCGA.BRCA.RNASeq.subset.16G.RData"
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

# load data
load ("./2 DATA/TCGA MA/BRCA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.Rdata")                                # no MA data for BRCA
load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_INES_MA.txt")                             # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,"Selected_by_DB"]==1),1])

# check availabilety of the genes in the dataset
available.genes.MA <- gene.list.selected[which(gene.list.selected %in% colnames(agilentData))]
unavailable.genes.MA <- gene.list.selected[-which(gene.list.selected %in% colnames(agilentData))]
available.genes.RNAseq <- gene.list.selected[which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]
unavailable.genes.RNAseq <- gene.list.selected[-which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]

## Subset data
MA.subset <- agilentData[,available.genes.MA]
RNASeq.subset <- t(RNASeq.NORM_Log2[available.genes.RNAseq,])

## report
print ("Genes selected for MA : ")
print(available.genes.MA)
print ("Genes missing for MA :")
print(unavailable.genes.MA)
print ("Genes selected for RNASeq : ")
print(available.genes.RNAseq)
print ("Genes missing for RNASeq :")
print(unavailable.genes.RNAseq)

# save subsetted data
save (MA.subset,file="./2 DATA/SUBSETS/BRCA/TCGA.BRCA.MA.subset.ISGS.RData")                #adjust output file names here !!!!!
save (RNASeq.subset,file="./2 DATA/SUBSETS/BRCA/TCGA.BRCA.RNASeq.subset.ISGS.RData")    #adjust output file names here !!!!!

