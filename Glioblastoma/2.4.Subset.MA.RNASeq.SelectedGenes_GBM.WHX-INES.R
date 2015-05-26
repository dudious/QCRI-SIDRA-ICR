#################################################################
###
### This Script creates a subset of the MA and RNAseq data 
### for a selected set  genes.(Gene_selection_XXX.txt)
### source data :
### "./2 DATA/TCGA GBM MA/GBM.MA.TCGA.ASSEMBLER.CLEANED.RData" (NO MA DATA FOR GBM in 2015)
### "./2 DATA/TCGA RNAseq/RNASeq_GBM_EDASeq/GBM.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"
### "./2 DATA/SUBSETS/Gene_selection_xxx.txt" (SELECTED GENES)
### Results are saved in
### ./2 DATA/SUBSETS/
### File to use :
### "TCGA.GBM.MA.subset.16G.RData"
### "TCGA.GBM.RNASeq.subset.16G.RData"
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

# load data
#load ("./2 DATA/TCGA MA/TCGA_GBM_ASSEMBLER/GBM.MA.Affy.TCGA.ASSEMBLER.QN.NORMALIZED.RData")                                # no MA data for GBM
load ("./2 DATA/TCGA RNAseq/RNASeq_GBM_EDASeq/GBM.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_INES_MA.txt")                             # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,"Selected_by_DB"]==1),1])

# check availabilety of the genes in the dataset
#available.genes.MA <- gene.list.selected[which(gene.list.selected %in% colnames(AffyData.NORM))]
#unavailable.genes.MA <- gene.list.selected[-which(gene.list.selected %in% colnames(AffyData.NORM))]#CD8A
available.genes.RNAseq <- gene.list.selected[which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]
unavailable.genes.RNAseq <- gene.list.selected[-which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]

## Subset data
# AffyData.NORM<-t(AffyData.NORM)
# MA.subset <- AffyData.NORM[available.genes.MA,]
# MA.subset<-t(MA.subset)
RNASeq.subset <- t(RNASeq.NORM_Log2[available.genes.RNAseq,])

## report
#print (paste0("Genes selected for MicroArray :",available.genes.MA))
#print (paste0("Genes missing for MicroArray  :",unavailable.genes.MA))
print ("Genes selected for RNASeq : ")
print(available.genes.RNAseq)
print ("Genes missing for RNASeq :")
print(unavailable.genes.RNAseq)

# save subsetted data
#save (MA.subset,file="./2 DATA/SUBSETS/GBM/MA.subset.ISGS.RData")                #adjust output file names here !!!!!
save (RNASeq.subset,file="./2 DATA/SUBSETS/GBM/TCGA.GBM.RNASeq.subset.ISGS.RData")    #adjust output file names here !!!!!

