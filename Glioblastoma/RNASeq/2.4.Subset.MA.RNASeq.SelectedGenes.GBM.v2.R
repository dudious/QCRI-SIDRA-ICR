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
setwd("~/Dropbox/BREAST_QATAR")

# Parameters
Cancerset <- "GBM"
Geneset <- "DBGS1"

# load data
load (paste0("./2 DATA/TCGA MA/",Cancerset,"/",Cancerset,".MA.TCGA.ASSEMBLER.CLEANED.Rdata"))                                # no MA data for ",Cancerset,"
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_v2.2.txt")                                 # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,Geneset]==1),1])


# check availabilety of the genes in the dataset
available.genes.MA <- gene.list.selected[which(gene.list.selected %in% colnames(AffyData))]
unavailable.genes.MA <- gene.list.selected[-which(gene.list.selected %in% colnames(AffyData))]
available.genes.RNAseq <- gene.list.selected[which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]
unavailable.genes.RNAseq <- gene.list.selected[-which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]

## Subset data
MA.subset <- AffyData[,available.genes.MA]
RNASeq.subset <- t(RNASeq.NORM_Log2[available.genes.RNAseq,])

## report
print (paste0("Geneset Selected : ",Geneset))
print ("Genes selected for MA : ")
print (available.genes.MA)
print ("Genes missing for MA :")
print (unavailable.genes.MA)
print ("Genes selected for RNASeq : ")
print (available.genes.RNAseq)
print ("Genes missing for RNASeq :")
print (unavailable.genes.RNAseq)

# save subsetted data
dir.create(paste0("./2 DATA/SUBSETS/",Cancerset,"/"), showWarnings = FALSE)
save (MA.subset,file=paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".MA.subset.",Geneset,".RData"))                #adjust output file names here !!!!!
save (RNASeq.subset,file=paste0("./2 DATA/SUBSETS/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".RData"))    #adjust output file names here !!!!!


