#################################################################
###
### This Script creates a subset RNAseq AND/OR MA data 
### for a selected set  genes.(Gene_selection_XXX.txt)
### source data :
### "./2 DATA/TCGA ",Cancerset," MA/",Cancerset,".MA.TCGA.ASSEMBLER.CLEANED.RData" 
### "./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".NORMALIZED.LOG2.RData"
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
DL.Method = "BIOLINKS" #Choose "ASSEMBLER" or "BIOLINKS"
Cancerset <- "BRCA"
Geneset <- "DBGS3"
Genedatabase <- "Gene_selection_v2.6.txt"
MA.Data <- "NO"

# Load data
gene.list <- read.csv (paste0("./2 DATA/SUBSETS/",Genedatabase))                                 # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,Geneset]==1),1])

# RNAseq
## load data
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".NORMALIZED.LOG2.RData"))


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
dir.create(paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/"), showWarnings = FALSE)
save (RNASeq.subset,file=paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.subset.",Geneset,".RData"))    #adjust output file names here !!!!!

# Micro Array
if (MA.Data == "YES"){
  ## load data
  load (paste0("./2 DATA/TCGA MA/",Cancerset,"/",Cancerset,".MA.TCGA.ASSEMBLER.agilent.CLEANED.Rdata"))                                # no MA data for ",Cancerset,"
  
  # check availabilety of the genes in the dataset
  available.genes.MA <- gene.list.selected[which(gene.list.selected %in% colnames(agilentData))]
  unavailable.genes.MA <- gene.list.selected[-which(gene.list.selected %in% colnames(agilentData))]
  
  ## Subset data
  MA.subset <- agilentData[,available.genes.MA]
  
  ## report
  print (paste0("Geneset Selected : ",Geneset))
  print ("Genes selected for MA : ")
  print (available.genes.MA)
  print ("Genes missing for MA :")
  print (unavailable.genes.MA)
  
  # save subsetted data
  dir.create(paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/"), showWarnings = FALSE)
  save (MA.subset,file=paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".MA.subset.",Geneset,".RData"))                #adjust output file names here !!!!!
  
}
