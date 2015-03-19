#################################################################
###
### This Script creates a subset of RNAseq data 
### for the N selected genes.
### source data :
### "./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.RData"
### "./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData"
### "./2 DATA/imunogenes.whx.TXT" (SELECTED GENES)
### Results are saved in
### ./2 DATA/SUBSETS/
### File to use :
### "MA.subset.16G.RData"
### "RNASeq.subset.16G.RData"
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

# load data
load ("./2 DATA/METABRIC/RNASEQ.DATA1.RData")
RNASEQ.DATA.1 <- RNASEQ.DATA
load ("./2 DATA/METABRIC/RNASEQ.DATA2.RData")
RNASEQ.DATA.2 <- RNASEQ.DATA
rm(RNASEQ.DATA)
gene.list <- read.csv ("./2 DATA/SUBSETS/Gene_selection_INES_MA.txt")                          # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,"Selected_by_DB"]==1),1])

# Subset RNASEQ.DATA.1
dim (RNASEQ.DATA.1) #997 25159
print ("Gene selection or Micro Array")
RNASEQ.DATA.1.subset <- as.data.frame(RNASEQ.DATA.1[,1]) #subset(RNASEQ.DATA.1,select=gene.list.selected[1])
for (i in 1:length(gene.list.selected))
  {
  if (gene.list.selected[i] %in% colnames(RNASEQ.DATA.1))
    {
    RNASEQ.DATA.1.subset <- cbind (RNASEQ.DATA.1.subset,subset (RNASEQ.DATA.1,select=gene.list.selected[i]))
    cat(i,"Gene added :",gene.list.selected[i],"\n")  
    }
  else
    { cat (i,"Gene NOT added :",gene.list.selected[i],"\n") }
  }
RNASEQ.DATA.1.subset[,1] <- NULL
dim (RNASEQ.DATA.1.subset)   #997  16

# Subset RNASEQ.DATA.2
dim (RNASEQ.DATA.2) #995 25159
print ("Gene selection or Micro Array")
RNASEQ.DATA.2.subset <- as.data.frame(RNASEQ.DATA.2[,1]) #subset(RNASEQ.DATA.2,select=gene.list.selected[1])
for (i in 1:length(gene.list.selected))
{
  if (gene.list.selected[i] %in% colnames(RNASEQ.DATA.2))
  {
    RNASEQ.DATA.2.subset <- cbind (RNASEQ.DATA.2.subset,subset (RNASEQ.DATA.2,select=gene.list.selected[i]))
    cat(i,"Gene added :",gene.list.selected[i],"\n")  
  }
  else
  { cat (i,"Gene NOT added :",gene.list.selected[i],"\n") }
}
RNASEQ.DATA.2.subset[,1] <- NULL
dim (RNASEQ.DATA.2.subset)   #995  16

# save subsetted data
save (RNASEQ.DATA.1.subset,file="./2 DATA/SUBSETS/METABRIC.RNASEQ.DATA.1.genesubset.16G.I.RData")    #adjust output file names here !!!!!
save (RNASEQ.DATA.2.subset,file="./2 DATA/SUBSETS/METABRIC.RNASEQ.DATA.2.genesubset.16G.I.RData")    #adjust output file names here !!!!!

