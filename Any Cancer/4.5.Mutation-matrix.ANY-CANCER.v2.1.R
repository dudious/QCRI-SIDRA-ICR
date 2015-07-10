
#################################################
## Creates the (Sample by Gene) Mutation Matrix
##    
## Input: 
## ./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata
## ./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv
##
## Output: Table (Row = Sample Name, Column = Gene, Cell = NA/Cluster)
## 
## Modify: Cancerset, Geneset
## Author: SA,WHX
#####################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

## Parameters
Cancerset   = "COAD"
Geneset     = "DBGS3.FLTR"
matrix.type = "Any"         # Alterantives "Any" , "Missense"

##load data
## Read the mutation .maf file 
load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
muts = maf.merged.table
muts$sample.name = substr(muts$Tumor_Sample_Barcode, 1, 12)
## Read the gene list (373 genes)
genes.list = read.table("./2 DATA/Frequently mutated cancer genes.csv")
## Load the mutation variation data
load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".VariationTables.RData"))
## RNASeq clustering
Consensus.class = read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class = Consensus.class[,-1]
colnames (Consensus.class) = c("Patient_ID","Cluster")
rownames(Consensus.class) = Consensus.class[,1]
cluster.assignment = Consensus.class[which(rownames(Consensus.class) %in% muts$sample.name),] # drop samples without any mutation data

## Merge the cluster info and Remove the mutations with samples not having cluster information
muts$cluster = cluster.assignment$Cluster[match(muts$sample.name, as.character(cluster.assignment$Patient_ID))]
muts =  (muts[-which(is.na(muts$cluster)), ])

## Pick the Missense Mutations only
if (matrix.type =="Missense") {muts = muts[which(muts$Variant_Classification=="Missense_Mutation"), ]}

## Get the unique list of genes and samples
all.genes = unique(as.character(muts$Hugo_Symbol))
all.samples = unique(substr(muts$Tumor_Sample_Barcode, 1, 12))

## Create the table with col = genes and row = samples
genes.mutations = data.frame(matrix(ncol=length(all.genes), nrow=length(all.samples)))
rownames(genes.mutations) = all.samples
colnames(genes.mutations) = all.genes

## Fill the table
for(i in 1:nrow(muts)){
  sample.name = as.character(muts$sample.name[i])
  gene = as.character(muts$Hugo_Symbol[i])
  genes.mutations[sample.name,gene] = 1
  #genes.mutations[sample,gene] =  cluster.assignment$cluster[which(cluster.assignment$sample.name==sample.name)]
}

## Pick the genes from the provided (373 genes)list or significant variation File
genes.mutations.373genes = genes.mutations[,colnames(genes.mutations) %in% as.character(genes.list$V1)]
genes.mutations.low = genes.mutations[,colnames(genes.mutations) %in% as.character(low.significant.variation.table$Gene)]
genes.mutations.high = genes.mutations[,colnames(genes.mutations) %in% as.character(high.significant.variation.table$Gene)]
save (genes.mutations,genes.mutations.373genes,genes.mutations.low,genes.mutations.high,
      file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))

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


