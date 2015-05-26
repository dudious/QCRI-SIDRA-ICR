
#################################################
## Creates the (Sample by Gene) Mutation Table
##    
## Input: Mutation .maf file, Cluster Assignment file (format: sample.name, cluster)
## Output: Table (Row = Sample Name, Column = Gene, Cell = NA/Cluster)
##
## Author: SA
#####################################

## Parameters
mutation.type = "Any"       # Select Any / Missense
Gene.selection = "List"     # Select List / min3pct / min5pct / min10pct
Matrix.Type = "Binary"      # Select Binary / Count / Cluster_Assignment

## load Data
# Read the mutation .maf file (LIHC)
muts = read.csv("./2 DATA/TCGA Mutations/LIHC/Somatic_Mutations/BCM__Mixed_DNASeq_curated/Level_2/hgsc.bcm.edu__Mixed_curated_DNA_sequencing_level2.maf", sep="\t" )
muts$sample.name = substr(muts$Tumor_Sample_Barcode, 1, 12)
# Ines RNASeq clustering
cluster.assignment = read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LIHC/LIHC.TCGA.EDASeq.k7.ISGS.reps5000/LIHC.TCGA.EDASeq.k7.ISGS.reps5000.k=4.consensusClass.ICR.csv", header =TRUE)
cluster.assignment$X = NULL
# Read the gene list (373 genes)
genes.list = read.csv("./2 DATA/Frequently mutated cancer genes.csv",header =FALSE)
colnames (genes.list) = "gene.name"

## Merge the cluster info and Remove the mutations with samples not having cluster information
muts$cluster = cluster.assignment$Group[match(muts$sample.name, cluster.assignment$PatientID)]
muts =  (muts[-which(is.na(muts$cluster)), ])

## Pick the Missense Mutations only
#mut = mut[which(mut$Variant_Classification=="Missense_Mutation"), ]

## Get the unique list of genes and samples
all.genes = unique(as.character(muts$Hugo_Symbol))
all.samples = unique(muts$sample.name)

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
genes.mutations[is.na(genes.mutations)] <- 0


## Pick the genes from the list provided (373 genes)
#genes.mutations = genes.mutations[,colnames(genes.mutations) %in% as.character(genes.list$V1)]

write.csv(genes.mutations, "./3 ANALISYS/Mutations/LIHC/Mutation.heatmap.matrix.csv", row.names = TRUE, quote = FALSE)

