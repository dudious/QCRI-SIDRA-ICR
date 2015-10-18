
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
# Dependencies
required.packages <- c("beepr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("beepr")

## Parameters
Cancerset      = "BRCA"
Geneset        = "DBGS3.FLTR"
BRCA.Filter    = "BSF2"
matrix.type    = "NonSilent"         # Alterantives "Any" , "Missense", "NonSilent"
IMS.filter     = "All"           # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
selected.genes = c("TP53","MAP2K4","MAP3K1","CTCF","FCGBP")

##load data
## Read the mutation .maf file and cluster assignments
if (Cancerset %in% c("COAD","READ","UCEC")) {
  #GA data
  Cancerset <- paste0(Cancerset,"-GA")
  load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
  maf.merged.table.GA <- maf.merged.table
  Consensus.class.GA <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class.GA <- Consensus.class.GA[,-1]
  colnames (Consensus.class.GA) <- c("Patient_ID","Cluster")
  rownames(Consensus.class.GA) <- Consensus.class.GA[,1]
  rm(maf.merged.table)
  Cancerset <- substring(Cancerset,1,4)
  #hiseq data
  Cancerset <- paste0(Cancerset,"-hiseq")
  load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
  maf.merged.table.hiseq <- maf.merged.table 
  Consensus.class.hiseq <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class.hiseq <- Consensus.class.hiseq[,-1]
  colnames (Consensus.class.hiseq) <- c("Patient_ID","Cluster")
  rownames(Consensus.class.hiseq) <- Consensus.class.hiseq[,1]
  rm(maf.merged.table)
  Cancerset <- substring(Cancerset,1,4)
  #merge GA-hiseq
  Consensus.class <- unique(rbind (Consensus.class.hiseq,Consensus.class.GA))
  maf.merged.table   <- unique(rbind (maf.merged.table.hiseq,maf.merged.table.GA))
} else {
  load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
  if (Cancerset == "BRCA"){
    if (substring(Geneset,7,10)=="FLTR"){
      Cancerset <- paste0(Cancerset,".",BRCA.Filter)
    }
  }
  Consensus.class = read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class = Consensus.class[,-1]
  colnames (Consensus.class) = c("Patient_ID","Cluster")
  rownames(Consensus.class) = Consensus.class[,1]
} 
muts = maf.merged.table
rm(maf.merged.table)
muts$sample.name = substr(muts$Tumor_Sample_Barcode, 1, 12)
cluster.assignment = Consensus.class[which(rownames(Consensus.class) %in% muts$sample.name),] # drop samples without any mutation data

## Load gistic data
gistic.ICR1 <- read.csv(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/OUTPUT/all_ThreshByGenesICR1.csv"),sep = ";")
gistic.ICR2 <- read.csv(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/OUTPUT/all_ThreshByGenesICR2.csv"),sep = ";")
gistic.ICR3 <- read.csv(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/OUTPUT/all_ThreshByGenesICR3.csv"),sep = ";")
gistic.ICR4 <- read.csv(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/OUTPUT/all_ThreshByGenesICR4.csv"),sep = ";")
gistic <- merge (gistic.ICR1[,-c(2,3)],gistic.ICR2[,-c(2,3)],by = "Gene.Symbol")
gistic <- merge (gistic,gistic.ICR3[,-c(2,3)],by = "Gene.Symbol")
gistic <- merge (gistic,gistic.ICR4[,-c(2,3)],by = "Gene.Symbol")
rownames(gistic) <- gistic$Gene.Symbol
gistic$Gene.Symbol <- NULL
gistic <- t(gistic)
gistic[gistic==1|gistic==2] <-"AMP"
gistic[gistic==-1|gistic==-2] <-"HOMDEL"
gistic[gistic==0] <- NA
rownames(gistic)<-gsub("\\.","-",rownames(gistic))

## Merge the cluster info and Remove the mutations with samples not having cluster information
muts$cluster = cluster.assignment$Cluster[match(muts$sample.name, as.character(cluster.assignment$Patient_ID))]
muts =  (muts[-which(is.na(muts$cluster)), ])

## Read the gene list (373 genes)
genes.list = read.table("./2 DATA/Frequently mutated cancer genes.csv")

## Load the mutation variation data
load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".",matrix.type,".VariationTables.RData"))

## Pick the Missense Mutations only
if (matrix.type =="Missense") {muts = muts[which(muts$Variant_Classification=="Missense_Mutation"), ]}

## Pick the NonSilent Mutations only
if (matrix.type =="NonSilent") {muts = muts[which(muts$Variant_Classification!="Silent"), ]}

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
  genes.mutations[sample.name,gene] = "MUT"
  #genes.mutations[sample,gene] =  cluster.assignment$cluster[which(cluster.assignment$sample.name==sample.name)]
}

##equalize matrixes
#samples
gistic <- as.data.frame(gistic)
gistic<-gistic[rownames(genes.mutations),]
#genes
overlap.genes <- colnames(gistic)[which(colnames(gistic) %in% colnames(genes.mutations))]
gistic.genes <- colnames(gistic)[-which(colnames(gistic) %in% colnames(genes.mutations))]
mutation.genes <- colnames(genes.mutations)[-which(colnames(genes.mutations) %in% colnames(gistic))]
temprow <- as.data.frame(matrix(c(rep.int(NA,length(mutation.genes))),nrow=nrow(gistic),ncol=length(mutation.genes)))
rownames(temprow) <- rownames(mutation.genes)
colnames(temprow) <-mutation.genes
gistic <- cbind(gistic,temprow)
temprow <- as.data.frame(matrix(c(rep.int(NA,length(gistic.genes))),nrow=nrow(genes.mutations),ncol=length(gistic.genes)))
rownames(temprow) <- rownames(gistic)
colnames(temprow) <-gistic.genes
genes.mutations<- cbind(genes.mutations,temprow)
#reorder gistic
gistic<-gistic[rownames(genes.mutations),colnames(genes.mutations)]
gistic[] <- lapply(gistic, as.character)

#merge gistic and mutation data
merged.matrix = as.data.frame(matrix(ncol=ncol(genes.mutations), nrow=nrow(genes.mutations)))
rownames(merged.matrix) = rownames(genes.mutations)
colnames(merged.matrix) = colnames(genes.mutations)
for (i in 1:nrow(gistic)) {
  merged.matrix[i,] <- paste0 (as.character(genes.mutations[i,]),";",as.character(gistic[i,]))
  print (i)
}
#cleanup NA data
merged.matrix <- gsub("NA;HOMDEL","HOMDEL",as.matrix(merged.matrix))
merged.matrix <- gsub("NA;AMP","AMP",as.matrix(merged.matrix))
merged.matrix <- gsub("NA;MUT","MUT",as.matrix(merged.matrix))
merged.matrix <- gsub("HOMDEL;NA","HOMDEL",as.matrix(merged.matrix))
merged.matrix <- gsub("AMP;NA","AMP",as.matrix(merged.matrix))
merged.matrix <- gsub("MUT;NA","MUT",as.matrix(merged.matrix))
merged.matrix[merged.matrix == "NA;NA"] <- NA

save (genes.mutations,merged.matrix,
      file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))

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
      merged.matrix,
      genes.mutations.373genes,
      genes.mutations.low,
      genes.mutations.high,
      genes.mutations.auto,
      genes.mutations.dbtest,
      genes.mutations.dbtest.strict,
      genes.mutations.chisqr,
      genes.mutations.selected,
      file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))

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


