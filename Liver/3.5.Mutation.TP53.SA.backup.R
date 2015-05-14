#########################
## Script to perform TP53 frequency analysis and barplot
## Input: Mutation .maf file, and the cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
##         Number of clusters (num.clusters)
##         Paths to mutation file, cluster assignment file, and output filename
## 
######


## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("ggplot2", "plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")

## Type of cancer
cancer = "SKCM"

## Read the mutation frequency file 
mut.freq <- read.csv("./3 ANALISYS/Mutations/LIHC/Mutations.TCGA.LIHC.Gene.by.Cluster.csv")

## Read the mutation .maf file (Melanome - BI Automated Mutation Calling Dataset)
muts.data = read.delim("./2 DATA/TCGA Mutations/LIHC/Somatic_Mutations/BCM__Mixed_DNASeq_curated/Level_2/hgsc.bcm.edu__Mixed_curated_DNA_sequencing_level2.maf")

## Read the cluster assignment file
cluster.assignment = read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LIHC/LIHC.TCGA.EDASeq.k7.ISGS.reps5000/LIHC.TCGA.EDASeq.k7.ISGS.reps5000.k=4.consensusClass.ICR.csv")
cluster.assignment = cluster.assignment[,-1]
colnames(cluster.assignment)= c("sample.name", "cluster")

## Pick the gene, mutation type, and sample name columns
muts = data.frame(gene = muts.data$Hugo_Symbol, mut.type = muts.data$Variant_Classification, sample.name = muts.data$Tumor_Sample_Barcode)
muts$sample.name = substr(muts$sample.name, 1, 12) #323324 mutations, 362 samples

## Ines RNASeq Clustering k = 4
num.clusters = 4
clusters = rep(paste0("ICR", 1:num.clusters))
cluster.assignment = data.frame(cluster.assignment) #469 samples

## Add the cluster assignment to the mutation table
muts$cluster = rep(NA, nrow(muts))
muts$cluster = cluster.assignment$cluster[match(muts$sample.name, cluster.assignment$sample.name)]

## Remove samples with no cluster information (1067 mutations and 2 samples)
muts =  (muts[-which(is.na(muts$cluster)), ]) #322257 mutations, 360 samples

## Count the number of TP53 genes per cluster
muts.TP53 = muts[which(muts$gene=="TP53" ), ] ## 67 mutations
muts.TP53 = muts.TP53[which(!duplicated(muts.TP53$sample.name)), ]  #56 unique samples, 11 duplicate gene,sample pair (different type of mutation)
TP53.cluster.count = as.data.frame(table(muts.TP53$gene=="TP53", muts.TP53$cluster))

## Count the number of samples in the mutation table per cluster (N)
muts.uniquesamples = muts[which(!duplicated(muts$sample.name)), ] 
sample.cluster.count = as.data.frame(table(muts.uniquesamples$cluster))

## Create data frame for plotting the percentage of TP53 patients per cluster
mut.df = data.frame(clusters, "muts" = rep("TP53", num.clusters))
mut.df$counts = TP53.cluster.count$Freq
mut.df$num.samples = sample.cluster.count$Freq
mut.df$percent = mut.df$counts/mut.df$num.samples*100
mut.df$label = paste0(sprintf("%.0f", mut.df$percent), "%")

## Plot TP53 percetage barplot
png(paste("./4 FIGURES/Mutation Plots/TP53.", cancer, ".k",num.clusters, ".png", sep=""))
cluster.order = rev(clusters)
gg = ggplot(mut.df, aes(x = clusters, y = percent, fill = muts))  +  geom_bar(stat = "identity", width = 0.8, position="dodge") +    geom_text( aes(y = percent+2, label = label), size = 7, position=position_dodge(width = 0.8)) 
gg = gg + scale_x_discrete(limits = cluster.order) + xlab("Clusters") + ylab("Percentage") + theme_bw()
gg = gg+ theme(legend.position="none", strip.text.x = element_text(size = 14)) + scale_fill_manual(values=c("gold", "lightsalmon4")) + ggtitle(("All Mutations - TP53"))
print(gg)
dev.off()
