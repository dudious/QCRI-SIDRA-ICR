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
cancer = "LIHC"

## Ines RNASeq Clustering k = 4
num.clusters = 4
clusters = rep(paste0("ICR", 1:num.clusters))

## Read the mutation frequency file 
load ("./3 ANALISYS/Mutations/LIHC/Mutation.Data.Frequencies.RDATA")
#load ("./3 ANALISYS/Mutations/LIHC/Mutation.Data.split.RDATA")

#Prepare Data for Boxplots

numMuts.SGall <- data.frame(count=Mutation.Frequency.Patient$Freq.All,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "All")

numMuts.Missense  <- data.frame(count=Mutation.Frequency.Patient$Freq.Missense,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Missense")
numMuts.Nonsense  <- data.frame(count=Mutation.Frequency.Patient$Freq.Nonsense,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Nonsense")
numMuts.Silent    <- data.frame(count=Mutation.Frequency.Patient$Freq.Silent,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Silent")
numMuts.Other     <- data.frame(count=Mutation.Frequency.Patient$Freq.Other,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Other")
numMuts.SGtype    <- rbind(numMuts.Missense,numMuts.Silent,numMuts.Nonsense,numMuts.Other)
rownames(numMuts.SGtype) <- NULL
colnames(numMuts.SGtype) = c( "count", "cluster", "mut.type")

# Combine
numMuts.SGall = rbind(numMuts.SGall, numMuts.SGtype)

meds <- ddply(numMuts.SGall, .(mut.type, cluster), summarize, med = median(count)) ## median
mean.n <- function(x){ return(c(y = 0 , label = round(mean(x),2))) } ## mean

png("./4 FIGURES/Mutation Plots/Mutations.TCGA.LIHC.All.SA.WHX.png", height = 1000, width= 1000)   #set filename
cluster.order = c("ICR4", "ICR3", "ICR2", "ICR1")
colors = c("blue", "green", "orange", "red")
gg = ggplot(numMuts.SGall, aes(cluster, count, fill=cluster)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  geom_jitter(position=position_jitter(width=0.1,height=0.1))
#, aes(color="lightgray")) 
gg = gg + ylab("Number of mutations per sample") +
  scale_x_discrete(limits=cluster.order) +
  facet_grid(.~mut.type,
             scales = "free",
             space="free") +
  xlab("Clusters") + theme_bw() 
gg = gg + scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = seq(0, 300, 100), limits = c(0, 250)) 
gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                legend.position = "none",
                axis.text.x = element_text(size = 12, vjust=1),
                axis.title.x = element_text(size = 18, vjust = -1),
                axis.text.y = element_text(size = 18, vjust=1),
                axis.title.y = element_text(size = 18, vjust = 1))
gg = gg + geom_text(data = meds,
                    aes(y = 0, label = round(med,2)),
                    size = 7, vjust = 1.2)
gg = gg + stat_summary(fun.data = mean.n,
                       geom = "text",
                       fun.y = mean,
                       colour = "black",
                       vjust = 3.7,
                       size = 7)
print(gg)
dev.off()

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
