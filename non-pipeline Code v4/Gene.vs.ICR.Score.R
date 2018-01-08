#########################
## Script to perform pecific gene barplot
## Input: cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
##         Number of clusters (num.clusters)
##         Paths to mutation file, cluster assignment file, and output filename
## 
######


## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")
#Dependencies
required.packages <- c("ggplot2", "plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")
library("heatmap3")
library("corrplot")

## Parameters
Cancerset <- "SKCM"           
GOF = c("HLA-E","HLA-A","HLA-DRA")
download.method = "TCGA_Assembler"  
ALL.HLA = c("HHLA1","HHLA2","HHLA3","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2",
            "HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DRB6","HLA-E","HLA-F","HLA-F-AS1","HLA-G","HLA-L")
code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/"

source(paste0(code_path, "R tools/heatmap.3.R"))

## Load Data
load (paste0(code_path, "Datalists/ICR_genes.RData"))
# Load normalised RNAseq data
Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancerset,"/RNASeqData")
if(Cancerset == "SKCM"){
  load(paste0(Cancer_path, "/", Cancerset, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
} else{
  load(paste0(Cancer_path, "/", Cancerset, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
}
# cluster assignment
load(paste0("./4_Analysis/TCGA_Assembler/",Cancerset,"/Clustering/",Cancerset,".TCGA_Assembler.EDASeq.ICR.reps5000/",Cancerset,"_ICR_cluster_assignment_k2-6.Rdata"))# select source data

#Select gene
RNASeq.Data <- as.data.frame(t(log(filtered.norm.RNAseqData[GOF,,drop=FALSE],2)))
HLA.heatmap <- as.data.frame(t(log(filtered.norm.RNAseqData[ALL.HLA,]+1,2)))
ICR.heatmap <- as.data.frame(t(log(filtered.norm.RNAseqData[ICR_genes,]+1,2)))

#Add Class to RNAseq data
RNASeq.Data$HML_Cluster <- factor(table_cluster_assignment$HML_cluster[match(rownames(RNASeq.Data),rownames(table_cluster_assignment))],levels = c("ICR Low","ICR Medium","ICR High"))
if (any(is.na(RNASeq.Data$HML_Cluster))){RNASeq.Data <-  RNASeq.Data[-which(is.na(RNASeq.Data$HML_Cluster)),]}
RNASeq.Data$ICR.score <-table_cluster_assignment$ICRscore[match(rownames(RNASeq.Data),rownames(table_cluster_assignment))]
colnames(RNASeq.Data) <- c("HLA.E","HLA.A","HLA.DRA","ICR.cluster","ICR.score")

#blots
GOF="HLA-E"
dev.new()
test.anova      = aov(HLA.E~ICR.cluster,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = ICR.cluster, y = HLA.E  )) +
       geom_boxplot() +
       labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "Cluster", y = "Log 2 Expression") +
       theme_classic()

corelation = round(cor(x = RNASeq.Data$HLA.E,y=RNASeq.Data$ICR.score,method = "spearman"),2)
dev.new()
ggplot(RNASeq.Data, aes(x = ICR.score, y = HLA.E  )) +
  geom_point() +
  labs(title = paste0(GOF," in ",Cancerset," (",corelation,"[Spearman])"), x = "ICS score", y = "Log 2 Expression") +
  geom_smooth(method='lm',formula=y~x) +
  theme_classic()

#blots
GOF="HLA-A"
dev.new()
test.anova      = aov(HLA.A~ICR.cluster,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = ICR.cluster, y = HLA.A  )) +
  geom_boxplot() +
  labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "Cluster", y = "Log 2 Expression") +
  theme_classic()

corelation = round(cor(x = RNASeq.Data$HLA.A,y=RNASeq.Data$ICR.score,method = "spearman"),2)
dev.new()
ggplot(RNASeq.Data, aes(x = ICR.score, y = HLA.A  )) +
  geom_point() +
  labs(title = paste0(GOF," in ",Cancerset," (",corelation,"[Spearman])"), x = "ICS score", y = "Log 2 Expression") +
  geom_smooth(method='lm',formula=y~x) +
  theme_classic()

#blots
GOF="HLA-DRA"
dev.new()
test.anova      = aov(HLA.DRA~ICR.cluster,data=RNASeq.Data)
p.value         = summary(test.anova)[[1]][["Pr(>F)"]][[1]]
if (p.value < 0.0001) {
  p.value.rounded = "p < 0.0001"
} else {
  p.value.rounded = paste0("p = ",round(p.value,4))
}
ggplot(RNASeq.Data, aes(x = ICR.cluster, y = HLA.DRA  )) +
  geom_boxplot() +
  labs(title = paste0(GOF," in ",Cancerset," (",p.value.rounded,"[ANOVA])"), x = "Cluster", y = "Log 2 Expression") +
  theme_classic()

corelation = round(cor(x = RNASeq.Data$HLA.DRA,y=RNASeq.Data$ICR.score,method = "spearman"),2)
dev.new()
ggplot(RNASeq.Data, aes(x = ICR.score, y = HLA.DRA  )) +
  geom_point() +
  labs(title = paste0(GOF," in ",Cancerset," (",corelation,"[Spearman])"), x = "ICS score", y = "Log 2 Expression") +
  geom_smooth(method='lm',formula=y~x) +
  theme_classic()

dev.new()
heatmap(x = as.matrix(t(RNASeq.Data[,c(1,2,3,5)])),
        scale = "row",
        Rowv = NA,
        labCol= FALSE)
title(main=paste0("HLA_heatmap_",Cancerset))

dev.new()
table_cluster_assignment = table_cluster_assignment[order(table_cluster_assignment$ICRscore),]
annotation = table_cluster_assignment[,c(1,6,8,10,13)]

levels(annotation$ICR_cluster_k3) = c("ICR Low", "ICR Medium1", "ICR High")
levels(annotation$ICR_cluster_k4) = c("ICR Low", "ICR Medium1", "ICR Medium2", "ICR High")
levels(annotation$ICR_cluster_k5) = c("ICR Low", "ICR Medium1", "ICR Medium2", "ICR Medium3", "ICR High")

annotation.blot = as.matrix(annotation[,-1])
annotation.blot[annotation.blot=="ICR Low"]="blue"
annotation.blot[annotation.blot=="ICR Medium1"]="green"
annotation.blot[annotation.blot=="ICR Medium2"]="yellow"
annotation.blot[annotation.blot=="ICR Medium3"]="orange"
annotation.blot[annotation.blot=="ICR High"]="red"
annotation.blot[annotation.blot=="ICR Medium"]="green"

my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
expression.matrix = t(ICR.heatmap)[,rownames(annotation.blot)]
heatmap.3(x = expression.matrix,
        scale = "row",
        col=my.palette,
        Rowv = NA,
        Colv = NA,
        labCol= FALSE)
title(main=paste0("HLA_heatmap_",Cancerset))
expression.matrix = t(HLA.heatmap)[,rownames(annotation.blot)]
heatmap.3(x = expression.matrix,
        scale =  "row",
        col=my.palette,
        Rowv = NA,
        Colv = NA,
        labCol= FALSE)
#title(main=paste0("HLA_heatmap_",Cancerset))

dev.new()
HLA.heatmap$ICR.score <- table_cluster_assignment$ICRscore[match(rownames(RNASeq.Data),rownames(table_cluster_assignment))]
corelation.matrix <- cor(HLA.heatmap) 
corrplot.mixed (corelation.matrix,
                #type="lower",
                #p.mat = Candidate_cor_sign[[1]],                                                                            # add significance to correlations
                #col =  colpattern,
                lower = "square",
                upper ="number",
                order="FPC",
                #cl.lim=lims,                                                                                               # only positive correlations
                tl.pos ="lt",
                #tl.col = as.character(annotation$color),
                insig= "pch",                                                                                              # remove insignificant correlations
                pch = "x",
                pch.cex= 1.5,
                tl.cex = 1,
                cl.cex = 1/par("cex"),
                cex.main = 1/par("cex"),
                #mar=c(6,4.1,7,5)
                )
