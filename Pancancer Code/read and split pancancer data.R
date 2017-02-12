#################################################################
###
### This Script ....
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")

#Dependencies
required.packages <- c("stringr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("stringr")

#parameters
Genedatabase = "Gene_selection_v2.6.txt"

#load data
clinical.data <- read.csv("./2 DATA/Clinical Information/PANCANCER/clinical_PANCAN_patient_with_followup.tsv",sep = "")
#RNAseq.data <- read.csv("./Expression/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "")
RNAseq.annot <- read.csv("./2 DATA/TCGA RNAseq/RNASeq_PANCANCER/EB++GeneExpAnnotation.tsv",sep = "")
load ("./2 DATA/TCGA RNAseq/RNASeq_PANCANCER/E.RData")
RNAseq.data <- Exp
rm(Exp)
rownames (RNAseq.data) <- RNAseq.data$gene_id
RNAseq.data$gene_id <-NULL
RNAseq.data <- as.matrix(RNAseq.data)
mode(RNAseq.data)
RNASeq.NORM_Log2 <- apply(RNAseq.data,1,function(x) log(x+1,2))
#save (RNASeq.NORM_Log2,RNAseq.annot,file="./2 DATA/TCGA RNAseq/RNASeq_PANCANCER/Expr.Log2.RData")
rm(RNAseq.data)

load ("./2 DATA/TCGA RNAseq/RNASeq_PANCANCER/Expr.Log2.RData")
rownames(RNASeq.NORM_Log2) <- gsub("\\.","-",rownames(RNASeq.NORM_Log2))
RNAseq.gene.annot <- as.data.frame(str_split_fixed(colnames(RNASeq.NORM_Log2),"\\|",2),stringsAsFactors=FALSE)
colnames(RNAseq.gene.annot) <- c("name","ID")
RNASeq.NORM_Log2_Filtered <- RNASeq.NORM_Log2[,-which(RNAseq.gene.annot$name=="?")]
RNAseq.gene.annot <- RNAseq.gene.annot[-which(RNAseq.gene.annot$name=="?"),]
RNAseq.gene.annot[RNAseq.gene.annot$ID=="728661","name"] <- "SLC35E2B"
colnames(RNASeq.NORM_Log2_Filtered) <- RNAseq.gene.annot$name
#rownames(RNASeq.NORM_Log2_Filtered) <- gsub("\\.","-",substring(rownames(RNASeq.NORM_Log2_Filtered),1,12))
save (RNASeq.NORM_Log2_Filtered,RNAseq.annot,file="./2 DATA/TCGA RNAseq/RNASeq_PANCANCER/Expr.Log2.Filtered.RData")

#Subset the Pancancer Matrix
gene.list <- read.csv (paste0("./2 DATA/SUBSETS/",Genedatabase))                                 # Select subset here !!!!! and change filename below !!!!
gene.list.ALL <- as.character(gene.list[which(gene.list[,"DBGS3"]==1),1])
gene.list.INH <- as.character(gene.list[which(gene.list[,"ImSuGS"]==1),1])
a<-which(gene.list[,"DBGS3"]==1)
b<-which(gene.list[,"ImSuGS"]==1)
gene.list.ACT <- as.character(gene.list[a[-which(a%in%b)],1])
RNAseq.ALL <- t(RNASeq.NORM_Log2_Filtered[,gene.list.ALL])
RNAseq.INH <- t(RNASeq.NORM_Log2_Filtered[,gene.list.INH])
RNAseq.ACT <- t(RNASeq.NORM_Log2_Filtered[,gene.list.ACT])
save (RNAseq.ALL,RNAseq.INH,RNAseq.ACT,file="./2 DATA/SUBSETS/PANCANCER/TCGA.PANCANCER.RNASeq.subset.DBGS3.RData")
write.csv()

#double patients ? 795
length(substring(colnames(RNAseq.data),1,12)) - length(unique(substring(colnames(RNAseq.data),1,12)))
#sample to patient barcode
sample.meta.data <- data.frame(sample.barcode=gsub("\\.","-",as.character(rownames(RNASeq.NORM_Log2))),
                               patient.barcode=gsub("\\.","-",substring(rownames(RNASeq.NORM_Log2),1,12)),
                               stringsAsFactors = FALSE)



# Calculate pancancer immune score
load ("./2 DATA/SUBSETS/PANCANCER/TCGA.PANCANCER.RNASeq.subset.DBGS3.RData")
IS.DBGS3.PANCANCER.unsplit <- data.frame (Patient_ID = colnames(RNAseq.ALL),
                                          ALL = colMeans(RNAseq.ALL),
                                          INH = colMeans(RNAseq.INH),
                                          ACT = colMeans(RNAseq.ACT))
write.csv (IS.DBGS3.PANCANCER.unsplit,"./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.UNSPLIT.DBGS3.csv")

#QC IS
IS.BLCA <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.BLCA.DBGS3.csv",stringsAsFactors = FALSE)
which(duplicated (IS.DBGS3.PANCANCER.unsplit$Patient_ID))
IS.DBGS3.PANCANCER.unsplit.QC <- IS.DBGS3.PANCANCER.unsplit[IS.DBGS3.PANCANCER.unsplit$Patient_ID %in% IS.BLCA$X,]
IS.DBGS3.PANCANCER.unsplit.QC <- IS.DBGS3.PANCANCER.unsplit.QC[-which(duplicated (IS.DBGS3.PANCANCER.unsplit.QC$Patient_ID)),]
IS.BLCA <- IS.BLCA[-which(duplicated (IS.BLCA$X)),]
row.names(IS.BLCA) <- IS.BLCA$X
IS.BLCA$X <- NULL
rownames(IS.DBGS3.PANCANCER.unsplit.QC) <- IS.DBGS3.PANCANCER.unsplit.QC$Patient_ID
IS.BLCA$QC <- IS.DBGS3.PANCANCER.unsplit.QC$ALL[match(rownames(IS.BLCA),rownames(IS.DBGS3.PANCANCER.unsplit.QC))]
load ("./2 DATA/SUBSETS/PANCANCER/BLCA/TCGA.BLCA.RNASeq.Split.subset.DBGS3.RData")
RNASeq.subset["TCGA-2F-A9KO",]
RNAseq.ALL[,"TCGA-2F-A9KO"]
mean(RNAseq.ALL[,"TCGA-2F-A9KO"])
mean(RNASeq.subset["TCGA-2F-A9KO",])
IS.BLCA["TCGA-2F-A9KO","unscaled.IS"]
IS.DBGS3.PANCANCER.unsplit.QC["TCGA-2F-A9KO","ALL"]

#correlation
IS.PANCANCER <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.UNSPLIT.DBGS3.csv")
cor(IS.PANCANCER$ALL,IS.PANCANCER$INH,method = "pearson")
cor(IS.PANCANCER$ALL,IS.PANCANCER$ACT,method = "pearson")
cor(IS.PANCANCER$INH,IS.PANCANCER$ACT,method = "pearson")

#test split
i = 2

## split by cancertype
for ( i in 1:length(levels (clinical.data$acronym))){
  cancer.type <- (levels (clinical.data$acronym)[i])
  patient.selection <- as.character(clinical.data[clinical.data$acronym==cancer.type,"bcr_patient_barcode"])
  sample.selection <- sample.meta.data[sample.meta.data$patient.barcode %in% patient.selection,]
  RNAseq.matrix <- t(RNASeq.NORM_Log2[sample.selection$sample.barcode,])
  RNAseq.gene.annot <- as.data.frame(str_split_fixed(rownames(RNAseq.matrix),"\\|",2),stringsAsFactors=FALSE)
  colnames(RNAseq.gene.annot) <- c("name","ID")
  RNAseq.matrix <- RNAseq.matrix[-which(RNAseq.gene.annot$name=="?"),]
  RNAseq.gene.annot <- RNAseq.gene.annot[-which(RNAseq.gene.annot$name=="?"),]
  RNAseq.gene.annot[RNAseq.gene.annot$ID=="728661","name"] <- "SLC35E2B"
  rownames(RNAseq.matrix) <- RNAseq.gene.annot$name
  #colnames(RNAseq.matrix) <- gsub("\\.","-",substring(colnames(RNAseq.matrix),1,12))
  #colnames(RNAseq.matrix) <- gsub("\\.","-",colnames(RNAseq.matrix))
  RNAseq.matrix <- t(RNAseq.matrix)
  mode(RNAseq.matrix) <- "numeric"
  clinical.matched <- clinical.data[clinical.data$bcr_patient_barcode %in% substring(rownames(RNAseq.matrix),1,12),]
  rownames(clinical.matched) <- clinical.matched$bcr_patient_barcode
  save (RNAseq.matrix,RNAseq.gene.annot,clinical.matched,file = paste0("./2 DATA/TCGA RNAseq/RNASeq_PANCANCER/",cancer.type,".RNASeq.TCGA.PANCANCER.SPLIT.DATA.RData"))
}

## set wd to my folder 
setwd("/mnt3/wouter/pancancer")

## hallmark pathways data
#hallmark.pathways <- scan("./h.all.v5.1.symbols.gmt", what="", sep="\n")
#hallmark.pathways <- strsplit(hallmark.pathways, "[[:space:]]+")
#names(hallmark.pathways) <- sapply(hallmark.pathways, `[[`, 1)
#hallmark.pathways <- lapply(hallmark.pathways, `[`, -c(1,2))
#save (hallmark.pathways,file = "./cancer.halmark.pathways.R")
