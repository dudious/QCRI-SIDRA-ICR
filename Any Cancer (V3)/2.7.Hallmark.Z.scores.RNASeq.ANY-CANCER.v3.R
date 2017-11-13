#################################################################
###
### This Script generates a zscore for the expression of each  
### cancer hallmark patway
###
#################################################################

# Setup environment
rm(list=ls())
#setwd("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/")
#setwd("/mnt3/wouter/BREAST-QATAR/")
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")
#Dependencies
required.packages <- c("corrplot")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (corrplot)
library (plyr)

# Set Parameters
DL.Method    = "BIOLINKS" #Choose "ASSEMBLER" or "BIOLINKS"
sample.types = "Selected" #Alternatives TP , TP_TM , Selected
Cancersets <- "ALL"           # do not use -GA or -hiseq (data is merged)
#BRCA.Filter <- "PCF"         # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3"
Genedatabase <- "Gene_selection_v2.7.txt"

# Load immune genes
gene.list <- read.csv (paste0("./2 DATA/SUBSETS/",Genedatabase))                                 # Select subset here !!!!! and change filename below !!!!
gene.list.selected <- as.character(gene.list[which(gene.list[,Geneset]==1),1])
# load hallmark pathways data
load("./2 DATA/Hallmark Cancer Pathways/cancer.halmark.pathways.R")
hallmark.selected <- read.csv("./2 DATA/Hallmark Cancer Pathways/cancer.hallmark.pathways.selection.csv",header = FALSE,stringsAsFactors = FALSE)[,c(1,2)]
colnames(hallmark.selected) <- c("name","selected")
# selected pathways only
hallmark.selected.pathways <- hallmark.selected[hallmark.selected$selected == "1","name"]
hallmark.pathways <- hallmark.pathways[hallmark.selected.pathways]
# add selected IPA pathways
IPA.pathways <- read.csv("./2 DATA/pathways_from_IPA_selection.csv",header = TRUE,stringsAsFactors = FALSE)
IPA.pathways.selected <- IPA.pathways[IPA.pathways$Selected=="1",]
IPA.pathways.selected$Ingenuity.Canonical.Pathways <- toupper(gsub(" ","_",IPA.pathways.selected$Ingenuity.Canonical.Pathways))
IPA.pathways.selected$Selected <- NULL
IPA.pathways.selected.list <- strsplit(IPA.pathways.selected$Molecules,split=",")
names(IPA.pathways.selected.list) <- IPA.pathways.selected$Ingenuity.Canonical.Pathways
# merge pathways List
Selected.pathways <- c(hallmark.pathways,IPA.pathways.selected.list)
# add barier genes
Selected.pathways$BARRIER_GENES <- c("FLG","TACSTD2","DSC3","DST","DSP","PPL","PKP3","JUP")
# add MAPK ICR1vs4 luminal mutation signature
Selected.pathways$MAPK_UP_GENES <- c("TAOK2","TP53","MAPK3","MAP3K1","MAPT","HSPA1A","FLNB","TAOK3","CRK","RPS6KA2","MAP2K4","DUSP5","CACNA1D","MAPK8","RASGRP1","CACNA1G")
Selected.pathways$MAPK_DOWN_GENES <- c("CACNG6","CACNA1B","CACNA2D3","FASLG","RASGRF1","JUN","JUND","DUSP16","PPM1B","SOS1","FGF12","RASGRP2","PRKCB","MAP4K1","PTPN7","GADD45G","DDIT3","DUSP8","DUSP10","FGFR4","FGF14","FGF13","MAP2K6","DUSP2")

#Save pathway file
save(Selected.pathways,file = "./PANCANCER/Selected.pathways.RData")
Selected.pathways.df <- data.frame(Name = names(Selected.pathways),
                                   Number = sapply(Selected.pathways, length),
                                   Genes = sapply(names(Selected.pathways),function(x) paste(Selected.pathways[[x]],collapse=" ")),
                                   row.names = NULL)
write.csv(Selected.pathways.df,"./PANCANCER/Selected.pathways.txt",row.names = F)


# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}
  Parent.Cancerset <- substring(Cancerset,1,4)

# load RNASeq data
if (DL.Method == "BIOLINKS") {
  load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".",sample.types,".NORMALIZED.TP_FILTERED_LOG2.RData"))
  RNASeq.NORM_Log2<-RNASeq.NORM.TP_Log2
  rm(RNASeq.NORM.TP_Log2)
  }
if (DL.Method == "ASSEMBLER") {load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_EDASeq/",Cancerset,".RNASeq.TCGA.",DL.Method,".NORMALIZED_LOG2.RData"))}
print (paste0(Cancerset," RNASeq data Loaded..."))
# check availabilety of the Immune genes in the dataset
available.genes.RNAseq <- gene.list.selected[which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]
unavailable.genes.RNAseq <- gene.list.selected[-which(gene.list.selected %in% rownames(RNASeq.NORM_Log2))]
Selected.pathways$IMMUNE_GENES_DBGS3 <- available.genes.RNAseq

#z-score matrix
Selected.pathways.Zscores <- data.frame(matrix(ncol=length(Selected.pathways), nrow=ncol(RNASeq.NORM_Log2)))
rownames(Selected.pathways.Zscores) = colnames(RNASeq.NORM_Log2)
colnames(Selected.pathways.Zscores) = names(Selected.pathways)

for(i in 1:length(Selected.pathways)){
  pathway.genes <- unlist(Selected.pathways[i],use.names = FALSE)
  present.genes <- pathway.genes[which(pathway.genes %in% rownames(RNASeq.NORM_Log2))]
  assign(paste(names(Selected.pathways[i]),".Expresion.matrix", sep=""),RNASeq.NORM_Log2[present.genes,])
  Selected.pathways.Zscores[,names(Selected.pathways[i])] <- colMeans(RNASeq.NORM_Log2[present.genes,],na.rm = TRUE)
}
print (paste0(Cancerset," Pathways Z-scores calculated..."))

# Correlation plot
Selected.pathways.spearcor <- cor (Selected.pathways.Zscores,method="spearman")

# cor significance
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
Selected.pathways.spearcor.sign <- cor.mtest(Selected.pathways.spearcor, 0.95)

png(paste0("./4 FIGURES/CORRELATION/Selected_Pathways/",DL.Method,".",sample.types,".correlation.",Cancerset,".",Geneset,".png"),res=600,height=6,width=6,unit="in")  #adjust output file names here !!!!!
cex.before <- par("cex")
par(cex = 0.35)
col1 = colorRampPalette(c("blue", "white", "#009900"))
#lims=c(-1,1)
#if (length(RNASeq.subset.cor[RNASeq.subset.cor<0]) == 0) {lims=c(0,1)} 
corrplot.mixed (Selected.pathways.spearcor,
                #type="lower",
                p.mat = Selected.pathways.spearcor.sign[[1]],    # add significance to correlations
                col = col1(100),
                lower = "square",
                upper = "square",
                order="FPC",
                #cl.lim=lims,                      # only positive correlations
                tl.pos ="lt",
                tl.col = "#c00000",
                insig="blank",                          # remove insignificant correlations
                tl.cex = 1,
                cl.cex = 1/par("cex"),
                title = paste0("Sign.Spearman TCGA-RNASeq (",Cancerset,".",Geneset,"vs SelectedPW)"),
                cex.main = 1/par("cex")
                #mar=c(5.1,4.1,4.1,2.1)
                )
#par(cex = cex.before)
dev.off()

#Save immune data
p.value.matrix <- as.matrix(Selected.pathways.spearcor.sign[[1]])
colnames(p.value.matrix)=colnames(Selected.pathways.spearcor)
rownames(p.value.matrix)=rownames(Selected.pathways.spearcor)
immune.data <- data.frame(correlation = Selected.pathways.spearcor["IMMUNE_GENES_DBGS3",],
                           p.value = p.value.matrix["IMMUNE_GENES_DBGS3",])
save(immune.data,file=paste0("./3 ANALISYS/Pathway correlation/",DL.Method,".",sample.types,"Immune.correlation.",Cancerset,".",Geneset,".RData"))
print (paste0(Cancerset," Correlation data saved."))
}
