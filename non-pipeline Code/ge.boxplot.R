#################################################################
###
### This Script Plots boxplot for genes based on 
### Consensus Clustering clustering of RNASeq Data 
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/RNAseq/...
### Data is saved :
### NO DATA
### Figures are saved :
### ./4 FIGURES/CNV
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
# Dependencies
required.packages <- c("ggplot2")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")
require(reshape2)

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "selected"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict", "chisqr"
IMS.filter     = "Luminal"         # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "All"         # Alternatives "1vs4" , "All"
gene.filter    = "FALSE"        # Alternatives "TRUE" , "FALSE"
genes.selection= c("MAP3K1","MAP2K4")

### Read other genes for plotting
genes.list = read.table("zscore-genes.txt", header = F, stringsAsFactors = F)
genes.list = genes.list$V1
genes.selection = genes.list

# Load Data
#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
#clustering Data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]
#RNAseq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")



RNASeq.table <- as.data.frame(t(RNASeq.NORM_Log2[c(genes.selection),]))
exp.table = RNASeq.table
exp.table$cluster = Consensus.class$Cluster[match(rownames(RNASeq.table),rownames(Consensus.class))]
exp.table <- exp.table[complete.cases(exp.table),]

############################
## Boxplot for genes expression
###########################

cluster.order = paste0("ICR", c(4:1))
colors = c("red", "orange", "green", "blue")
data.table = melt(exp.table)
data.table$cluster <- factor(data.table$cluster,levels = cluster.order,ordered = TRUE)

gg<-ggplot(data.table , aes(x=cluster, y=value, fill = cluster))+
  geom_boxplot(notch = T) + 
  stat_boxplot(geom ='errorbar') +
  geom_point(position=position_dodge(width=0.75, height = 0.5) )+
  facet_grid(. ~ variable, scales="free_x", space="free_x") +
  scale_fill_manual(values = colors) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5))+
  ggtitle("Up-regulated genes")
print(gg)
###################