#################
##
## MAP MUT vs WT
##
#################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
# Dependencies
required.packages <- c("gplots","plyr","Hmisc","data.table")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("gplots")
library("plyr")
library("Hmisc")
library("data.table")
source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "selected"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict", "chisqr"
IMS.filter     = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "All"         # Alternatives "1vs4" , "All"
gene.filter    = "FALSE"        # Alternatives "TRUE" , "FALSE"
selected.genes = c("TP53","MAP3K1","MAP2K4","CTCF","FCGBP","CXCL9")

# Load Data
# Mutation and amplification data
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
# clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
# cluster data
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#load and subset RNAseq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")
RNASeq.subset <- as.data.frame(t(RNASeq.NORM[,])) #select all / no subset used

#sELECT RELEVANT dATA
subset.merged.matrix <- as.data.frame(merged.matrix[,c("MAP3K1","MAP2K4")])
subset.merged.matrix$Subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(rownames(subset.merged.matrix),rownames(ClinicalData.subset))]
subset.merged.matrix$Cluster <- Consensus.class$Cluster[match(rownames(subset.merged.matrix),Consensus.class$Patient_ID)]  
subset.merged.matrix$MAP3K1.expression <- RNASeq.subset$MAP3K1[match(rownames(subset.merged.matrix),rownames(RNASeq.subset))]
subset.merged.matrix$MAP2K4.expression <- RNASeq.subset$MAP2K4[match(rownames(subset.merged.matrix),rownames(RNASeq.subset))]

#FILTER SUBTYPE
subset.merged.matrix.filtered <- subset.merged.matrix[subset.merged.matrix$Subtype == "Luminal B",]
  
#count
table (subset.merged.matrix.filtered$MAP3K1,useNA = "always")
table (subset.merged.matrix.filtered$MAP2K4,useNA = "always")


  
subset.merged.matrix.dt <- as.data.table(subset.merged.matrix.filtered)
subset.merged.matrix.dt$ID <- rownames(subset.merged.matrix.filtered)
setkey(subset.merged.matrix.dt,ID)
subset.merged.matrix.dt <- subset.merged.matrix.dt[-1,]


avg.MAP3K1.counts <- subset.merged.matrix.dt[, mean(MAP3K1.expression), by=MAP3K1]  
avg.MAP2K4.counts <- subset.merged.matrix.dt[, mean(MAP2K4.expression), by=MAP2K4]  

levels(subset.merged.matrix.dt$MAP3K1)<- c(levels(subset.merged.matrix.dt$MAP3K1),"NO_CNVorMUT")
levels(subset.merged.matrix.dt$MAP2K4)<- c(levels(subset.merged.matrix.dt$MAP2K4),"NO_CNVorMUT")
subset.merged.matrix.dt[is.na(subset.merged.matrix.dt)] <- "NO_CNVorMUT"


class(subset.merged.matrix.dt$MAP3K1.expression) = "numeric"
class(subset.merged.matrix.dt$MAP2K4.expression) = "numeric"

t.test(subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1=="AMP"], subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1=="HOMDEL"], subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1=="MUT"], subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1=="MUT;AMP"], subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1=="MUT;HOMDEL"], subset.merged.matrix.dt$MAP3K1.expression[subset.merged.matrix.dt$MAP3K1 == "NO_CNVorMUT"])

t.test(subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4=="AMP"], subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4=="HOMDEL"], subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4=="MUT"], subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4=="MUT;AMP"], subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4 == "NO_CNVorMUT"])
t.test(subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4=="MUT;HOMDEL"], subset.merged.matrix.dt$MAP2K4.expression[subset.merged.matrix.dt$MAP2K4 == "NO_CNVorMUT"])




