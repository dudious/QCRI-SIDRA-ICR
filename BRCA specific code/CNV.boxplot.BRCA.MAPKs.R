#################################################################
###
### This Script PLots boxplot for MAPKs based on 
### Consensus Clustering clustering of RNASeq Data and mutation data
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

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "selected"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict", "chisqr"
IMS.filter     = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "All"         # Alternatives "1vs4" , "All"
gene.filter    = "FALSE"        # Alternatives "TRUE" , "FALSE"
genes.selection= c("MAP3K1","MAP2K4")
#selected.genes = c("TP53","MAP3K1","MAP2K4","CTCF","FCGBP","CXCL9")

# Load Data
#combined mutation and CNV datamatix
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
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
#mutation variation data
#load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".",matrix.type,".VariationTables.RData"))
mutstat.data <- read.csv("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
mutstat.data$X  <-NULL
#DEG stats
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsall.rdata")
DEGsMAPKs.AllIMS<-DEGs
rm(DEGs)
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.rdata")
DEGsMAPKs.Luminal<-DEGs
rm(DEGs)
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGinICR4vs1.RDATA")
#load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/CNV/DEGin.MAPKsDELvsMAPKsAMP.RDATA")
#load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/CNV/DEGin.MAPKsDELvsMAPKsAMP.Luminal.RDATA")

# cluster selection
if (cluster.select == "1vs4") {
  Consensus.class <- Consensus.class[Consensus.class$Cluster %in% c("ICR1","ICR4"),]
}
# subtype selection
if (IMS.filter == "Luminal") {
  ClinicalData.subset <- ClinicalData.subset[ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq %in% c("Luminal A","Luminal B"),]
}

# select data to plot
selected.data <- merged.matrix[,genes.selection]
selected.data[selected.data=="MUT"] <- NA
selected.data <- gsub("MUT;","",as.matrix(selected.data))
selected.data[selected.data=="HOMDEL"] <- "Deleted"
selected.data[selected.data=="AMP"] <- "Amplified"
selected.data <- as.matrix(selected.data)
selected.data[is.na(selected.data)] <- "Normal"

#append meta data
selected.data <- as.data.frame(selected.data)
selected.data$subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(rownames(selected.data),rownames(ClinicalData.subset))]
selected.data$cluster <- Consensus.class$Cluster[match(rownames(selected.data),rownames(Consensus.class))]
selected.data$MUTSTAT <- mutstat.data$MAP2K4.MAP3K1[match(rownames(selected.data),mutstat.data$Sample)]
selected.data <- selected.data[order(factor(selected.data$cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                     factor(selected.data$subtype,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"))
                                     ),]
#append expression data
RNASeq.table <- as.data.frame(t(RNASeq.NORM_Log2[c(genes.selection),]))
selected.data$MAP3K1.exp <- RNASeq.table$MAP3K1[match(rownames(selected.data),rownames(RNASeq.table))]
selected.data$MAP2K4.exp <- RNASeq.table$MAP2K4[match(rownames(selected.data),rownames(RNASeq.table))]

#drop incomplete data (apply filter)
selected.data <- selected.data[complete.cases(selected.data),]

#expression vs CNV barchart
ICR.MAP3K1.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","MAP3K1.exp")],mean)
ICR.MAP2K4.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","MAP2K4.exp")],mean)
MUTSTAT.MAP3K1.exp <- aggregate(.~MUTSTAT,data=selected.data[,c("MUTSTAT","MAP3K1.exp")],mean)
MUTSTAT.MAP2K4.exp <- aggregate(.~MUTSTAT,data=selected.data[,c("MUTSTAT","MAP2K4.exp")],mean)
CNV.MAP3K1.exp <- aggregate(.~MAP3K1,data=selected.data[,c("MAP3K1","MAP3K1.exp")],mean)
CNV.MAP2K4.exp <- aggregate(.~MAP2K4,data=selected.data[,c("MAP2K4","MAP2K4.exp")],mean)

T1<-cbind("MAP3K1.exp",ICR.MAP3K1.exp)
T2<-cbind("MAP2K4.exp",ICR.MAP2K4.exp)
T3<-cbind("MAP3K1.exp",MUTSTAT.MAP3K1.exp)
T4<-cbind("MAP2K4.exp",MUTSTAT.MAP2K4.exp)
T5<-cbind("MAP3K1.exp",CNV.MAP3K1.exp)
T6<-cbind("MAP2K4.exp",CNV.MAP2K4.exp)

colnames(T1) = c("MAPK","Group","Expression")
colnames(T2) = c("MAPK","Group","Expression")
colnames(T3) = c("MAPK","Group","Expression")
colnames(T4) = c("MAPK","Group","Expression")
colnames(T5) = c("MAPK","Group","Expression")
colnames(T6) = c("MAPK","Group","Expression")

expression.table<-rbind (T1,T2,T3,T4,T5,T6)

gg<-ggplot(expression.table,aes(x=Group,y=Expression,fill=MAPK))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("MAPK expression") +
  ylab("Expression (Log2)") +
  theme_bw()
print(gg)
  