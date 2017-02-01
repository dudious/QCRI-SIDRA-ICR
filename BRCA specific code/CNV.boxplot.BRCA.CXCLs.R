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
#setwd("~/Dropbox/BREAST_QATAR/")
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")
# Dependencies
required.packages <- c("ggplot2")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "selected"    # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict", "chisqr"
IMS.filter     = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "All"         # Alternatives "1vs4" , "All"
gene.filter    = "FALSE"       # Alternatives "TRUE" , "FALSE"
genes.selection= c("CXCL9","CXCL10","CXCL11")


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
selected.data <- selected.data[order(factor(selected.data$cluster,levels = c("ICR4","ICR3","ICR2","ICR1")),
                                     factor(selected.data$subtype,levels = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like"))
                                     ),]
#append expression data
RNASeq.table <- as.data.frame(t(RNASeq.NORM_Log2[c(genes.selection),]))
selected.data$CXCL9.exp <- RNASeq.table$CXCL9[match(rownames(selected.data),rownames(RNASeq.table))]
selected.data$CXCL10.exp <- RNASeq.table$CXCL10[match(rownames(selected.data),rownames(RNASeq.table))]
selected.data$CXCL11.exp <- RNASeq.table$CXCL11[match(rownames(selected.data),rownames(RNASeq.table))]
#drop incomplete data (apply filter)
selected.data <- selected.data[complete.cases(selected.data),]

#expression vs CNV barchart
ICR.CXCL9.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","CXCL9.exp")],mean)
CNV.CXCL9.exp <- aggregate(.~CXCL9,data=selected.data[,c("CXCL9","CXCL9.exp")],mean)
ICR.CXCL10.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","CXCL10.exp")],mean)
CNV.CXCL10.exp <- aggregate(.~CXCL10,data=selected.data[,c("CXCL10","CXCL10.exp")],mean)
ICR.CXCL11.exp <- aggregate(.~cluster,data=selected.data[,c("cluster","CXCL11.exp")],mean)
CNV.CXCL11.exp <- aggregate(.~CXCL11,data=selected.data[,c("CXCL11","CXCL11.exp")],mean)

T1<-cbind("CXCL9.exp",ICR.CXCL9.exp)
T2<-cbind("CXCL10.exp",ICR.CXCL10.exp)
T3<-cbind("CXCL11.exp",ICR.CXCL11.exp)
T4<-cbind("CXCL9.exp",CNV.CXCL9.exp)
T5<-cbind("CXCL10.exp",CNV.CXCL10.exp)
T6<-cbind("CXCL11.exp",CNV.CXCL11.exp)

colnames(T1) = c("CXCL","Group","Expression")
colnames(T2) = c("CXCL","Group","Expression")
colnames(T3) = c("CXCL","Group","Expression")
colnames(T4) = c("CXCL","Group","Expression")
colnames(T5) = c("CXCL","Group","Expression")
colnames(T6) = c("CXCL","Group","Expression")

expression.table<-rbind (T1,T2,T3,T4,T5,T6)

gg<-ggplot(expression.table,aes(x=Group,y=Expression,fill=CXCL))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("CXCL expression") +
  ylab("Expression (Log2)") +
  theme_bw()
print(gg)
  
t.test(data=selected.data[selected.data$CXCL9 %in% c("Amplified","Normal"),],CXCL9.exp~CXCL9)
t.test(data=selected.data[selected.data$CXCL9 %in% c("Amplified","Deleted"),],CXCL9.exp~CXCL9)
t.test(data=selected.data[selected.data$CXCL9 %in% c("Normal","Deleted"),],CXCL9.exp~CXCL9)

t.test(data=selected.data[selected.data$CXCL10 %in% c("Amplified","Normal"),],CXCL10.exp~CXCL10)
t.test(data=selected.data[selected.data$CXCL10 %in% c("Amplified","Deleted"),],CXCL10.exp~CXCL10)
t.test(data=selected.data[selected.data$CXCL10 %in% c("Normal","Deleted"),],CXCL10.exp~CXCL10)

t.test(data=selected.data[selected.data$CXCL11 %in% c("Amplified","Normal"),],CXCL11.exp~CXCL11)
t.test(data=selected.data[selected.data$CXCL11 %in% c("Amplified","Deleted"),],CXCL11.exp~CXCL11)
t.test(data=selected.data[selected.data$CXCL11 %in% c("Normal","Deleted"),],CXCL11.exp~CXCL11)

dev.new()
cluster.order = c("Amplified", "Deleted" ,"Normal" )
colors = c("red", "blue","grey")
gg = ggplot(selected.data, aes(CXCL11, CXCL11.exp, fill=CXCL11)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(limits=cluster.order)
print (gg)

dev.new()
cluster.order = c("ICR4", "ICR3", "ICR2", "ICR1")
colors = rev(c("red", "orange","green", "blue"))
gg = ggplot(selected.data, aes(cluster, CXCL11.exp, fill=cluster)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(limits=cluster.order)
print (gg)
