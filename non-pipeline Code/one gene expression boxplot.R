# Setup environment
setwd("~/Dropbox (TBI-Lab)/BREAST_QATAR/")
rm(list=ls())
## dependencies
library (ggplot2)

#parameters
exp.gene = "SLFN11"
mut.gene = "TP53"
Cancerset = "BRCA.PCF"
Parent.Cancerset = Cancerset
Parent.Cancerset = substring(Cancerset,1,4)
Geneset = "DBGS3.FLTR"

# Load data files
load (paste0("./2 DATA/TCGA RNAseq/RNASeq_",Parent.Cancerset,"_EDASeq/",Parent.Cancerset,".RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"),header=TRUE)
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",Geneset,".split.RDATA"))
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",
                                   Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
rownames(Consensus.class) <- Consensus.class$PatientID

#expression of gene
Exp.data <- as.data.frame(t(RNASeq.NORM_Log2[exp.gene,,drop=FALSE]))
rm(RNASeq.NORM_Log2)

#mutation status for gene
MUT.samples <- as.character(unique (Mutation.All[Mutation.All$Hugo_Symbol==mut.gene,"Patient_ID"]))

#gather plotting data
Plot.data <- Exp.data
if (Parent.Cancerset == "BRCA") {Plot.data$subtype <- Clinical.data$TCGA.PAM50.RMethod.RNASeq[match(rownames(Plot.data),rownames(Clinical.data))]
                                Plot.data$stage <- Clinical.data$ajcc_pathologic_tumor_stage[match(rownames(Plot.data),rownames(Clinical.data))]}
if (Cancerset == "OV") {Plot.data$stage <- Clinical.data$clinical_stage[match(rownames(Plot.data),rownames(Clinical.data))]}
if (Cancerset == "SKCM") {Plot.data$stage <- Clinical.data$ajcc_pathologic_tumor_stage[match(rownames(Plot.data),rownames(Clinical.data))]}
if (Parent.Cancerset == "COAD") {Plot.data$stage <- Clinical.data$ajcc_pathologic_tumor_stage[match(rownames(Plot.data),rownames(Clinical.data))]}
Plot.data[Plot.data$stage=="[Not Available]" & !is.na(Plot.data$stage),][,"stage"] <- NA
Plot.data$MUTSTAT <- "WT"
Plot.data[rownames(Plot.data)%in% MUT.samples,"MUTSTAT"] <- "MUT"
Plot.data$cluster <- Consensus.class$Group[match(rownames(Plot.data),rownames(Consensus.class))]

#filter IMS
#Plot.data<-Plot.data[Plot.data$subtype=="HER2-enriched",]

#plot
if (Parent.Cancerset == "BRCA") {
dev.new()
subtype.order = c("HER2-enriched", "Basal-like" , "Luminal B" ,"Luminal A","Normal-like")
colors = c("#da70d6","#daa520","#eaff00","#00c0ff","#d3d3d3")
gg = ggplot(Plot.data, aes_string("subtype", exp.gene, fill="subtype")) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  scale_x_discrete(limits=subtype.order) +
  scale_fill_manual(values = colors) +
  theme_bw()
print (gg)
}

dev.new()
#stage.order = c("HER2-enriched", "Basal-like" , "Luminal B" ,"Luminal A","Normal-like")
#colors = c("#da70d6","#daa520","#eaff00","#00c0ff","#d3d3d3")
gg = ggplot(Plot.data, aes_string("stage", exp.gene , fill="stage")) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  #scale_x_discrete(limits=subtype.order) +
  #scale_fill_manual(values = colors) +
  theme_bw()
print (gg)

dev.new()
#stage.order = c("HER2-enriched", "Basal-like" , "Luminal B" ,"Luminal A","Normal-like")
#colors = c("#da70d6","#daa520","#eaff00","#00c0ff","#d3d3d3")
gg = ggplot(Plot.data, aes_string("MUTSTAT", exp.gene , fill="MUTSTAT")) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  #scale_x_discrete(limits=subtype.order) +
  #scale_fill_manual(values = colors) +
  theme_bw()
print (gg)

dev.new()
cluster.order = c("ICR1", "ICR2" , "ICR3" ,"ICR4")
colors = c("blue","green","orange","red")
gg = ggplot(Plot.data, aes_string("cluster", exp.gene , fill="cluster")) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  scale_x_discrete(limits=cluster.order) +
  scale_fill_manual(values = colors) +
  theme_bw()
print (gg)
