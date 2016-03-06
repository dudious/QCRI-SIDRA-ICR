# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
## dependencies
required.packages <- c("ggplot2")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library(ggplot2)

#load data
load("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.All.DBGS3.FLTR.Mutation.Matrixes.NonSilent.oncoplot.Rdata")
load("./3 ANALISYS/Mutations/BRCA.BSF2/Mutation.Data.TCGA.BRCA.BSF2.All.DBGS3.FLTR.Frequencies.RDATA")

#mutation state of selected gene or geneset
selected.genes = c("FCGBP")
selected.gene.mutations <- genes.mutations[,selected.genes,drop=FALSE]
selected.gene.mutations$DUMMY <- NA
WT.samples <- rownames(selected.gene.mutations[which(apply(selected.gene.mutations, 1, function(x) all(is.na(x)))),])
MUT.samples <- rownames(selected.gene.mutations[-which(apply(selected.gene.mutations, 1, function(x) all(is.na(x)))),])
MUT.state <- rbind(data.frame(Sample_ID=WT.samples,MUT_State="WT"),data.frame(Sample_ID=MUT.samples,MUT_State="MUT"))
table (MUT.state$MUT_State)
#Add mutation load
MUT.state$MUT_Load <- Mutation.Frequency.Patient$Freq.NonSilent[match(MUT.state$Sample_ID,Mutation.Frequency.Patient$Patient_ID)]
#Blot
dev.new(width=6, height=6)
gg = ggplot(MUT.state,aes(x=MUT_State,y=log(MUT_Load)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.1,height=0.1)) +
  theme_bw() +
  ggtitle(paste(selected.genes,sep=" & ",collapse=" & "))
print(gg)

all(is.na(selected.gene.mutations))

