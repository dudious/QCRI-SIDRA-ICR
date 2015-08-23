#########################
## Script to perform pecific gene barplot
## Input: Mutation .Frequencies.RDATA file, and the cluster assignment file (sample name, cluster assignment)
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

## Parameters
Cancerset <- "BRCA"           # do not use -GA or -hiseq (data is merged)
BRCA.Filter <- "BSF"          # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
GOF = "TP53"

## Load Mutaion data
load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
Mutation.selected.data <- data.frame(Hugo_Symbol = maf.merged.table$Hugo_Symbol, Variant_Classification = maf.merged.table$Variant_Classification, Patient_ID = substr(maf.merged.table$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins

#clinical data
Cancerset = paste0(Cancerset,".",BRCA.Filter)
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
#add subtype to Mutation.Frequency.Patient Table
Mutation.selected.data$Subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(Mutation.selected.data$Patient_ID,rownames(ClinicalData.subset))]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Subtype)),]

#cluster assignment
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]
# Add Class to mutation data
Mutation.selected.data$Cluster <- Consensus.class$Cluster [match(Mutation.selected.data$Patient_ID,Consensus.class$Patient_ID)]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Cluster)),]


#collapse IMS and mutation type
# Luminal A + Luminal B = Luminal
# Frame_Shift_Ins + Frame_Shift_Del = Frame_Shift
levels (Mutation.selected.data$Variant_Classification) <- c(levels (Mutation.selected.data$Variant_Classification),"Frame_Shift")
Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Ins","Frame_Shift_Del"),"Variant_Classification"] <- "Frame_Shift"
Mutation.selected.data.luminal <- Mutation.selected.data
levels (Mutation.selected.data.luminal$Subtype) <- c(levels (Mutation.selected.data.luminal$Subtype),"Luminal")
Mutation.selected.data.luminal[Mutation.selected.data.luminal$Subtype %in% c("Luminal A","Luminal B"),"Subtype"] <- "Luminal"
Mutation.selected.data.luminal <- Mutation.selected.data.luminal[Mutation.selected.data.luminal$Subtype=="Luminal",]
Mutation.selected.data <- rbind(Mutation.selected.data,Mutation.selected.data.luminal)

#count patients / Cluster-IMS
Cluster_IMS.counts <- count(unique(Mutation.selected.data[,c("Patient_ID","Subtype","Cluster")]),vars = c("Subtype","Cluster"))
Cluster_IMS.counts$Group <- paste0(Cluster_IMS.counts$Cluster,".",Cluster_IMS.counts$Subtype)

# TEST : Number of IMS/cluster from Clindata.subset+Consensus.class vs mutationselected.data
test.df <- ClinicalData.subset[,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
test.df$Cluster <- Consensus.class$Cluster[match(rownames(test.df),Consensus.class$Patient_ID)]
test.df <- test.df[complete.cases(test.df),]
test.count <- count(test.df)
test.count$Group <- paste0(test.count$Cluster,".",test.count$TCGA.PAM50.RMethod.RNASeq)
test.count$Matched <- Cluster_IMS.counts$freq[match(test.count$Group,Cluster_IMS.counts$Group)]
#NOT EQUAL ???

#Filter by GOF
Mutation.selected.data.gof <- Mutation.selected.data [Mutation.selected.data$Hugo_Symbol == GOF,]

#count
#count.Any <- count (Mutation.selected.data.gof.Any[,-c(1,3)])
count.Separate <- count (Mutation.selected.data.gof,vars=c("Variant_Classification","Subtype","Cluster"))
count.NonSilent <- data.frame(Variant_Classification = "NonSilent", 
                              count (unique(Mutation.selected.data.gof[-which(Mutation.selected.data.gof$Variant_Classification=="Silent"),]),
                                     vars=c("Subtype","Cluster")))


#count.All.Missense <- count.All[count.All$Variant_Classification=="Missense_Mutation",] 
blot.df <- rbind(count.Separate,count.NonSilent)
colnames(blot.df) <- c("Mutation_Type","Molecular_Subtype","Cluster_Assignment","Mutation_Count")
#blot.df <- blot.df[-which (blot.df$Variant_Classification %in% c("In_Frame_Del","Silent")),]
#blot.df <- rbind(blot.df,c("Nonsense_Mutation","HER2-enriched","ICR1",0))
#blot.df <- rbind(blot.df,c("Nonsense_Mutation","Luminal","ICR1",0))
#blot.df <- rbind(blot.df,c("Splice_Site","Luminal","ICR4",0))
#blot.df <- rbind(blot.df,c("Frame_Shift","HER2-enriched","ICR1",0))
#blot.df <- rbind(blot.df,c("Frame_Shift","HER2-enriched","ICR3",0))
#blot.df <- blot.df[-27,]

class (blot.df$Mutation_Count) <-"numeric"
blot.df$Group <- paste0(blot.df$Cluster_Assignment,".",blot.df$Molecular_Subtype)
blot.df$Group_Count <- Cluster_IMS.counts$freq[match(blot.df$Group,Cluster_IMS.counts$Group)]
#blot.df$Group_Count <- test.count$freq[match(blot.df$Group,test.count$Group)] #alternative group count
class (blot.df$Group_Count) <-"numeric"
blot.df$Mutation_Frequency <- round((blot.df$Mutation_Count / blot.df$Group_Count) *100,1)



png(paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"By.IMS",".png", sep=""), height = 1000, width= 1500)
 gg = ggplot(blot.df, aes(x = Cluster_Assignment, y = Mutation_Frequency  )) +
              geom_bar(stat="identity",position="dodge",drop=FALSE,ylim=c(0,100)) +
              facet_grid(Molecular_Subtype~Mutation_Type, space="free")
 print(gg)
dev.off()

write.csv (blot.df,file=paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"By.IMS",".csv", sep=""))



