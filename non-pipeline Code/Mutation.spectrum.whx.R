#########################
## Script to perform mutation spectrum stacked barplot
## Input: Mutation .maf file, and the cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
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
Cancerset  = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF ,Dont use -GA or -hiseq
Geneset    = "DBGS3.FLTR"  # SET GENESET HERE !!!!!!!!!!!!!!
IMS.filter = "All"   # Alterantives "All" , "Luminal" , "Basal", "Her2" "Luminal A" "Luminal B"
stats      = ""            # Alterantives : "stats" ""
GOF        = "FALSE"

## Load Data
# Load MAF file
load ("./2 DATA/TCGA Mutations/BRCA/Somatic_Mutations/BRCA.TCGA.combined.Mutation.Data.maf.Rdata")
mutated.allele.data <- maf.merged.table[,c("Hugo_Symbol","Tumor_Sample_Barcode","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")]
mutated.allele.data$Patient_ID <- substr(mutated.allele.data$Tumor_Sample_Barcode,1,12)
#cluster assignment
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]
#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL


## check allele mutation 
mutated.allele.data$allele1.pair <- paste0(mutated.allele.data$Match_Norm_Seq_Allele1,">",mutated.allele.data$Tumor_Seq_Allele1)
mutated.allele.data$allele2.pair <- paste0(mutated.allele.data$Match_Norm_Seq_Allele2,">",mutated.allele.data$Tumor_Seq_Allele2)
mutated.allele.data[mutated.allele.data$allele1.pair == paste0(mutated.allele.data$Match_Norm_Seq_Allele1,">",mutated.allele.data$Match_Norm_Seq_Allele1),"allele1.pair"] <- ""
mutated.allele.data[mutated.allele.data$allele2.pair == paste0(mutated.allele.data$Match_Norm_Seq_Allele2,">",mutated.allele.data$Match_Norm_Seq_Allele2),"allele2.pair"] <- ""
mutated.allele.data$allele.both <- paste0(mutated.allele.data$allele1.pair,mutated.allele.data$allele2.pair)
#remove deletions
mutated.allele.data[grep("-",mutated.allele.data$allele.both),"allele.both"] <- ""
mutated.allele.data <- mutated.allele.data[mutated.allele.data$allele.both != "", ]
#copse to pyrimidines only (CT)
mutated.allele.data[mutated.allele.data$allele2.pair == "A>T",] <- "T>A"
mutated.allele.data[mutated.allele.data$allele2.pair == "A>C",] <- "T>G"
mutated.allele.data[mutated.allele.data$allele2.pair == "G>T",] <- "C>A"
mutated.allele.data[mutated.allele.data$allele2.pair == "G>C",] <- "C>G"
mutated.allele.data[mutated.allele.data$allele2.pair == "A>G",] <- "T>C" 
mutated.allele.data[mutated.allele.data$allele2.pair == "G>A",] <- "C>T"

## Prepare data for blotting
# Add Class to mutation data
mutated.allele.data$Cluster <- Consensus.class$Cluster [match(mutated.allele.data$Patient_ID,Consensus.class$Patient_ID)]
mutated.allele.data <- mutated.allele.data [-which(is.na(mutated.allele.data$Cluster)),]
#add subtype to mutation data
mutated.allele.data$Subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(mutated.allele.data$Patient_ID,rownames(ClinicalData.subset))]
#Filter  mutation data Table by IMS
if (IMS.filter == "Luminal") {
  mutated.allele.data.filtered <- mutated.allele.data[mutated.allele.data$Subtype %in% c("Luminal A","Luminal B"),]
} else if (IMS.filter == "Basal"){
  mutated.allele.data.filtered <- mutated.allele.data[mutated.allele.data$Subtype %in% c("Basal-like"),]
} else if (IMS.filter == "Her2"){
  mutated.allele.data.filtered <- mutated.allele.data[mutated.allele.data$Subtype %in% c("HER2-enriched"),]
} else if (IMS.filter == "Luminal A"){
  mutated.allele.data.filtered <- mutated.allele.data[mutated.allele.data$Subtype %in% c("Luminal A"),]
} else if (IMS.filter == "Luminal B"){
  mutated.allele.data.filtered <- mutated.allele.data[mutated.allele.data$Subtype %in% c("Luminal B"),]
} else{
  mutated.allele.data.filtered <- mutated.allele.data
}

#Filter  Mmutation data Table by GOF
if (GOF != "FALSE") {
  mutated.allele.data.filtered <- mutated.allele.data.filtered[mutated.allele.data.filtered$Hugo_Symbol == GOF,]
} else { GOF = ""}


# Select blot data
mutated.allele.blot <- mutated.allele.data.filtered[,c("Patient_ID","Cluster","Subtype","Hugo_Symbol","allele.both")]
mutated.allele.blot.count <- count(mutated.allele.blot,vars=c("allele.both","Cluster"))
mutations.per.cluster <- as.data.frame(table(mutated.allele.blot$Cluster))
colnames(mutations.per.cluster) <- c("Cluster","Count")
mutated.allele.blot.count$cluster.Count <- mutations.per.cluster$Count[match(mutated.allele.blot.count$Cluster,mutations.per.cluster$Cluster)]
mutated.allele.blot.count$frequency <- mutated.allele.blot.count$freq / mutated.allele.blot.count$cluster.Count * 100


# plot
dir.create("./4 FIGURES/mutationspectrum/",showWarnings = FALSE)
tittle = paste0(GOF,".mutation.spectrum.from.",IMS.filter,".bycluster")
png (filename = paste0("./4 FIGURES/mutationspectrum/",tittle,".stacked.png") , height = 600, width= 800)
gg = ggplot(mutated.allele.blot.count, aes(x = Cluster,y=frequency ,fill = allele.both))  +
  geom_bar(stat = "identity", width = 0.8)
gg = gg + ggtitle(tittle) +
  theme(plot.title = element_text(size = 15, lineheight=5, face="bold"))
gg = gg + geom_text(data = mutated.allele.blot.count,
                    aes(y = max(frequency)/10, label = paste0("Tot.Count = ",cluster.Count)),
                    size = 4, vjust = 1.2)
print(gg)
dev.off()

# by sample
mutated.allele.blot.bysamp <- count(mutated.allele.blot,vars=c("allele.both","Patient_ID","Cluster"))
mutations.per.sample <- as.data.frame(table(mutated.allele.blot$Patient_ID))
colnames(mutations.per.sample) <- c("Patien_ID","Count")
mutated.allele.blot.bysamp$sample.Count <- mutations.per.sample$Count[match(mutated.allele.blot.bysamp$Patient_ID,mutations.per.sample$Patien_ID)]
mutated.allele.blot.bysamp$frequency <- mutated.allele.blot.bysamp$freq / mutated.allele.blot.bysamp$sample.Count * 100
# ICR 1, 4 only
mutated.allele.blot.bysamp <- mutated.allele.blot.bysamp[mutated.allele.blot.bysamp$Cluster== "ICR4",]

#order by mutation load
mutated.allele.blot.bysamp <- mutated.allele.blot.bysamp[order(-mutated.allele.blot.bysamp$sample.Count),]

# plot
dir.create("./4 FIGURES/mutationspectrum/",showWarnings = FALSE)
tittle = paste0(GOF,".mutation.spectrum.from.",IMS.filter,".frequency.bysample.ICR4.")
png (filename = paste0("./4 FIGURES/mutationspectrum/",tittle,".stacked.png") , height = 600, width= 2000)
gg = ggplot(mutated.allele.blot.bysamp, aes(x = reorder(Patient_ID,-sample.Count),y=frequency ,fill = allele.both))  +
  geom_bar(stat = "identity", width = 0.95)
gg = gg + ggtitle(tittle) +
  theme(plot.title = element_text(size = 15, lineheight=5, face="bold")) +
  theme(axis.text=element_blank())
#gg = gg + geom_text(data = mutated.allele.blot.bysamp,
#                    aes(y = max(frequency)/10, label = paste0("Tot.Count = ",sample.Count)),
#                   size = 4, vjust = 1.2)
#gg= gg + facet_grid(.~Cluster) +
print(gg)
dev.off()


