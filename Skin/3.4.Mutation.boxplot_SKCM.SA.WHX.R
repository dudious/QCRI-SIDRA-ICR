#################################################################
###
### This Script PLots boxplots of mutation frequency by ICR cluster
### based on RNASeq Data Consensus Clustering
### 
### Input data :
### ./3 ANALISYS/CLUSTERING/RNAseq/...
### ./2 DATA/TCGA Mutations/...
### Data is saved :
### ./3 ANALISYS/Mutations/Mutations.Gene.by.Cluster.csv
### ./3 ANALISYS/Mutations/Mutations.Patient.by.Cluster.csv
### Figures are saved :
### ./4 FIGURES/Mutation Plots
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
 #Dependencies
 required.packages <- c("ggplot2", "plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")
#dir.create("./3 ANALISYS/Mutations/SKCM")
## Load Data
## download data from TCGA site (what dataset)
Mutation.data <- read.csv ("./2 DATA/TCGA Mutations/SKCM/Somatic_Mutations/BI__IlluminaGA_DNASeq_automated/Level_2/broad.mit.edu__IlluminaGA_automated_DNA_sequencing_level2.maf",sep ="\t")
Mutation.selected.data <- data.frame(Hugo_Symbol = Mutation.data$Hugo_Symbol, Variant_Classification = Mutation.data$Variant_Classification, Patient_ID = substr(Mutation.data$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins
Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/SKCM/SKCM.TCGA.EDASeq.k7.ISGS.reps5000/SKCM.TCGA.EDASeq.k7.ISGS.reps5000.k=4.consensusClass.ICR.csv",header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

# Add Class to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,Consensus.class,by="Patient_ID",all.x=TRUE, all.y=FALSE)
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Cluster)),]

# split mutation types
Mutation.Missense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Missense_Mutation",]
Mutation.Silent <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Silent",]
Mutation.Nonsense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Nonsense_Mutation",]
Mutation.Other <- Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","RNA","Splice_Site"),]
Mutation.Any <- unique (Mutation.selected.data[c("Patient_ID","Hugo_Symbol","Cluster")])
Mutation.All <- Mutation.selected.data
save (Mutation.All,Mutation.Missense,Mutation.Silent,Mutation.Nonsense,Mutation.Other,Mutation.Any,file="./3 ANALISYS/Mutations/SKCM/Mutation.Data.split.RDATA")

#Gene mutation frequency by Cluster table
count.gene.bycluster <- function (Mutation.x){
  Count.Gene <- count(Mutation.x[,c("Hugo_Symbol","Cluster")])
  Count.Gene <- Count.Gene [order(-as.numeric(Count.Gene$freq)),]
  rownames(Count.Gene) <- paste0(Count.Gene$Hugo_Symbol,"_",Count.Gene$Cluster)
  return (Count.Gene)
}
Count.Missense.Gene <- count.gene.bycluster (Mutation.Missense)
Count.Silent.Gene <- count.gene.bycluster (Mutation.Silent)
Count.Nonsense.Gene <- count.gene.bycluster (Mutation.Nonsense)
Count.Other.Gene <- count.gene.bycluster (Mutation.Other)
Count.Any.Gene <- count.gene.bycluster (Mutation.Any)           # unique variants classification
Count.All.Gene <- count.gene.bycluster (Mutation.All) # raw count 

Mutation.Frequency.Gene <- Count.All.Gene
colnames(Mutation.Frequency.Gene)[3] <- "Freq.All"
Mutation.Frequency.Gene$Freq.Any <- Count.Any.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Any.Gene))]
Mutation.Frequency.Gene$Freq.Missense <- Count.Missense.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Missense.Gene))]
Mutation.Frequency.Gene$Freq.Nonsense <- Count.Nonsense.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Nonsense.Gene))]
Mutation.Frequency.Gene$Freq.Silent <- Count.Silent.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Silent.Gene))]
Mutation.Frequency.Gene$Freq.Other <- Count.Other.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Other.Gene))]

write.csv (Mutation.Frequency.Gene,file="./3 ANALISYS/Mutations/SKCM/Mutations.TCGA.SKCM.Gene.by.Cluster.csv")

#Patient mutation frequency by Cluster table
count.Patient.bycluster <- function (Mutation.x){
  Count.Patient <- count(Mutation.x[,c("Patient_ID","Cluster")])
  Count.Patient <- Count.Patient [order(-as.numeric(Count.Patient$freq)),]
  rownames(Count.Patient) <- paste0(Count.Patient$Patient_ID,"_",Count.Patient$Cluster)
  return (Count.Patient)
}
Count.Missense.Patient <- count.Patient.bycluster (Mutation.Missense)
Count.Silent.Patient <- count.Patient.bycluster (Mutation.Silent)
Count.Nonsense.Patient <- count.Patient.bycluster (Mutation.Nonsense)
Count.Other.Patient <- count.Patient.bycluster (Mutation.Other)
Count.Any.Patient <- count.Patient.bycluster (Mutation.Any)           # unique variants classification
Count.All.Patient <- count.Patient.bycluster (Mutation.All) # raw count 

Mutation.Frequency.Patient <- Count.All.Patient
colnames(Mutation.Frequency.Patient)[3] <- "Freq.All"
Mutation.Frequency.Patient$Freq.Any <- Count.Any.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Any.Patient))]
Mutation.Frequency.Patient$Freq.Missense <- Count.Missense.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Missense.Patient))]
Mutation.Frequency.Patient$Freq.Nonsense <- Count.Nonsense.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Nonsense.Patient))]
Mutation.Frequency.Patient$Freq.Silent <- Count.Silent.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Silent.Patient))]
Mutation.Frequency.Patient$Freq.Other <- Count.Other.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Other.Patient))]

write.csv (Mutation.Frequency.Patient,file="./3 ANALISYS/Mutations/SKCM/Mutations.TCGA.SKCM.Patient.by.Cluster.csv")

#Prepare Data for Boxplots

numMuts.SGall <- data.frame(count=Count.All.Patient$freq,cluster=Count.All.Patient$Cluster,mut.type = "All")

numMuts.Missense <- cbind(Count.Missense.Patient[,c("freq","Cluster")],mut.type = "Missense")
numMuts.Nonsense <- cbind(Count.Nonsense.Patient[,c("freq","Cluster")],mut.type = "Nonsense")
numMuts.Silent <- cbind(Count.Silent.Patient[,c("freq","Cluster")],mut.type = "Silent")
numMuts.Other <- cbind(Count.Other.Patient[,c("freq","Cluster")],mut.type = "Other")
numMuts.SGtype <- rbind(numMuts.Missense,numMuts.Silent,numMuts.Nonsense,numMuts.Other)
rownames(numMuts.SGtype) <- NULL
colnames(numMuts.SGtype) = c( "count", "cluster", "mut.type")


# Choose the types of mutations
#muts.count.type = count(Mutation.selected.data[,c("Patient_ID","Variant_Classification","Cluster")])
#numMuts.SGtype = data.frame(count=muts.count.type$freq, cluster= muts.count.type$Cluster, mut.type = muts.count.type$Variant_Classification)
#numMuts.SGtype = numMuts.SGtype[numMuts.SGtype$mut.type %in% c("Missense","Silent","Nonsense","Other"),]

# Combine
numMuts.SGall = rbind(numMuts.SGall, numMuts.SGtype)

meds <- ddply(numMuts.SGall, .(mut.type, cluster), summarize, med = median(count)) ## median
mean.n <- function(x){ return(c(y = 0 , label = round(mean(x),2))) } ## mean

png("./4 FIGURES/Mutation Plots/Mutations.TCGA.SKCM.All.SA.WHX.png", height = 1000, width= 1000)   #set filename
  cluster.order = c("ICR4", "ICR3", "ICR2", "ICR1")
  colors = c("blue", "green", "orange", "red")
  gg = ggplot(numMuts.SGall, aes(cluster, count, fill=cluster)) +
        stat_boxplot(geom ='errorbar') +
        geom_boxplot() +
        geom_jitter(position=position_jitter(width=0.1,height=0.1))
        #, aes(color="lightgray")) 
  gg = gg + ylab("Number of mutations per sample") +
            scale_x_discrete(limits=cluster.order) +
            facet_grid(.~mut.type,
                      scales = "free",
                      space="free") +
           xlab("Clusters") + theme_bw() 
  gg = gg + scale_fill_manual(values = colors) +
            scale_y_continuous(breaks = seq(0, 300, 100), limits = c(0, 250)) 
  gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                  legend.position = "none",
                  axis.text.x = element_text(size = 12, vjust=1),
                  axis.title.x = element_text(size = 18, vjust = -1),
                  axis.text.y = element_text(size = 18, vjust=1),
                  axis.title.y = element_text(size = 18, vjust = 1))
  gg = gg + geom_text(data = meds,
                      aes(y = 0, label = round(med,2)),
                      size = 7, vjust = 1.2)
  gg = gg + stat_summary(fun.data = mean.n,
                         geom = "text",
                         fun.y = mean,
                         colour = "black",
                         vjust = 3.7,
                         size = 7)
print(gg)
dev.off()


