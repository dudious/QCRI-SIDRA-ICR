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

## Load Data

Mutation.data <- read.csv ("./2 DATA/TCGA Mutations WUSM/LIHC/LIHC.mutation.TCGA.txt",sep ="\t")
Patient_ID <- substr(Mutation.data$Tumor_Sample_Barcode,1,12)
Mutation.selected.data <- cbind(Patient_ID,Mutation.data[,c("Hugo_Symbol","Variant_Classification","Variant_Type")])
dim (Mutation.selected.data) #33963 mutations
Mutation.selected.data <- unique (Mutation.selected.data) #33399 mutations
Mutation_ID <- paste0(Mutation.selected.data$Patient_ID,"_",Mutation.selected.data$Hugo_Symbol)
Mutation.selected.data <- cbind(Mutation_ID,Mutation.selected.data)
Mutation.selected.data <- Mutation.selected.data[order(Mutation.selected.data$Mutation_ID),]
  
dim(Mutation.selected.data[which(duplicated (Mutation.selected.data$Mutation_ID)),]) #390 doubles
#Mutation.selected.data[Mutation.selected.data$Mutation_ID %in% Mutation.selected.data[which(duplicated (Mutation.selected.data$Mutation_ID)),"Mutation_ID"],]

Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/LIHC.TCGA.EDASeq.k4.ISGS.reps2000/LIHC.TCGA.EDASeq.k4.ISGS.reps2000.k=4.consensusClass.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]

# Add Class to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,Consensus.class["Group"],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.selected.data) <- Mutation.selected.data$Row.names
Mutation.selected.data$Row.names <- NULL

# split mutation types
Mutation.Missense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Missense_Mutation",]
#Mutation.Missense[Mutation.Missense$Mutation_ID %in% Mutation.Missense[which(duplicated (Mutation.Missense$Mutation_ID)),"Mutation_ID"],]
Mutation.Missense$Variant_Classification <- NULL
Mutation.Missense$Variant_Type <- NULL
Mutation.Missense <- unique(Mutation.Missense)
rownames(Mutation.Missense) <- Mutation.Missense$Mutation_ID
Mutation.Missense$Mutation_ID <- NULL

Mutation.Silent <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Silent",]
rownames(Mutation.Silent) <- Mutation.Silent$Mutation_ID
Mutation.Silent$Mutation_ID <- NULL
Mutation.Silent$Variant_Classification <- NULL
Mutation.Silent$Variant_Type <- NULL

Mutation.Nonsense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Nonsense_Mutation",]
rownames(Mutation.Nonsense) <- Mutation.Nonsense$Mutation_ID
Mutation.Nonsense$Mutation_ID <- NULL
Mutation.Nonsense$Variant_Classification <- NULL
Mutation.Nonsense$Variant_Type <- NULL

Mutation.other <- Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","RNA","Splice_Site"),]
#Mutation.other[Mutation.other$Mutation_ID %in% Mutation.other[which(duplicated (Mutation.other$Mutation_ID)),"Mutation_ID"],]
#rownames(Mutation.other) <- Mutation.other$Mutation_ID
Mutation.other$Variant_Classification <- NULL
Mutation.other$Variant_Type <- NULL
Mutation.other <- unique (Mutation.other)
rownames(Mutation.other) <- Mutation.other$Mutation_ID
Mutation.other$Mutation_ID <- NULL

Mutation.any <- Mutation.selected.data
dim (Mutation.any) #27482
Mutation.any$Variant_Classification <- NULL
Mutation.any$Variant_Type <- NULL
Mutation.any <- unique (Mutation.any)
rownames(Mutation.any) <- Mutation.any$Mutation_ID
dim (Mutation.any) #27092
Mutation.any$Mutation_ID <- NULL

#Gene mutation frequency by Cluster table
Count.Missense.Gene <- count(Mutation.Missense[,c("Hugo_Symbol","Group")])
Count.Missense.Gene <- Count.Missense.Gene [order(-as.numeric(Count.Missense.Gene$freq)),]
rownames(Count.Missense.Gene) <- paste0(Count.Missense.Gene$Hugo_Symbol,"_",Count.Missense.Gene$Group)

Count.Silent.Gene <- count(Mutation.Silent[,c("Hugo_Symbol","Group")])
Count.Silent.Gene <- Count.Silent.Gene [order(-as.numeric(Count.Silent.Gene$freq)),]
rownames(Count.Silent.Gene) <- paste0(Count.Silent.Gene$Hugo_Symbol,"_",Count.Silent.Gene$Group)

Count.Nonsense.Gene <- count(Mutation.Nonsense[,c("Hugo_Symbol","Group")])
Count.Nonsense.Gene <- Count.Nonsense.Gene [order(-as.numeric(Count.Nonsense.Gene$freq)),]
rownames(Count.Nonsense.Gene) <- paste0(Count.Nonsense.Gene$Hugo_Symbol,"_",Count.Nonsense.Gene$Group)

Count.All.Gene <- count(Mutation.any[,c("Hugo_Symbol","Group")])
Count.All.Gene <- Count.All.Gene [order(-as.numeric(Count.All.Gene$freq)),]
rownames(Count.All.Gene) <- paste0(Count.All.Gene$Hugo_Symbol,"_",Count.All.Gene$Group)

Mutation.Frequency.Gene <- Count.All.Gene
Mutation.Frequency.Gene <- merge (Mutation.Frequency.Gene,Count.Missense.Gene[,"freq",drop=FALSE],by="row.names",all.x=TRUE,all.y=FALSE)
rownames(Mutation.Frequency.Gene) <- Mutation.Frequency.Gene$Row.names
Mutation.Frequency.Gene$Row.names <- NULL
Mutation.Frequency.Gene <- merge (Mutation.Frequency.Gene,Count.Nonsense.Gene[,"freq",drop=FALSE],by="row.names",all.x=TRUE,all.y=FALSE)
rownames(Mutation.Frequency.Gene) <- Mutation.Frequency.Gene$Row.names
Mutation.Frequency.Gene$Row.names <- NULL
Mutation.Frequency.Gene <- merge (Mutation.Frequency.Gene,Count.Silent.Gene[,"freq",drop=FALSE],by="row.names",all.x=TRUE,all.y=FALSE)
rownames(Mutation.Frequency.Gene) <- Mutation.Frequency.Gene$Row.names
Mutation.Frequency.Gene$Row.names <- NULL
colnames (Mutation.Frequency.Gene) = c("HUGO.Symbol","RNASeq.Cluster","Freq.all","Freq.Missense","Freq.Nonsense","Freq.Silent")
Mutation.Frequency.Gene <- Mutation.Frequency.Gene[order(-as.numeric(Mutation.Frequency.Gene$Freq.all)),]

write.csv (Mutation.Frequency.Gene,file="./3 ANALISYS/Mutations/Mutations.TCGA.LIHC.Gene.by.Cluster.csv")

#Patient mutation frequency by Cluster table
Count.Missense.Patient <- count(Mutation.Missense[,c("Patient_ID","Group")])
Count.Missense.Patient <- Count.Missense.Patient [order(-as.numeric(Count.Missense.Patient$freq)),]
rownames(Count.Missense.Patient) <- paste0(Count.Missense.Patient$Patient_ID,"_",Count.Missense.Patient$Group)

Count.Silent.Patient <- count(Mutation.Silent[,c("Patient_ID","Group")])
Count.Silent.Patient <- Count.Silent.Patient [order(-as.numeric(Count.Silent.Patient$freq)),]
rownames(Count.Silent.Patient) <- paste0(Count.Silent.Patient$Patient_ID,"_",Count.Silent.Patient$Group)

Count.Nonsense.Patient <- count(Mutation.Nonsense[,c("Patient_ID","Group")])
Count.Nonsense.Patient <- Count.Nonsense.Patient [order(-as.numeric(Count.Nonsense.Patient$freq)),]
rownames(Count.Nonsense.Patient) <- paste0(Count.Nonsense.Patient$Patient_ID,"_",Count.Nonsense.Patient$Group)

Count.Other.Patient <- count(Mutation.other[,c("Patient_ID","Group")])
Count.Other.Patient <- Count.Other.Patient [order(-as.numeric(Count.Other.Patient$freq)),]
rownames(Count.Other.Patient) <- paste0(Count.Other.Patient$Patient_ID,"_",Count.Other.Patient$Group)

Count.All.Patient <- count(Mutation.any[,c("Patient_ID","Group")])
Count.All.Patient <- Count.All.Patient [order(-as.numeric(Count.All.Patient$freq)),]
rownames(Count.All.Patient) <- paste0(Count.All.Patient$Patient_ID,"_",Count.All.Patient$Group)

Mutation.Frequency.Patient <- Count.All.Patient
Mutation.Frequency.Patient <- merge (Mutation.Frequency.Patient,Count.Missense.Patient[,"freq",drop=FALSE],by="row.names",all.x=TRUE,all.y=FALSE)
rownames(Mutation.Frequency.Patient) <- Mutation.Frequency.Patient$Row.names
Mutation.Frequency.Patient$Row.names <- NULL
Mutation.Frequency.Patient <- merge (Mutation.Frequency.Patient,Count.Nonsense.Patient[,"freq",drop=FALSE],by="row.names",all.x=TRUE,all.y=FALSE)
rownames(Mutation.Frequency.Patient) <- Mutation.Frequency.Patient$Row.names
Mutation.Frequency.Patient$Row.names <- NULL
Mutation.Frequency.Patient <- merge (Mutation.Frequency.Patient,Count.Silent.Patient[,"freq",drop=FALSE],by="row.names",all.x=TRUE,all.y=FALSE)
rownames(Mutation.Frequency.Patient) <- Mutation.Frequency.Patient$Row.names
Mutation.Frequency.Patient$Row.names <- NULL
colnames (Mutation.Frequency.Patient) = c("Patient_ID","RNASeq.Cluster","Freq.all","Freq.Missense","Freq.Nonsense","Freq.Silent")
Mutation.Frequency.Patient <- Mutation.Frequency.Patient[order(-as.numeric(Mutation.Frequency.Patient$Freq.all)),]

write.csv (Mutation.Frequency.Patient,file="./3 ANALISYS/Mutations/Mutations.TCGA.LIHC.Patient.by.Cluster.csv")

#Prepare Data for Boxplots

numMuts.SGall <- cbind(Count.All.Patient[,c("freq","Group")],mut.type = "All")
rownames(numMuts.SGall) <- NULL
colnames(numMuts.SGall) = c( "count", "cluster", "mut.type")

numMuts.Missense <- cbind(Count.Missense.Patient[,c("freq","Group")],mut.type = "Missense")
numMuts.Nonsense <- cbind(Count.Nonsense.Patient[,c("freq","Group")],mut.type = "Nonsense")
numMuts.Silent <- cbind(Count.Silent.Patient[,c("freq","Group")],mut.type = "Silent")
numMuts.Other <- cbind(Count.Other.Patient[,c("freq","Group")],mut.type = "Other")
numMuts.SGtype <- rbind(numMuts.Missense,numMuts.Silent,numMuts.Nonsense,numMuts.Other)
rownames(numMuts.SGtype) <- NULL
colnames(numMuts.SGtype) = c( "count", "cluster", "mut.type")


# Choose the types of mutations
numMuts.SGtype = numMuts.SGtype[numMuts.SGtype$mut.type %in% c("Missense","Silent","Nonsense","Other"),]

# Combine
numMuts.SGall = rbind(numMuts.SGall, numMuts.SGtype)

meds <- ddply(numMuts.SGall, .(mut.type, cluster), summarize, med = median(count)) ## median
mean.n <- function(x){ return(c(y = 0 , label = round(mean(x),2))) } ## mean

png("./4 FIGURES/Mutation Plots/Mutations.TCGA.LIHC.All.SA.WHX.png", height = 1000, width= 1000)   #set filename
  cluster.order = c("ICR4", "ICR3", "ICR2", "ICR1")
  colors = c("blue", "green", "orange", "red")
  gg = ggplot(numMuts.SGall, aes(cluster, count, fill=cluster)) +
        stat_boxplot(geom ='errorbar') +
        geom_boxplot() +
        geom_jitter(position=position_jitter())
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


