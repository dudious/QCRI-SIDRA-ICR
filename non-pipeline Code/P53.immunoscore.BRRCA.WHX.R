#################################################################
###
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
 #Dependencies
 required.packages <- c("ggplot2", "plyr","Hmisc")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")
library("Hmisc")

## Load Data

Mutation.data <- read.csv ("./2 DATA/TCGA Mutations WUSM/BRCA/BRCA.mutation.TCGA.txt",sep ="\t")
Patient_ID <- substr(Mutation.data$Tumor_Sample_Barcode,1,12)
Mutation.selected.data <- cbind(Patient_ID,Mutation.data[,c("Hugo_Symbol","Variant_Classification","Variant_Type")])
Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Hugo_Symbol == "TP53",]
dim (Mutation.selected.data) #192 mutations
Mutation.selected.data <- unique (Mutation.selected.data) #190 mutations
  
Consensus.class <- read.csv("./3 ANALISYS/CLUSTERING/RNAseq/INES/ClusterAssigmentsk4RNASeq.ICR.csv",header=FALSE) # select source data
colnames (Consensus.class) <- c("PatientID","Group")
rownames(Consensus.class) <- Consensus.class[,1]

load (paste0("./2 DATA/SUBSETS/RNASeq.subset.ISGS.RData"))
RNASeq.subset <- as.matrix(RNASeq.subset)

ClinicalData.subset <- read.csv ("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.RNASeq_subset_clinicaldata.csv")                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL

# Calculate Immunposcore
RNASeq.subset.scaled <- scale (RNASeq.subset,scale=FALSE)
RNASeq.subset.scaled <- cbind(RNASeq.subset.scaled,rowMeans(RNASeq.subset.scaled[,-ncol(RNASeq.subset.scaled)]))
colnames(RNASeq.subset.scaled)[ncol(RNASeq.subset.scaled)] <- c("avg")
immunoscore <- RNASeq.subset.scaled[,c("avg"),drop=FALSE]

# Add Class to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,Consensus.class["Group"],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.selected.data) <- Mutation.selected.data$Row.names
Mutation.selected.data$Row.names <- NULL

# add immunoscore to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,immunoscore[,"avg",drop=FALSE],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.selected.data) <- Mutation.selected.data$Row.names
Mutation.selected.data$Row.names <- NULL

# add IMS to mutation data
Mutation.selected.data <- merge(Mutation.selected.data,ClinicalData.subset[,c("TCGA.PAM50.RMethod.MA","TCGA.PAM50.RMethod.RNASeq")],by.x="Patient_ID",by.y="row.names",all.x=TRUE, all.y=FALSE)
row.names(Mutation.selected.data) <- Mutation.selected.data$Row.names
Mutation.selected.data$Row.names <- NULL

# order by ICR group and immunoscore 
Mutation.selected.data <- Mutation.selected.data[order(factor(Mutation.selected.data$Group,levels = c("ICR4","ICR3","ICR2","ICR1")),-Mutation.selected.data$avg),]    # if sorting withing cluster add ,RNASeq.subset$avg 

# split mutation types
Mutation.Missense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Missense_Mutation",]
rownames(Mutation.Missense) <- Mutation.Missense$Mutation_ID
Mutation.Missense$Mutation_ID <- NULL
Mutation.Missense$Variant_Classification <- NULL
Mutation.Missense$Variant_Type <- NULL

#Mutation.Silent <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Silent",]
#rownames(Mutation.Silent) <- Mutation.Silent$Mutation_ID
#Mutation.Silent$Mutation_ID <- NULL
#Mutation.Silent$Variant_Classification <- NULL
#Mutation.Silent$Variant_Type <- NULL

Mutation.Nonsense <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Nonsense_Mutation",]
rownames(Mutation.Nonsense) <- Mutation.Nonsense$Mutation_ID
Mutation.Nonsense$Mutation_ID <- NULL
Mutation.Nonsense$Variant_Classification <- NULL
Mutation.Nonsense$Variant_Type <- NULL

Mutation.other <- Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","RNA","Splice_Site","Silent"),]
#Mutation.other[Mutation.other$Mutation_ID %in% Mutation.other[which(duplicated (Mutation.other$Mutation_ID)),"Mutation_ID"],]
#rownames(Mutation.other) <- Mutation.other$Mutation_ID
Mutation.other$Variant_Classification <- NULL
Mutation.other$Variant_Type <- NULL
Mutation.other <- unique (Mutation.other)
rownames(Mutation.other) <- Mutation.other$Mutation_ID
Mutation.other$Mutation_ID <- NULL

Mutation.any <- Mutation.selected.data
dim (Mutation.any) #192
Mutation.any$Variant_Classification <- NULL
Mutation.any$Variant_Type <- NULL
Mutation.any <- unique (Mutation.any)
rownames(Mutation.any) <- Mutation.any$Mutation_ID
dim (Mutation.any) #185
Mutation.any$Mutation_ID <- NULL

#TP53 imunoscore by cluster and by type of mutation

TP53.IS.any <- aggregate (Mutation.any[,c("avg"),drop=FALSE],by=list(Mutation.any$Group),FUN=mean)
colnames (TP53.IS.any) <- c("Cluster","TP53.Any")
TP53.IS.Missense <- aggregate (Mutation.Missense[,c("avg"),drop=FALSE],by=list(Mutation.Missense$Group),FUN=mean)
colnames (TP53.IS.Missense) <- c("Cluster","TP53.Missense")
TP53.IS.Nonsense <- aggregate (Mutation.Nonsense[,c("avg"),drop=FALSE],by=list(Mutation.Nonsense$Group),FUN=mean)
colnames (TP53.IS.Nonsense) <- c("Cluster","TP53.Nonsense")
#TP53.IS.Silent <- aggregate (Mutation.Silent[,c("avg"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean) (only 3 silent mutations)
#colnames (TP53.IS.Silent) <- c("Cluster","TP53.Silent")
TP53.IS.other <- aggregate (Mutation.other[,c("avg"),drop=FALSE],by=list(Mutation.other$Group),FUN=mean)
colnames (TP53.IS.other) <- c("Cluster","TP53.other")

#Wild type IS

WT.IS <- immunoscore[-which(rownames(immunoscore)  %in% Mutation.any$Patient_ID),,drop=FALSE]
WT.IS <- merge(WT.IS,Consensus.class["Group"],by.x="row.names",by.y="row.names",all.x=TRUE, all.y=FALSE)
rownames (WT.IS) <- WT.IS$Row.names
WT.IS$Row.names <- NULL
WT.IS.avg <- aggregate (WT.IS[,1,drop=FALSE],by=list(WT.IS$Group),FUN=mean)
colnames (WT.IS.avg) <- c("Cluster","TP53.WT")

#NON-Missense IS

NMM.IS <- immunoscore[-which(rownames(immunoscore)  %in% Mutation.Missense$Patient_ID),,drop=FALSE]
NMM.IS <- merge(NMM.IS,Consensus.class["Group"],by.x="row.names",by.y="row.names",all.x=TRUE, all.y=FALSE)
rownames (NMM.IS) <- NMM.IS$Row.names
NMM.IS$Row.names <- NULL
NMM.IS.avg <- aggregate (NMM.IS[,1,drop=FALSE],by=list(NMM.IS$Group),FUN=mean)
colnames (NMM.IS.avg) <- c("Cluster","TP53.NMM")

TP53.IS <- cbind (TP53.IS.any,TP53.IS.Missense,TP53.IS.Nonsense,TP53.IS.other,WT.IS.avg,NMM.IS.avg)
TP53.IS <- TP53.IS[,-c(3,5,7,9,11)]

#TP53 imunoscore by IMS and by type of mutation

TP53.byIMS.IS.any <- aggregate (Mutation.any[,c("avg"),drop=FALSE],by=list(Mutation.any$TCGA.PAM50.RMethod.MA),FUN=mean)
colnames (TP53.byIMS.IS.any) <- c("Subtype","TP53.byIMS.Any")
TP53.byIMS.IS.Missense <- aggregate (Mutation.Missense[,c("avg"),drop=FALSE],by=list(Mutation.Missense$TCGA.PAM50.RMethod.MA),FUN=mean)
colnames (TP53.byIMS.IS.Missense) <- c("Subtype","TP53.byIMS.Missense")
TP53.byIMS.IS.Nonsense <- aggregate (Mutation.Nonsense[,c("avg"),drop=FALSE],by=list(Mutation.Nonsense$TCGA.PAM50.RMethod.MA),FUN=mean)
colnames (TP53.byIMS.IS.Nonsense) <- c("Subtype","TP53.byIMS.Nonsense")
#TP53.byIMS.IS.Silent <- aggregate (Mutation.Silent[,c("avg"),drop=FALSE],by=list(Mutation.other$TCGA.PAM50.RMethod.MA),FUN=mean) (only 3 silent mutations)
#colnames (TP53.byIMS.IS.Silent) <- c("ubtype","TP53.byIMS.Silent")
TP53.byIMS.IS.other <- aggregate (Mutation.other[,c("avg"),drop=FALSE],by=list(Mutation.other$TCGA.PAM50.RMethod.MA),FUN=mean)
colnames (TP53.byIMS.IS.other) <- c("Subtype","TP53.byIMS.other")

#Wild type IS byIMS

WT.byIMS.IS <- immunoscore[-which(rownames(immunoscore)  %in% Mutation.any$Patient_ID),,drop=FALSE]
WT.byIMS.IS <- merge(WT.byIMS.IS,ClinicalData.subset[,c("TCGA.PAM50.RMethod.MA","TCGA.PAM50.RMethod.RNASeq")],by.x="row.names",by.y="row.names",all.x=TRUE, all.y=FALSE)
rownames (WT.byIMS.IS) <- WT.byIMS.IS$Row.names
WT.byIMS.IS$Row.names <- NULL
WT.byIMS.IS.avg <- aggregate (WT.byIMS.IS[,1,drop=FALSE],by=list(WT.byIMS.IS$TCGA.PAM50.RMethod.MA),FUN=mean)
colnames (WT.byIMS.IS.avg) <- c("Subtype","TP53.WT.byIMS")

#NON-Missense IS byIMS

NMM.byIMS.IS <- immunoscore[-which(rownames(immunoscore)  %in% Mutation.Missense$Patient_ID),,drop=FALSE]
NMM.byIMS.IS <- merge(NMM.byIMS.IS,ClinicalData.subset[,c("TCGA.PAM50.RMethod.MA","TCGA.PAM50.RMethod.RNASeq")],by.x="row.names",by.y="row.names",all.x=TRUE, all.y=FALSE)
rownames (NMM.byIMS.IS) <- NMM.byIMS.IS$Row.names
NMM.byIMS.IS$Row.names <- NULL
NMM.byIMS.IS.avg <- aggregate (NMM.byIMS.IS[,1,drop=FALSE],by=list(NMM.byIMS.IS$TCGA.PAM50.RMethod.MA),FUN=mean)
colnames (NMM.byIMS.IS.avg) <- c("Subtype","TP53.NMM.byIMS")

TP53.byIMS.IS.other <- rbind(TP53.byIMS.IS.other,c("Normal-like","NA"))
TP53.byIMS.IS <- cbind (TP53.byIMS.IS.any,TP53.byIMS.IS.Missense,TP53.byIMS.IS.Nonsense,TP53.byIMS.IS.other,WT.byIMS.IS.avg,NMM.byIMS.IS.avg)
TP53.byIMS.IS <- TP53.byIMS.IS[,-c(3,5,7,9,11)]


#Export Patient w/ P53 mutation list
Mutation.any <- cbind (Mutation.any,Mutated="TRUE")
Mutation.Missense <- cbind (Mutation.Missense,Missense="TRUE")
TP53.patients <- merge (Mutation.any[,c("Patient_ID","Mutated")],Mutation.Missense[,c("Patient_ID","Missense")],by="Patient_ID", all.x=TRUE)
mutation.freq <- read.csv ("./3 ANALISYS/Mutations/Mutations.TCGA.BRCA.Patient.by.Cluster.csv")                               # mutation frequencies
rownames(mutation.freq) <- mutation.freq$Patient_ID
mutation.freq <-  mutation.freq[,c("Freq.all","Freq.Missense")]
#TP53.patients <-TP53.patients[,-c(2,8,9,10,11,12)]
TP53.patients <- merge (mutation.freq,TP53.patients,by.x="row.names",by.y="Patient_ID",all.x=TRUE)
row.names(TP53.patients) <- TP53.patients$Row.names
TP53.patients$Row.names <- NULL
levels (TP53.patients$Mutated) <- c("TRUE","FALSE")
TP53.patients$Mutated[is.na(TP53.patients$Mutated)] <- "FALSE"
levels (TP53.patients$Missense) <- c("TRUE","FALSE")
TP53.patients$Missense[is.na(TP53.patients$Missense)] <- "FALSE"
TP53.patients <- merge(TP53.patients,immunoscore,by="row.names",all.x=TRUE)
row.names(TP53.patients) <- TP53.patients$Row.names
TP53.patients$Row.names <- NULL
#write.csv (TP53.patients,file=("./3 ANALISYS/Mutations/TP53.patients.csv"))

#export immunoscore
immunoscore <- merge(immunoscore,ClinicalData.subset[,c("TCGA.PAM50.RMethod.MA","TCGA.PAM50.RMethod.RNASeq")],by="row.names",all.x=TRUE, all.y=FALSE)
row.names(immunoscore) <- immunoscore$Row.names
immunoscore$Row.names <- NULL
#write.csv (immunoscore,file=("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BRCA.ISGS.csv"))

IS.byIMS <- aggregate (immunoscore[,c("avg"),drop=FALSE],by=list(immunoscore$TCGA.PAM50.RMethod.RNASeq),FUN=mean)
colnames (IS.byIMS) <- c("Subtype","IS")

IS.byIMS.MUT.subset <- aggregate (immunoscore[unique(Patient_ID),"avg",drop=FALSE],by=list(immunoscore[unique(Patient_ID),"TCGA.PAM50.RMethod.RNASeq"]),FUN=mean)
colnames (IS.byIMS) <- c("Subtype","IS")

WT.Patients <- Patient_ID[which(unique(Patient_ID) %nin% Mutation.any$Patient_ID)]

IS.byTP53.AMT <- aggregate (TP53.patients[,c("avg"),drop=FALSE],by=list(TP53.patients$Mutated),FUN=mean)
colnames (IS.byTP53.AMT) <- c("Mutated","IS")

IS.byTP53.MMT <- aggregate (TP53.patients[,c("avg"),drop=FALSE],by=list(TP53.patients$Missense),FUN=mean)
colnames (IS.byTP53.MMT) <- c("Mutated","IS")



