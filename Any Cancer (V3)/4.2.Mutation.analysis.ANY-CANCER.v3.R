#################################################################
###
### This Script PLots boxplots of mutation frequency by ICR cluster
### based on ",Cancerset," RNASeq Data Consensus Clustering
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
#setwd("~/Dropbox/BREAST_QATAR/")
setwd("/mnt3/wouter/BREAST-QATAR/")
 #Dependencies
 required.packages <- c("plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
library("plyr")

# Set Parameters
Cancersets <- "ALL"             # do not use -GA or -hiseq (data is merged)
BRCA.Filter <- "BSF2"         # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
IMS.filter = "All"            # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
K <- 4                        # SET K here
DL.Method    = "BIOLINKS"     #Choose "ASSEMBLER" or "BIOLINKS"

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}
  Parent.Cancerset <- substring(Cancerset,1,4)

## Load Data
if (DL.Method == "ASSEMBLER") {
### ASSEMBLER
if (Cancerset %in% c("COAD","READ","UCEC")) {
  #GA data
  Cancerset <- paste0(Cancerset,"-GA")
  load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
  Mutation.selected.data.GA <- data.frame(Hugo_Symbol = maf.merged.table$Hugo_Symbol, Variant_Classification = maf.merged.table$Variant_Classification, Patient_ID = substr(maf.merged.table$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins
  Consensus.class.GA <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class.GA <- Consensus.class.GA[,-1]
  colnames (Consensus.class.GA) <- c("Patient_ID","Cluster")
  rownames(Consensus.class.GA) <- Consensus.class.GA[,1]
  rm(maf.merged.table)
  Cancerset <- substring(Cancerset,1,4)
  #hiseq data
  Cancerset <- paste0(Cancerset,"-hiseq")
  load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
  Mutation.selected.data.hiseq <- data.frame(Hugo_Symbol = maf.merged.table$Hugo_Symbol, Variant_Classification = maf.merged.table$Variant_Classification, Patient_ID = substr(maf.merged.table$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins
  Consensus.class.hiseq <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
  Consensus.class.hiseq <- Consensus.class.hiseq[,-1]
  colnames (Consensus.class.hiseq) <- c("Patient_ID","Cluster")
  rownames(Consensus.class.hiseq) <- Consensus.class.hiseq[,1]
  rm(maf.merged.table)
  Cancerset <- substring(Cancerset,1,4)
  #merge GA-hiseq
  Consensus.class <- rbind (Consensus.class.hiseq,Consensus.class.GA)
  Mutation.selected.data <- rbind (Mutation.selected.data.hiseq,Mutation.selected.data.GA)
  Consensus.class <- Consensus.class[which(rownames(Consensus.class) %in% as.character(Mutation.selected.data$Patient_ID)),] # drop samples without any mutation data
} else {
  load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
  Mutation.selected.data <- data.frame(Hugo_Symbol = maf.merged.table$Hugo_Symbol, Variant_Classification = maf.merged.table$Variant_Classification, Patient_ID = substr(maf.merged.table$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins
  if (Cancerset == "BRCA"){
    if (substring(Geneset,7,10)=="FLTR"){
    Cancerset <- paste0(Cancerset,".",BRCA.Filter)
    }
  }

 Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
 Consensus.class <- Consensus.class[,-1]
 colnames (Consensus.class) <- c("Patient_ID","Cluster")
 rownames(Consensus.class) <- Consensus.class[,1]
 Consensus.class <- Consensus.class[which(rownames(Consensus.class) %in% as.character(Mutation.selected.data$Patient_ID)),] # drop samples without any mutation data
 rm(maf.merged.table)
}}

### BIOLINKS
if (DL.Method == "BIOLINKS") {
load (file=paste0("./2 DATA/TCGA Mutations/BIOLINKS/",Cancerset,".MAF.RData"))
Mutation.selected.data <- data.frame(Hugo_Symbol = mut$Hugo_Symbol, Variant_Classification = mut$Variant_Classification, Patient_ID = substr(mut$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins
Consensus.class <- read.csv(file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]
Consensus.class <- Consensus.class[which(rownames(Consensus.class) %in% as.character(Mutation.selected.data$Patient_ID)),] # drop samples without any mutation data
}

dir.create (paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Cancerset,"/"),showWarnings=FALSE)

# Add Class to mutation data
Mutation.selected.data$Cluster <- Consensus.class$Cluster [match(Mutation.selected.data$Patient_ID,Consensus.class$Patient_ID)]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Cluster)),]

Mutation.gene.freq <- count(unique(Mutation.selected.data[,c("Patient_ID","Hugo_Symbol")]),vars=c("Hugo_Symbol"))
Mutation.gene.freq$PCTofN <- round(Mutation.gene.freq$freq/length(unique(Mutation.selected.data$Patient_ID))*100,2)
nrow(Mutation.gene.freq[Mutation.gene.freq$PCTofN>0.5,])
write.csv(Mutation.gene.freq,file=paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".Gene.mutation.frequency.csv"),row.names = FALSE)

# split mutation types
# Raw count
Mutation.All       <- Mutation.selected.data
Mutation.Missense  <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Missense_Mutation",]
Mutation.Silent    <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Silent",]
Mutation.Nonsense  <- Mutation.selected.data[Mutation.selected.data$Variant_Classification == "Nonsense_Mutation",]
Mutation.Other     <- Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","RNA","Splice_Site"),]
Mutation.NonSilent <- Mutation.selected.data[-which(Mutation.selected.data$Variant_Classification == "Silent"),]
# One mutation / gene
Mutation.Any          <- unique (Mutation.All[c("Patient_ID","Hugo_Symbol","Cluster")])
Mutation.Missense.Any <- unique (Mutation.Missense[c("Patient_ID","Hugo_Symbol","Cluster")])
Mutation.Silent.Any   <- unique (Mutation.Silent[c("Patient_ID","Hugo_Symbol","Cluster")])
Mutation.NonSilent.Any<- unique (Mutation.NonSilent[c("Patient_ID","Hugo_Symbol","Cluster")])

save (Mutation.All,Mutation.Missense,Mutation.Silent,Mutation.Nonsense,Mutation.Other,Mutation.NonSilent,
      Mutation.Any,Mutation.Missense.Any,Mutation.Silent.Any,Mutation.NonSilent.Any,
      file=paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",Geneset,".split.RDATA"))

#### ByGene ####

#Gene mutation frequency by Cluster table
count.gene.bycluster <- function (Mutation.x){
  Count.Gene <- count(Mutation.x[,c("Hugo_Symbol","Cluster")])
  Count.Gene <- Count.Gene [order(-as.numeric(Count.Gene$freq)),]
  rownames(Count.Gene) <- paste0(Count.Gene$Hugo_Symbol,"_",Count.Gene$Cluster)
  return (Count.Gene)
}
Count.All.Gene <- count.gene.bycluster (Mutation.All)           # raw count
Count.Missense.Gene <- count.gene.bycluster (Mutation.Missense)
Count.Silent.Gene <- count.gene.bycluster (Mutation.Silent)
Count.Nonsense.Gene <- count.gene.bycluster (Mutation.Nonsense)
Count.Other.Gene <- count.gene.bycluster (Mutation.Other)
Count.NonSilent.Gene <- count.gene.bycluster (Mutation.NonSilent)
Count.Any.Gene <- count.gene.bycluster (Mutation.Any)           # unique variants classification
Count.Missense.Any.Gene <- count.gene.bycluster (Mutation.Missense.Any)
Count.Silent.Any.Gene <- count.gene.bycluster (Mutation.Silent.Any)    
Count.NonSilent.Any.Gene <- count.gene.bycluster (Mutation.NonSilent.Any)   

#merge the byGene Counts data
Mutation.Frequency.Gene <- Count.All.Gene
colnames(Mutation.Frequency.Gene)[3]  <- "Freq.All"
Mutation.Frequency.Gene$Freq.Missense <- Count.Missense.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Missense.Gene))]
Mutation.Frequency.Gene$Freq.Silent   <- Count.Silent.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Silent.Gene))]
Mutation.Frequency.Gene$Freq.Nonsense <- Count.Nonsense.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Nonsense.Gene))]
Mutation.Frequency.Gene$Freq.Other    <- Count.Other.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Other.Gene))]
Mutation.Frequency.Gene$Freq.NonSilent<- Count.NonSilent.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.NonSilent.Gene))]

Mutation.Frequency.Gene$Freq.Any          <- Count.Any.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Any.Gene))]
Mutation.Frequency.Gene$Freq.Missense.Any <- Count.Missense.Any.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Missense.Any.Gene))]
Mutation.Frequency.Gene$Freq.Silent.Any   <- Count.Silent.Any.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.Silent.Any.Gene))]
Mutation.Frequency.Gene$Freq.NonSilent.Any<- Count.NonSilent.Any.Gene$freq [match(rownames(Mutation.Frequency.Gene),rownames(Count.NonSilent.Any.Gene))]

# Count the number of samples in the mutation table per cluster (N)
muts.uniquesamples <- Mutation.All[which(!duplicated(Mutation.All$Patient_ID)),c("Patient_ID","Cluster") ] 
sample.cluster.count <- as.data.frame(table(muts.uniquesamples$Cluster))
colnames (sample.cluster.count) <- c("Cluster","N")

#count(unique(Mutation.All[,c("Patient_ID","Cluster")]),vars ="Cluster")

# Gene mutation freq in percent
Mutation.Frequency.Gene$N <- sample.cluster.count$N [match(Mutation.Frequency.Gene$Cluster,sample.cluster.count$Cluster)]
Mutation.Frequency.Gene$Freq.Any.pct <- round((Mutation.Frequency.Gene$Freq.Any / Mutation.Frequency.Gene$N)*100,1)
Mutation.Frequency.Gene$Freq.Missense.Any.pct <- round((Mutation.Frequency.Gene$Freq.Missense.Any / Mutation.Frequency.Gene$N)*100,1)
Mutation.Frequency.Gene$Freq.NonSilent.Any.pct <- round((Mutation.Frequency.Gene$Freq.NonSilent.Any / Mutation.Frequency.Gene$N)*100,1)

#### ByPAtient ####

#Patient mutation frequency by Cluster table
count.Patient.bycluster <- function (Mutation.x){
  Count.Patient <- count(Mutation.x[,c("Patient_ID","Cluster")])
  Count.Patient <- Count.Patient [order(-as.numeric(Count.Patient$freq)),]
  rownames(Count.Patient) <- paste0(Count.Patient$Patient_ID,"_",Count.Patient$Cluster)
  return (Count.Patient)
}
Count.All.Patient <- count.Patient.bycluster (Mutation.All)           # raw count
Count.Missense.Patient <- count.Patient.bycluster (Mutation.Missense)
Count.Silent.Patient <- count.Patient.bycluster (Mutation.Silent)
Count.Nonsense.Patient <- count.Patient.bycluster (Mutation.Nonsense)
Count.Other.Patient <- count.Patient.bycluster (Mutation.Other)
Count.NonSilent.Patient <- count.Patient.bycluster (Mutation.NonSilent)
Count.Any.Patient <- count.Patient.bycluster (Mutation.Any)           # unique variants classification
Count.Missense.Any.Patient <- count.Patient.bycluster (Mutation.Missense.Any)
Count.Silent.Any.Patient <- count.Patient.bycluster (Mutation.Silent.Any)
Count.NonSilent.Any.Patient <- count.Patient.bycluster (Mutation.NonSilent.Any) 

#merge the byPatient Count data
Mutation.Frequency.Patient <- Count.All.Patient
colnames(Mutation.Frequency.Patient)[3]   <- "Freq.All"
Mutation.Frequency.Patient$Freq.Missense  <- Count.Missense.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Missense.Patient))]
Mutation.Frequency.Patient$Freq.Nonsense  <- Count.Nonsense.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Nonsense.Patient))]
Mutation.Frequency.Patient$Freq.Silent    <- Count.Silent.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Silent.Patient))]
Mutation.Frequency.Patient$Freq.Other     <- Count.Other.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Other.Patient))]
Mutation.Frequency.Patient$Freq.NonSilent <- Count.NonSilent.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.NonSilent.Patient))]

Mutation.Frequency.Patient$Freq.Any           <- Count.Any.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Any.Patient))]
Mutation.Frequency.Patient$Freq.Missense.Any  <- Count.Missense.Any.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Missense.Any.Patient))]
Mutation.Frequency.Patient$Freq.Silent.Any    <- Count.Silent.Any.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.Silent.Any.Patient))]
Mutation.Frequency.Patient$Freq.NonSilent.Any <- Count.NonSilent.Any.Patient$freq [match(rownames(Mutation.Frequency.Patient),rownames(Count.NonSilent.Any.Patient))]


### Save as R object or CVS ###
write.csv (Mutation.Frequency.Gene,file=paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Cancerset,"/Mutations.TCGA.",Cancerset,".",IMS.filter,".",Geneset,".Gene.by.Cluster.csv"))
write.csv (Mutation.Frequency.Patient,file=paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Cancerset,"/Mutations.TCGA.",Cancerset,".",IMS.filter,".",Geneset,".Patient.by.Cluster.csv"))
save (Mutation.Frequency.Gene,Mutation.Frequency.Patient,sample.cluster.count,file=paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",IMS.filter,".",Geneset,".Frequencies.RDATA"))
print(paste0("MAF file for ",Cancerset," Processed..."))

}

