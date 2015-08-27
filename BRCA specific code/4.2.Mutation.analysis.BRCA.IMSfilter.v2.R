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
setwd("~/Dropbox/BREAST_QATAR/")
 #Dependencies
 required.packages <- c("plyr")
 missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
 if(length(missing.packages)) install.packages(missing.packages)
library("plyr")

# Set Parameters
Cancerset <- "BRCA"           # do not use -GA or -hiseq (data is merged)
BRCA.Filter <- "BSF2"          # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
K <- 4                        # SET K here
IMS.filter = "Her2"           # Alterantives "All" , "Luminal" , "Basal", "Her2","LumA","LumB"

## Load Data
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
 Clinical.data <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"),header=TRUE)
 rownames(Clinical.data) <- Clinical.data[,1]
 Clinical.data[,1] <-NULL
 rm(maf.merged.table)
}

dir.create (paste0("./3 ANALISYS/Mutations/",Cancerset,"/"),showWarnings=FALSE)

# Add Class to mutation data
Mutation.selected.data$Cluster <- Consensus.class$Cluster [match(Mutation.selected.data$Patient_ID,Consensus.class$Patient_ID)]
Mutation.selected.data$Subtype <- Clinical.data$TCGA.PAM50.RMethod.RNASeq [match(Mutation.selected.data$Patient_ID,rownames(Clinical.data))]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Cluster)),]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Subtype)),]


#Filer  Mutation.selected.data Table by IMS
if (IMS.filter == "Luminal") {
  Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Subtype %in% c("Luminal A","Luminal B"),]
} else if (IMS.filter == "Basal"){
  Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Subtype %in% c("Basal-like"),]
} else if (IMS.filter == "Her2"){
  Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Subtype %in% c("HER2-enriched"),]
} else if (IMS.filter == "LumA"){
  Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Subtype %in% c("Luminal A"),]
} else if (IMS.filter == "LumB"){
  Mutation.selected.data <- Mutation.selected.data[Mutation.selected.data$Subtype %in% c("Luminal B"),]
}

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
      file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",IMS.filter,".",Geneset,".split.RDATA"))

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

# Gene mutation freq in percent
Mutation.Frequency.Gene$N <- sample.cluster.count$N [match(Mutation.Frequency.Gene$Cluster,sample.cluster.count$Cluster)]
Mutation.Frequency.Gene$Freq.Any.pct <- floor((Mutation.Frequency.Gene$Freq.Any / Mutation.Frequency.Gene$N)*100)
Mutation.Frequency.Gene$Freq.Missense.Any.pct <- floor((Mutation.Frequency.Gene$Freq.Missense.Any / Mutation.Frequency.Gene$N)*100)
Mutation.Frequency.Gene$Freq.NonSilent.Any.pct <- floor((Mutation.Frequency.Gene$Freq.NonSilent.Any / Mutation.Frequency.Gene$N)*100)

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
write.csv (Mutation.Frequency.Gene,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutations.TCGA.",Cancerset,".",IMS.filter ,".",Geneset,".Gene.by.Cluster.csv"))
write.csv (Mutation.Frequency.Patient,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutations.TCGA.",Cancerset,".",IMS.filter ,".",Geneset,".Patient.by.Cluster.csv"))
save (Mutation.Frequency.Gene,Mutation.Frequency.Patient,sample.cluster.count,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",IMS.filter,".",Geneset,".Frequencies.RDATA"))




