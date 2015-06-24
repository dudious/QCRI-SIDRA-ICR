#################################################################
###
### This script compatres the scaled Immunoscore vs TP53 mutation 
###
###
#################################################################


## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

## Load Data
immunoscore <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.SKCM.ISGS.csv")
rownames(immunoscore) <- immunoscore$X
immunoscore$X <- NULL
load ("./3 ANALISYS/Mutations/SKCM/Mutation.Data.split.RDATA")

# add immunoscore to mutation data
Mutation.All$scaled.IS <- immunoscore$scaled.IS[match(Mutation.All$Patient_ID,rownames(immunoscore))]
row.names(Mutation.All) <- NULL
Mutation.Any$scaled.IS <- immunoscore$scaled.IS[match(Mutation.Any$Patient_ID,rownames(immunoscore))]
row.names(Mutation.Any) <- NULL
Mutation.Missense$scaled.IS <- immunoscore$scaled.IS[match(Mutation.Missense$Patient_ID,rownames(immunoscore))]
row.names(Mutation.Missense) <- NULL
Mutation.Nonsense$scaled.IS <- immunoscore$scaled.IS[match(Mutation.Nonsense$Patient_ID,rownames(immunoscore))]
row.names(Mutation.Nonsense) <- NULL
Mutation.Silent$scaled.IS <- immunoscore$scaled.IS[match(Mutation.Silent$Patient_ID,rownames(immunoscore))]
row.names(Mutation.Silent) <- NULL
Mutation.Other$scaled.IS <- immunoscore$scaled.IS[match(Mutation.Other$Patient_ID,rownames(immunoscore))]
row.names(Mutation.Other) <- NULL


# Select TP53 mutations
Mutation.All.TP53 <- Mutation.All[Mutation.All$Hugo_Symbol == "TP53",]
dim (Mutation.All.TP53) #64 TP53 mutations
Mutation.All.TP53 <- unique (Mutation.All.TP53) #64 mutations

Mutation.Any.TP53 <- Mutation.Any[Mutation.Any$Hugo_Symbol == "TP53",]
dim (Mutation.Any.TP53) #64 TP53 mutations
Mutation.Any.TP53 <- unique (Mutation.Any.TP53) #64 mutations

Mutation.Missense.TP53 <- Mutation.Missense[Mutation.Missense$Hugo_Symbol == "TP53",]
dim (Mutation.Missense.TP53) #64 TP53 mutations
Mutation.Missense.TP53 <- unique (Mutation.Missense.TP53) #64 mutations

Mutation.Nonsense.TP53 <- Mutation.Nonsense[Mutation.Nonsense$Hugo_Symbol == "TP53",]
dim (Mutation.Nonsense.TP53) #64 TP53 mutations
Mutation.Nonsense.TP53 <- unique (Mutation.Nonsense.TP53) #64 mutations

Mutation.Silent.TP53 <- Mutation.Silent[Mutation.Silent$Hugo_Symbol == "TP53",]
dim (Mutation.Silent.TP53) #64 TP53 mutations
Mutation.Silent.TP53 <- unique (Mutation.Silent.TP53) #64 mutations

Mutation.Other.TP53 <- Mutation.Other[Mutation.Other$Hugo_Symbol == "TP53",]
dim (Mutation.Other.TP53) #64 TP53 mutations
Mutation.Other.TP53 <- unique (Mutation.Other.TP53) #64 mutations

TP53mutated.patients <- data.frame (Patient_ID=Mutation.All.TP53$Patient_ID,P53_mutation = "TRUE")
TP53mutated.patients <- merge (TP53mutated.patients,data.frame (Patient_ID=Mutation.Missense.TP53$Patient_ID,P53_missense = "TRUE"),by="Patient_ID",all=TRUE)
TP53mutated.patients <- merge (TP53mutated.patients,unique(Mutation.All[,"Patient_ID",drop=FALSE]),by="Patient_ID",all=TRUE)
TP53mutated.patients <- unique (TP53mutated.patients)
rownames (TP53mutated.patients) <- NULL
levels(TP53mutated.patients$P53_mutation) <- c("TRUE","FALSE")
levels(TP53mutated.patients$P53_missense) <- c("TRUE","FALSE")
TP53mutated.patients[is.na(TP53mutated.patients)] <- "FALSE"
write.csv (TP53mutated.patients,file = ("./3 ANALISYS/Mutations/SKCM/TP53mutations.bypatient.csv"))

#TP53 imunoscore by cluster and by type of mutation

TP53.IS.All <- aggregate (Mutation.All.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.All.TP53$Cluster),FUN=mean)
colnames (TP53.IS.All) <- c("Cluster","TP53.All")
TP53.IS.Any <- aggregate (Mutation.Any.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Any.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Any) <- c("Cluster","TP53.Any")
TP53.IS.Missense <- aggregate (Mutation.Missense.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Missense.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Missense) <- c("Cluster","TP53.Missense")
TP53.IS.Nonsense <- aggregate (Mutation.Nonsense.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Nonsense.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Nonsense) <- c("Cluster","TP53.Nonsense")
TP53.IS.Silent <- aggregate (Mutation.Silent.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Silent.TP53$Cluster),FUN=mean) 
colnames (TP53.IS.Silent) <- c("Cluster","TP53.Silent")
TP53.IS.Other <- aggregate (Mutation.Other.TP53[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Other.TP53$Cluster),FUN=mean)
colnames (TP53.IS.Other) <- c("Cluster","TP53.Other")

TP53.IS <- cbind (TP53.IS.All,TP53.IS.Any,TP53.IS.Missense,TP53.IS.Nonsense,TP53.IS.Other) #drop Nonsense Silent  as not all clusters are available
TP53.IS <- TP53.IS[,-c(3,5,7,9)]

print (TP53.IS)
