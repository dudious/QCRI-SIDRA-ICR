#################################################################################
###
### This script compatres the scaled Immunoscore vs Gene of FOCUS (GOF) mutation 
###
###
################################################################################


## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

# Set Parameters
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE
Parent.Geneset <- substring(Geneset,1,5)
Cancerset <- "BRCA.BSF"
Parent.Cancerset <- substring(Cancerset,1,4)
K <- 4
GOF <- "TP53"

## Load Data
immunoscore <- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.",Parent.Cancerset,".",Parent.Geneset,".csv"))
rownames(immunoscore) <- immunoscore$X
immunoscore$X <- NULL
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",Geneset,".split.RDATA"))

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


# Select GOF mutations
Mutation.All.GOF <- Mutation.All[Mutation.All$Hugo_Symbol == GOF,]
dim (Mutation.All.GOF) #64 GOF mutations
Mutation.All.GOF <- unique (Mutation.All.GOF) #64 mutations

Mutation.Any.GOF <- Mutation.Any[Mutation.Any$Hugo_Symbol == GOF,]
dim (Mutation.Any.GOF) #64 GOF mutations
Mutation.Any.GOF <- unique (Mutation.Any.GOF) #64 mutations

Mutation.Missense.GOF <- Mutation.Missense[Mutation.Missense$Hugo_Symbol == GOF,]
dim (Mutation.Missense.GOF) #64 GOF mutations
Mutation.Missense.GOF <- unique (Mutation.Missense.GOF) #64 mutations

Mutation.Nonsense.GOF <- Mutation.Nonsense[Mutation.Nonsense$Hugo_Symbol == GOF,]
dim (Mutation.Nonsense.GOF) #64 GOF mutations
Mutation.Nonsense.GOF <- unique (Mutation.Nonsense.GOF) #64 mutations

Mutation.Silent.GOF <- Mutation.Silent[Mutation.Silent$Hugo_Symbol == GOF,]
dim (Mutation.Silent.GOF) #64 GOF mutations
Mutation.Silent.GOF <- unique (Mutation.Silent.GOF) #64 mutations

Mutation.Other.GOF <- Mutation.Other[Mutation.Other$Hugo_Symbol == GOF,]
dim (Mutation.Other.GOF) #64 GOF mutations
Mutation.Other.GOF <- unique (Mutation.Other.GOF) #64 mutations

GOFmutated.patients <- data.frame (Patient_ID=Mutation.All.GOF$Patient_ID,GOF_mutation = "TRUE")
GOFmutated.patients <- merge (GOFmutated.patients,data.frame (Patient_ID=Mutation.Missense.GOF$Patient_ID,GOF_missense = "TRUE"),by="Patient_ID",all=TRUE)
GOFmutated.patients <- merge (GOFmutated.patients,unique(Mutation.All[,"Patient_ID",drop=FALSE]),by="Patient_ID",all=TRUE)
GOFmutated.patients <- unique (GOFmutated.patients)
rownames (GOFmutated.patients) <- NULL
levels(GOFmutated.patients$GOF_mutation) <- c("TRUE","FALSE")
levels(GOFmutated.patients$GOF_missense) <- c("TRUE","FALSE")
GOFmutated.patients[is.na(GOFmutated.patients)] <- "FALSE"
dir.create(paste0("./3 ANALISYS/Mutations/",Cancerset,"/"), showWarnings = FALSE)
write.csv (GOFmutated.patients,file = paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",GOF,".mutations.",Geneset,".bypatient.csv"))

#GOF imunoscore by cluster and by type of mutation

GOF.IS.All <- aggregate (Mutation.All.GOF[,c("scaled.IS"),drop=FALSE],by=list(Mutation.All.GOF$Cluster),FUN=mean)
colnames (GOF.IS.All) <- c("Cluster","GOF.All")
GOF.IS.Any <- aggregate (Mutation.Any.GOF[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Any.GOF$Cluster),FUN=mean)
colnames (GOF.IS.Any) <- c("Cluster","GOF.Any")
GOF.IS.Missense <- aggregate (Mutation.Missense.GOF[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Missense.GOF$Cluster),FUN=mean)
colnames (GOF.IS.Missense) <- c("Cluster","GOF.Missense")
GOF.IS.Nonsense <- aggregate (Mutation.Nonsense.GOF[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Nonsense.GOF$Cluster),FUN=mean)
colnames (GOF.IS.Nonsense) <- c("Cluster","GOF.Nonsense")
GOF.IS.Silent <- aggregate (Mutation.Silent.GOF[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Silent.GOF$Cluster),FUN=mean) 
colnames (GOF.IS.Silent) <- c("Cluster","GOF.Silent")
GOF.IS.Other <- aggregate (Mutation.Other.GOF[,c("scaled.IS"),drop=FALSE],by=list(Mutation.Other.GOF$Cluster),FUN=mean)
colnames (GOF.IS.Other) <- c("Cluster","GOF.Other")

GOF.IS <- cbind (GOF.IS.All,GOF.IS.Any,GOF.IS.Missense,GOF.IS.Nonsense,GOF.IS.Other) #drop Nonsense or Silent or other if not all clusters are available
GOF.IS <- GOF.IS[,-c(3,5,7,9)]

print (GOF.IS)
save(GOF.IS,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",GOF,".IS.",Geneset,".RData"))
write.table(GOF.IS,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",GOF,".IS.",Geneset,".txt"),sep="\t")
