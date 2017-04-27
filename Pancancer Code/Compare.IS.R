#################################################################
###
### This script calculated the correlation between Immunoscores 
###
### 
###
#################################################################

## Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR/")

# Set Parameters
Cancersets        = "ALL"
Geneset           = "DBGS3"   
DL.Method         = "PANCANCER.CLEAN.v2"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Split"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer



# TEST : "Cancerset = "BLCA"
# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}

clinical.data <- read.csv("./2 DATA/Clinical Information/PANCANCER/clinical_PANCAN_patient_with_followup.tsv",sep = "")

N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}

#load Data
  
IS.BIOLINKS <- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BIOLINKS.",Cancerset,".DBGS3.csv"))
rownames(IS.BIOLINKS) <- IS.BIOLINKS$X
IS.BIOLINKS$X <- NULL
IS.PANCANCER<- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.CLEAN.v2.",Cancerset,".DBGS3.csv"))
rownames(IS.PANCANCER) <- gsub("\\.","-",make.names(IS.PANCANCER$X,unique = TRUE)) #DUPLICATES IN PANCANCER MATRIX
IS.PANCANCER$X <- NULL

print (paste0(length(which(substring(rownames(IS.PANCANCER),1,12) %in% rownames(IS.BIOLINKS)))," Overlapping ",Cancerset," samples"))
IS.PANCANCER <- IS.PANCANCER[rownames(IS.BIOLINKS),]

print (paste0 ("Mean unscaled IS using BIOLINKS Data is ",round (mean(IS.BIOLINKS$unscaled.IS,NA.rm=TRUE),1),
               " While using the PANCANCER matrix it is ", round (mean(IS.PANCANCER$unscaled.IS,na.rm=TRUE),1)
))
scaled.cor <- cor.test(IS.BIOLINKS$scaled.IS,IS.PANCANCER$scaled.IS)
unscaled.cor <- cor.test(IS.BIOLINKS$unscaled.IS,IS.PANCANCER$unscaled.IS)
print (paste0( "with correlation of ", round(unscaled.cor$estimate,2)))
}

#IS Master file
Immunescore.master <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.UNSPLIT.DBGS3.csv",row.names = 1)
Immunescore.master$Patient_ID <- substring(Immunescore.master$Patient_ID,1,12)

N.sets = length(Cancersets)
Immunescore.master$Scaled.By.Cancer <- "NA"
Immunescore.master$Scaled.By.Cancer.ACT <- "NA"
Immunescore.master$Scaled.By.Cancer.INH <- "NA"
Immunescore.master$Cancer.Type <- "NA"
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}
  Scaled.by.cancer <- read.csv(paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.CLEAN.v2.",Cancerset,".DBGS3.csv"),row.names = 1)
  new.data <- Scaled.by.cancer$scaled.IS[na.omit(match(rownames(Immunescore.master[rownames(Scaled.by.cancer),]),rownames(Scaled.by.cancer)))]
  Immunescore.master[rownames(Scaled.by.cancer),"Scaled.By.Cancer"] <- new.data
  new.data <- Scaled.by.cancer$scaled.IS.ACT[na.omit(match(rownames(Immunescore.master[rownames(Scaled.by.cancer),]),rownames(Scaled.by.cancer)))]
  Immunescore.master[rownames(Scaled.by.cancer),"Scaled.By.Cancer.ACT"] <- new.data
  new.data <- Scaled.by.cancer$scaled.IS.INH[na.omit(match(rownames(Immunescore.master[rownames(Scaled.by.cancer),]),rownames(Scaled.by.cancer)))]
  Immunescore.master[rownames(Scaled.by.cancer),"Scaled.By.Cancer.INH"] <- new.data
}
Immunescore.master$Scaled.By.Cancer <- as.numeric(Immunescore.master$Scaled.By.Cancer)
Immunescore.master$Scaled.By.Cancer.ACT <- as.numeric(Immunescore.master$Scaled.By.Cancer.ACT)
Immunescore.master$Scaled.By.Cancer.INH <- as.numeric(Immunescore.master$Scaled.By.Cancer.INH)
Immunescore.master$Cancer.Type <- clinical.data$acronym[match(Immunescore.master$Patient_ID,clinical.data$bcr_patient_barcode)]
Immunescore.master.clean <- Immunescore.master[complete.cases(Immunescore.master),]
write.csv(Immunescore.master.clean,file = "./3 ANALISYS/IMMUNOSCORE/immunoscore.Master.clean.v2.csv")
cor.test(Immunescore.master.clean$ALL,Immunescore.master.clean$Scaled.By.Cancer)
Immunescore.master.clean$Scaled.By.Cancer
summary(Immunescore.master.clean$Cancer.Type)
which(duplicated(Immunescore.master.clean$Patient_ID))

