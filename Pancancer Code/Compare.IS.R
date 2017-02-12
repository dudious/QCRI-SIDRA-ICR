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
DL.Method         = "PANCANCER.CLEAN"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Split"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer



# TEST : "Cancerset = "BLCA"
# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}

#load Data
IS.BIOLINKS <- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.BIOLINKS.",Cancerset,".DBGS3.csv"))
rownames(IS.BIOLINKS) <- IS.BIOLINKS$X
IS.BIOLINKS$X <- NULL
IS.PANCANCER<- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.CLEAN.",Cancerset,".DBGS3.csv"))
rownames(IS.PANCANCER) <- gsub("\\.","-",make.names(IS.PANCANCER$X,unique = TRUE)) #DUPLICATES IN PANCANCER MATRIX
IS.PANCANCER$X <- NULL

print (paste0(length(which(rownames(IS.PANCANCER) %in% rownames(IS.BIOLINKS)))," Overlapping ",Cancerset," samples"))
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
Immunescore.master$Cancer.Type <- "NA"
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}
  Scaled.by.cancer <- read.csv(paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.PANCANCER.CLEAN.",Cancerset,".DBGS3.csv"),row.names = 1)
  new.data <- Scaled.by.cancer$scaled.IS[na.omit(match(rownames(Immunescore.master[rownames(Scaled.by.cancer),]),rownames(Scaled.by.cancer)))]
  Immunescore.master[rownames(Scaled.by.cancer),"Scaled.By.Cancer"] <- new.data
  Immunescore.master[rownames(Scaled.by.cancer),"Cancer.Type"] <- Cancerset
}
Immunescore.master$Scaled.By.Cancer <- as.numeric(Immunescore.master$Scaled.By.Cancer)
write.csv(Immunescore.master,file = "./3 ANALISYS/IMMUNOSCORE/immunoscore.Master.clean.csv")
cor.test(Immunescore.master$ALL,Immunescore.master$Scaled.By.Cancer)
Immunescore.master$Scaled.By.Cancer

