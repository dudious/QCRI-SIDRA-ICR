#################################################################
###
### This script analyses the duplicates in the PANCANCER RNASeq Matrix 
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
DL.Method         = "PANCANCER"     #Choose "ASSEMBLER" or "BIOLINKS" or "PANCANCER"
sample.types      = "Split"     #Alternatives TP , TP_TM , Selected or "Split" for Pancancer

#load data
White.list <- read.csv("./2 DATA/TCGA Whitelists/Data Freeze 1.3.1.csv")
clinical.data <- read.csv("./2 DATA/Clinical Information/PANCANCER/clinical_PANCAN_patient_with_followup.tsv",sep = "\t")
levels(clinical.data$history_of_neoadjuvant_treatment)
Neoadjuvant.Patients <- as.character(clinical.data$bcr_patient_barcode[which(clinical.data$history_of_neoadjuvant_treatment != "No")])

#test
Cancerset = "SKCM"

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
Cleanup.table <- data.frame(Cancerset = Cancersets,
                            Available.samples = NA,
                            Blacklisted.samples_removed = NA,
                            Normal.control.samples_removed = NA,
                            Metastatic.samples_removed = NA,
                            Recurence.samples_removed = NA,
                            Non.unique.samples_left = NA,
                            Primary.samples = NA,
                            Non.Primary.samples = NA,
                            Neoadjuvant.samples = NA)
rownames(Cleanup.table) <-Cleanup.table$Cancerset
Cleanup.table$Cancerset <- NULL
for (i in 1:N.sets) {
  Cancerset = Cancersets[i]
  if (Cancerset %in% c("LAML","FPPP")) {next}
  
  ## Load Data
  load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
  RNASeq.subset <- as.matrix(RNASeq.subset)
  available.samples <- nrow(RNASeq.subset)
  print(paste0(Cancerset," has ",available.samples," total samples available"))
  doubles <- substring(rownames(RNASeq.subset),1,12)[which(duplicated(substring(rownames(RNASeq.subset),1,12)))]
  print(paste0(Cancerset," has ",length(rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles])," non-unique samples."))
  rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles]
  #remove non-whitelisted samples
  non.whitelisted.samples <- rownames(RNASeq.subset)[-which(rownames(RNASeq.subset)%in% White.list$aliquot_barcode)]
  if (length(non.whitelisted.samples) > 0) {RNASeq.subset <- RNASeq.subset[-which(rownames(RNASeq.subset)%in% non.whitelisted.samples),]}
  if (length(non.whitelisted.samples)>0){print(paste0("WARNING ",length(non.whitelisted.samples)," BLACK LISTED SAMPLES REMOVED"))}
  #remove normal controls
  normal.control.samples <- which(substring(rownames(RNASeq.subset),14,15) == "11")
  if(length(normal.control.samples)>0){RNASeq.subset <- RNASeq.subset[-normal.control.samples,]}
  doubles <- substring(rownames(RNASeq.subset),1,12)[which(duplicated(substring(rownames(RNASeq.subset),1,12)))]
  print(paste0("Afer removal of ",length(normal.control.samples)," normal control samples ",Cancerset," has ",length(rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles])," non-unique samples."))
  #remove metastatic samples
  metastatic.samples <- NULL
  if (Cancerset != "SKCM") {
    metastatic.samples <- which(substring(rownames(RNASeq.subset),14,15) == "06")
    if(length(metastatic.samples)>0){RNASeq.subset <- RNASeq.subset[-metastatic.samples,]}
    doubles <- substring(rownames(RNASeq.subset),1,12)[which(duplicated(substring(rownames(RNASeq.subset),1,12)))]
    print(paste0("Afer removal of ",length(metastatic.samples)," metastatic samples ",Cancerset," has ",length(rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles])," non-unique samples."))
  }
  #manual.remove.SKCM <- c("TCGA-D3-A1QA-07A-11R-A37K-07","TCGA-ER-A19T-06A-11R-A18U-07","TCGA-ER-A2NF-06A-11R-A18T-07")
  #manual.remove.SKCM.remove <- which(rownames(RNASeq.subset)%in% manual.remove.SKCM)
  #if (length(manual.remove.SKCM.remove) > 0) {RNASeq.subset <- RNASeq.subset[-manual.remove.SKCM.remove,]}
  #remove recurence samples
  recurence.samples <- which(substring(rownames(RNASeq.subset),14,15) %in% c("02","05"))
  if(length(recurence.samples)>0){RNASeq.subset <- RNASeq.subset[-recurence.samples,]}
  doubles <- substring(rownames(RNASeq.subset),1,12)[which(duplicated(substring(rownames(RNASeq.subset),1,12)))]
  print(paste0("Afer removal of ",length(recurence.samples)," recurence samples ",Cancerset," has ",length(rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles])," non-unique samples."))
  
  #check left over samples
  if(Cancerset != "SKCM") {
  if(length(doubles)>0){
    duplicate.primary.samples <- rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles]
    print(paste0("Duplicate sample : " ,duplicate.primary.samples))
    #remove "B" vials
    B.vial.samples <- which(substring(rownames(RNASeq.subset),16,16) == "B" & rownames(RNASeq.subset)%in% duplicate.primary.samples)
    non02.portion.samples <- which(substring(rownames(RNASeq.subset),18,19) == "02" & rownames(RNASeq.subset)%in% duplicate.primary.samples)
    print ("B vials and non 02 portions removed for duplicates only.")
    RNASeq.subset <- RNASeq.subset[-c(B.vial.samples,non02.portion.samples),]
    doubles <- substring(rownames(RNASeq.subset),1,12)[which(duplicated(substring(rownames(RNASeq.subset),1,12)))]
  }}
  #remove samples that dont cluster
  not.clutering.samples <- c("TCGA-23-1023-01A-02R-1564-13","TCGA-23-1023-01R-01R-1564-13","TCGA-06-0211-01A-01R-1849-01","TCGA-06-0211-01B-01R-1849-01")
  not.clutering.samples.to.remove <- which(rownames(RNASeq.subset)%in% not.clutering.samples)
  if (length(not.clutering.samples.to.remove) > 0) {RNASeq.subset <- RNASeq.subset[-not.clutering.samples.to.remove,]}
  if(length(doubles)==0) {print(paste0(Cancerset," has ",nrow(RNASeq.subset)," unique samples."))}
  primary.samples <- which(substring(rownames(RNASeq.subset),14,15) == "01")
  print(paste0(length(primary.samples)," ",Cancerset," samples are from primary tumor origin"))
  non.primary.samples <- nrow(RNASeq.subset)-length(primary.samples)
  #Neoadjuvant samples
  Neoadjuvant.samples <- rownames(RNASeq.subset)[which(substring(rownames(RNASeq.subset),1,12) %in% Neoadjuvant.Patients)]
  print(paste0(length(Neoadjuvant.samples)," ",Cancerset," NON BLACKLISTED Neoadjuvant samples present"))
  print("_____________________________________________________________________")
  Output <- c(available.samples,
              length(non.whitelisted.samples),
              length(normal.control.samples),
              length(metastatic.samples),
              length(recurence.samples),
              length(rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles]),
              length(primary.samples),
              non.primary.samples,
              length(Neoadjuvant.samples))
  Cleanup.table[Cancerset,] <- Output
  dir.create(paste0("./2 DATA/SUBSETS/",DL.Method,".CLEAN/",Cancerset,"/"), showWarnings = FALSE)
  save (RNASeq.subset , file=paste0("./2 DATA/SUBSETS/",DL.Method,".CLEAN/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
}
write.csv(Cleanup.table,file = "./PANCANCER/Cleanup.table.pancancer.matrix.csv")

