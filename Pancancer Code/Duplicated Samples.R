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
White.list <- read.csv("./2 DATA/TCGA Whitelists/TCGA_QC_ by_Panimmune_AWG/Whitelist.20.03.17.csv")
QC.table <- read.csv("./2 DATA/TCGA Whitelists/TCGA_QC_ by_Panimmune_AWG/Sample.Quality.Annotation.TCGA.PC.Synapse.txt",sep = "\t",stringsAsFactors = FALSE)
clinical.data <- read.csv("./2 DATA/Clinical Information/PANCANCER/clinical_PANCAN_patient_with_followup.tsv",sep = "\t")
levels(clinical.data$history_of_neoadjuvant_treatment)
Neoadjuvant.Patients <- as.character(clinical.data$bcr_patient_barcode[which(clinical.data$history_of_neoadjuvant_treatment != "No")])

#test
Cancerset = "ALL"

# DO ALL
TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
if (Cancersets == "ALL") { 
  Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
}
N.sets = length(Cancersets)
Cleanup.table <- data.frame(Cancerset = Cancersets,
                            Available.samples = NA,
                            Blacklisted.samples_removed = NA,
                            QC.table.samples_removed = NA,
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
  non.whitelisted.samples <- rownames(RNASeq.subset)[-which(substring(rownames(RNASeq.subset),1,12)%in% White.list$patient_barcode)]
  if (length(non.whitelisted.samples) > 0) {RNASeq.subset <- RNASeq.subset[-which(rownames(RNASeq.subset)%in% non.whitelisted.samples),]}
  if (length(non.whitelisted.samples) > 0){print(paste0("WARNING ",length(non.whitelisted.samples)," BLACK LISTED SAMPLES REMOVED"))}
  #remove QC do not use of AWG pathology samples
  QC.table[is.na(QC.table$AWG_excluded_because_of_pathology),"AWG_excluded_because_of_pathology"] <- "x"
  non.QC.samples <- QC.table[(QC.table$AWG_excluded_because_of_pathology=="1" | QC.table$Do_not_use == "TRUE")& QC.table$cancer.type==Cancerset,"aliquot_barcode",]
  non.QC.samples <- non.QC.samples[non.QC.samples %in% rownames(RNASeq.subset)]
  if (length(non.QC.samples) > 0) {RNASeq.subset <- RNASeq.subset[-which(rownames(RNASeq.subset)%in% non.QC.samples),]}
  if (length(non.QC.samples) > 0){print(paste0("WARNING ",length(non.QC.samples)," QC SAMPLES REMOVED"))}
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
              length(non.QC.samples),
              length(normal.control.samples),
              length(metastatic.samples),
              length(recurence.samples),
              length(rownames(RNASeq.subset)[substring(rownames(RNASeq.subset),1,12)%in% doubles]),
              length(primary.samples),
              non.primary.samples,
              length(Neoadjuvant.samples))
  Cleanup.table[Cancerset,] <- Output
  dir.create(paste0("./2 DATA/SUBSETS/",DL.Method,".CLEAN.v2.1/"), showWarnings = FALSE)
  dir.create(paste0("./2 DATA/SUBSETS/",DL.Method,".CLEAN.v2.1/",Cancerset,"/"), showWarnings = FALSE)
  save (RNASeq.subset , file=paste0("./2 DATA/SUBSETS/",DL.Method,".CLEAN.v2.1/",Cancerset,"/TCGA.",Cancerset,".RNASeq.",sample.types,".subset.",Geneset,".RData"))
}
write.csv(Cleanup.table,file = "./PANCANCER/Cleanup.table.v2.1.pancancer.matrix.csv")


#cleanup leukocyte data
setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA pancancer-germline/")
infil.data <- read.csv("./Data from AWG/Immune Response AWG/TCGA_all_leuk_estimate.masked.tsv",sep = "\t",col.names = c("Cancer.Type","SampleID","Leuk.Estimate"),stringsAsFactors = FALSE)
leuk.doubles <- substring(infil.data$SampleID,1,12)[which(duplicated(substring(infil.data$SampleID,1,12)))]
leuk.doubles.data <- infil.data[substring(infil.data$SampleID,1,12)%in% leuk.doubles,]
#337
White.list <- read.csv("./Data/TCGA_QC_ by_Panimmune_AWG/Whitelist.20.03.17.csv",stringsAsFactors = FALSE)
non.whitelisted.samples <- infil.data$SampleID[-which(substring(infil.data$SampleID,1,12) %in% White.list$patient_barcode)]
infil.data <- infil.data [-which(infil.data$SampleID %in% non.whitelisted.samples),]
QC.table <- read.csv("./Data/TCGA_QC_ by_Panimmune_AWG/Sample.Quality.Annotation.TCGA.PC.Synapse.txt",sep = "\t",stringsAsFactors = FALSE)
non.QC.samples <- QC.table[QC.table$AWG_excluded_because_of_pathology=="1" | QC.table$Do_not_use == "TRUE","aliquot_barcode",]
non.QC.samples <- non.QC.samples[non.QC.samples %in% leuk.doubles.data$SampleID]
infil.data <- infil.data [-which(infil.data$SampleID %in% non.QC.samples),]
leuk.doubles <- substring(infil.data$SampleID,1,12)[which(duplicated(substring(infil.data$SampleID,1,12)))]
leuk.doubles.data <- infil.data[substring(infil.data$SampleID,1,12)%in% leuk.doubles,]
#268
metastatic.sampels <- leuk.doubles.data$SampleID[which(substring (leuk.doubles.data$SampleID,14,15 ) == "06" & leuk.doubles.data$Cancer.Type != "SKCM")] #metastatic
recurence.sampels <- leuk.doubles.data$SampleID[which(substring (leuk.doubles.data$SampleID,14,15) %in% c("02","05"))] #recurence
infil.data <- infil.data [-which(infil.data$SampleID %in% c(metastatic.sampels,recurence.sampels)),]
leuk.doubles <- substring(infil.data$SampleID,1,12)[which(duplicated(substring(infil.data$SampleID,1,12)))]
leuk.doubles.data <- infil.data[substring(infil.data$SampleID,1,12)%in% leuk.doubles,]
#92
manual.drop.list <- c("TCGA-B6-A1KC-01B-11D-A161-05","TCGA-06-0137-01B-02D-0218-05","TCGA-ER-A19T-06A-11D-A19D-05","TCGA-ER-A2NF-06A-11D-A19D-05","TCGA-D3-A1QA-07A-11D-A373-05",
                      "TCGA-23-1023-01R-01D-0807-05","TCGA-21-1076-01A-02D-0689-05","TCGA-B2-3924-01A-02D-A27A-05","TCGA-AK-3453-01A-02D-1275-05","TCGA-AK-3440-01A-02D-1275-05",
                      "TCGA-06-0137-01A-02D-0218-05","TCGA-06-0137-01A-03D-0218-05","TCGA-06-0137-01B-02D-0218-05","TCGA-06-0145-01A-02D-0218-05","TCGA-06-0145-01A-03D-0218-05",
                      "TCGA-06-0145-01A-04D-0218-05","TCGA-06-0145-01A-05D-0218-05","TCGA-06-0145-01A-06D-0218-05")
infil.data <- infil.data [-which(infil.data$SampleID %in% manual.drop.list),]
leuk.doubles <- substring(infil.data$SampleID,1,12)[which(duplicated(substring(infil.data$SampleID,1,12)))]
leuk.doubles.data <- infil.data[substring(infil.data$SampleID,1,12)%in% leuk.doubles,]
#64
leuk.doubles.data <- leuk.doubles.data[order(leuk.doubles.data$SampleID),]
double.patients <- unique(substring(leuk.doubles.data$SampleID,1,12))
for (i in 1:length(double.patients)) {
  samples.by.patient <- leuk.doubles.data$SampleID[substring(leuk.doubles.data$SampleID,1,12) == double.patients[i]]
  infil.data <- infil.data [-which(infil.data$SampleID == samples.by.patient[2]),]
}
leuk.doubles <- substring(infil.data$SampleID,1,12)[which(duplicated(substring(infil.data$SampleID,1,12)))]
leuk.doubles.data <- infil.data[substring(infil.data$SampleID,1,12)%in% leuk.doubles,]
#0

write.csv(infil.data,"./Data/Leuk.infil.data.clean.csv")
