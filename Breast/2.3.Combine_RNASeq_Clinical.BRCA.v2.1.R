#################################################################
###
### This Script combines the clinical data relating to patients 
### that have RNASeq Data available with IMS predictions based 
### based on both RNASeq and MicroArray data.
###
#################################################################

# Setup environment
  rm(list=ls())
  ## dependencies
  ## install java for xlsx export
  ## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
  required.packages <- c("xlsx","Hmisc","HGNChelper","gridExtra")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  library (xlsx) #xlsx needs java installed
  library (Hmisc)
  setwd("~/Dropbox/BREAST_QATAR/")
  source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
  source("./1 CODE/R tools/TCGA-Assembler/Module_B.r")
  
# Load data files 
  ## RNASeq DAta from TCGA assembler
  load ("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")
  PatientIDs <- unique(substr(colnames(RNASeq.NORM),1,12)) 
  rm(RNASeq.NORM)
  ## clinical data retrieved through TCGA assembler
  clinicalData <- read.csv ("./2 DATA/Clinical Information/BRCA/selected_clinicaldata.txt", header = TRUE)
  
  ## IMS prediction data from TCGA.PAM50.METHOD using Agilent data through ASSEMBLER
  IMS.TCGA.METHOD<-read.csv("./3 ANALISYS/IMS/TCGA IMS/MA/BRCA.TCGA.MA.IMS.OUTPUT_pam50scores.txt",sep="\t",header = TRUE)
  IMS.TCGA.METHOD <- IMS.TCGA.METHOD[complete.cases(IMS.TCGA.METHOD),] ## remove incomlete cases
  rownames (IMS.TCGA.METHOD) <- IMS.TCGA.METHOD[,1]
  IMS.TCGA.METHOD[,1] <- NULL
  ### rename PAM50 names
  levels(IMS.TCGA.METHOD$Call) <- c(levels(IMS.TCGA.METHOD$Call), c("Basal-like","HER2-enriched","Luminal A","Luminal B","Normal-like")) 
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="Basal"]  <- "Basal-like"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="Her2"]   <- "HER2-enriched"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="LumA"]   <- "Luminal A"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="LumB"]   <- "Luminal B"
  IMS.TCGA.METHOD$Call[IMS.TCGA.METHOD$Call=="Normal"] <- "Normal-like"
  IMS.TCGA.METHOD$Call <- droplevels(IMS.TCGA.METHOD$Call)
  
  ## IMS  prediction data from TCGA.PAM50.METHOD using RNASeq data through ASSEMBLER and EDASeq normalized (with Log2 transformation)
  IMS.TCGA.METHOD.RNASeq<-read.csv("./3 ANALISYS/IMS/TCGA IMS/RNASeq/EDASEQ_LOG2/BRCA.TCGA.EDASeq.RNASeq.IMS.OUTPUT_pam50scores.txt",sep="\t",header = TRUE)
  IMS.TCGA.METHOD.RNASeq <- IMS.TCGA.METHOD.RNASeq[complete.cases(IMS.TCGA.METHOD.RNASeq),] ## remove incomplete cases
  rownames (IMS.TCGA.METHOD.RNASeq) <- IMS.TCGA.METHOD.RNASeq[,1]
  IMS.TCGA.METHOD.RNASeq[,1] <- NULL
  ### rename PAM50 names
  levels(IMS.TCGA.METHOD.RNASeq$Call) <- c(levels(IMS.TCGA.METHOD.RNASeq$Call), c("Basal-like","HER2-enriched","Luminal A","Luminal B","Normal-like")) 
  IMS.TCGA.METHOD.RNASeq$Call[IMS.TCGA.METHOD.RNASeq$Call=="Basal"]  <- "Basal-like"
  IMS.TCGA.METHOD.RNASeq$Call[IMS.TCGA.METHOD.RNASeq$Call=="Her2"]   <- "HER2-enriched"
  IMS.TCGA.METHOD.RNASeq$Call[IMS.TCGA.METHOD.RNASeq$Call=="LumA"]   <- "Luminal A"
  IMS.TCGA.METHOD.RNASeq$Call[IMS.TCGA.METHOD.RNASeq$Call=="LumB"]   <- "Luminal B"
  IMS.TCGA.METHOD.RNASeq$Call[IMS.TCGA.METHOD.RNASeq$Call=="Normal"] <- "Normal-like"
  IMS.TCGA.METHOD.RNASeq$Call <- droplevels(IMS.TCGA.METHOD.RNASeq$Call)
  
# subset clinical data
  ClinicalData.subset <- subset (clinicalData,is.element (clinicalData$bcr_patient_barcode,PatientIDs))
  row.names(ClinicalData.subset) <- ClinicalData.subset$bcr_patient_barcode
  ClinicalData.subset$bcr_patient_barcode <- NULL
# merge clinical data with IMS data
  ClinicalData.subset <- merge(ClinicalData.subset,IMS.TCGA.METHOD["Call"],by="row.names",all.x=TRUE, all.y=FALSE)
  row.names(ClinicalData.subset) <- ClinicalData.subset$Row.names
  ClinicalData.subset$Row.names <- NULL
  colnames(ClinicalData.subset)[colnames(ClinicalData.subset) == "Call"] <- "TCGA.PAM50.RMethod.MA"
  ClinicalData.subset <- merge(ClinicalData.subset,IMS.TCGA.METHOD.RNASeq["Call"],by="row.names",all.x=TRUE, all.y=FALSE)
  row.names(ClinicalData.subset) <- ClinicalData.subset$Row.names
  ClinicalData.subset$Row.names <- NULL
  colnames(ClinicalData.subset)[colnames(ClinicalData.subset) == "Call"] <- "TCGA.PAM50.RMethod.RNASeq" 
  print ("Data merged...")

# test match in PAM50 IMS
  ClinicalData.subset <- cbind(ClinicalData.subset,ClinicalData.subset$TCGA.PAM50.RMethod.MA == ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq  )
  col.names <- colnames(ClinicalData.subset)
  col.names[length(col.names)] <- "MA.PAM50vsRNASeq.PAM50"
  colnames(ClinicalData.subset) <- col.names
  ClinicalData.subset.complete <- ClinicalData.subset[complete.cases(ClinicalData.subset),]
  
  print ("MA.PAM50")
  print (table (ClinicalData.subset.complete$TCGA.PAM50.RMethod.MA))
  print ("RNASeq.PAM50")
  print (table (ClinicalData.subset.complete$TCGA.PAM50.RMethod.RNASeq))
  print (table (ClinicalData.subset$MA.PAM50vsRNASeq.PAM50))

  nomatch <- ClinicalData.subset[which(ClinicalData.subset$MA.PAM50vsRNASeq.PAM50 == FALSE),27:29]
  nomatch <- merge(nomatch,IMS.TCGA.METHOD.RNASeq[,1:7],by="row.names",all.x=TRUE, all.y=FALSE)
  write.csv (nomatch, file = "./3 ANALISYS/IMS/TCGA IMS/RNASeq/nomatchwithMA.csv",row.names = TRUE);

# append exclusion parameters
  exclude.samples.males <- which(ClinicalData.subset$gender == "MALE")
  exclude.samples.nat <- which(ClinicalData.subset$history_neoadjuvant_treatment == "Yes")
  exclude.samples.histo <- which(ClinicalData.subset$histological_type %nin% c("Infiltrating Ductal Carcinoma","Infiltrating Lobular Carcinoma"))
  exclude.samples.ims <- which(ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq == "Normal-like")
  exclude.samples.history <- which(ClinicalData.subset$history_other_malignancy == "Yes")
  exclude.samples.preclust <- unique(c(exclude.samples.males,exclude.samples.nat,exclude.samples.ims)) #exclude.samples.males,exclude.samples.nat,exclude.samples.histo,exclude.samples.ims,exclude.samples.history
  exclude.samples.postclust <- unique(c(exclude.samples.history))
  ClinicalData.subset$exclude.pre <- "No"
  ClinicalData.subset$exclude.post <- "No"
  ClinicalData.subset[exclude.samples.preclust,"exclude.pre"] <-"Yes"
  ClinicalData.subset[exclude.samples.postclust,"exclude.post"] <-"Yes"
  print ("exclusion parameter added...")
  
# export data to txt and excel
  write.csv (ClinicalData.subset, file = "./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.RNASeq_subset_clinicaldata.csv",row.names = TRUE);
  write.xlsx (ClinicalData.subset, file = "./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.RNASeq_subset_clinicaldata.xlsx", sheetName ="RNASeq subset clinical data", row.names=TRUE);
  print ("Data on all Samples are saved in RNASeq_subset_clinicaldata.xlsx and RNASeq_subset_clinicaldata.txt.")
 
