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
  required.packages <- c("xlsx","Hmisc","HGNChelper")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  library (xlsx) #xlsx needs java installed
  library (Hmisc)
  setwd("~/Dropbox/BREAST_QATAR/")
 
  
# Load data files 
  ClinicalData.subset <- read.csv ("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.RNASeq_subset_clinicaldata.csv")                       # Clinical data including IMS
  rownames(ClinicalData.subset) <- ClinicalData.subset$X 
  ClinicalData.subset$X <-NULL
  load ("./2 DATA/SUBSETS/RNASeq.subset.ISGS.RData")
  CC.INES.RNASeq <- read.csv ("~/Dropbox/BREAST_QATAR/3 ANALISYS/CLUSTERING/RNAseq/INES/ClusterAssigmentsk4RNASeq.ICR.csv")     # Cluster assignment
  colnames(CC.INES.RNASeq) <- c("PatientID","Cluster.RNSeq.INES")
  rownames(CC.INES.RNASeq) <- CC.INES.RNASeq$PatientID
  CC.INES.RNASeq$PatientID <-NULL
  immunoscore <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.ISGS.csv")                                                     # Immunoscore
  rownames(immunoscore) <- immunoscore$X
  immunoscore$X <- NULL
  colnames(immunoscore) <- c("IS")
  mutation.freq <- read.csv ("./3 ANALISYS/Mutations/Mutations.TCGA.BRCA.Patient.by.Cluster.csv")                               # mutation frequencies
  rownames(mutation.freq) <- mutation.freq$Patient_ID
  mutation.freq <-  mutation.freq[,c("Freq.all","Freq.Missense")]
  TP53.patients <- read.csv ("./3 ANALISYS/Mutations/TP53.patients.csv")                                                        # TP53 
  rownames(TP53.patients) <- TP53.patients$Patient_ID
  TP53.patients$Patient_ID <- NULL
  TP53.patients$X <- NULL

# merge data
  Master.file <- merge (ClinicalData.subset,CC.INES.RNASeq,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,RNASeq.subset,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,immunoscore,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,mutation.freq,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,TP53.patients,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL

# export data to txt and excel
write.csv (Master.file, file = "./3 ANALISYS/MASTER FILES/TCGA.BRCA.RNASeq_subset_ISGS.Master.csv",row.names = TRUE);
write.xlsx (Master.file, file = "./3 ANALISYS/MASTER FILES/TCGA.BRCA.RNASeq_subset_ISGS.Master.xlsx", sheetName ="RNASeq ISGS.Master data", row.names=TRUE)