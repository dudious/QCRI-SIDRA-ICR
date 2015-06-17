#################################################################
###
### This Script combines the clinical data relating to patients 
### that have MA Data available with IMS predictions based 
### based on both MA and MicroArray data.
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
  ClinicalData.subset <- read.csv ("./3 ANALISYS/CLINICAL DATA/TCGA.GBM.MA.Affy_subset_clinicaldata.csv")                       # Clinical data including IMS
  rownames(ClinicalData.subset) <- ClinicalData.subset$X 
  ClinicalData.subset$X <-NULL
  load ("./2 DATA/SUBSETS/GBM/MA.subset.ISGS.RData")
  CC.MA <- read.csv ("~/Dropbox/BREAST_QATAR/3 ANALISYS/CLUSTERING/MA/GBM/GBM.TCGA.MA.Affy.k7.ISGS.reps5000/GBM.TCGA.MA.Affy.k7.ISGS.reps5000.k=4.consensusClass.ICR.csv")     # Cluster assignment
  rownames(CC.MA) <- CC.MA$X 
  CC.MA$X <- NULL
  colnames(CC.MA) <- c("PatientID","Cluster.ISGS.RNSeq")
  CC.MA$PatientID <-NULL
  immunoscore <- read.csv ("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.MA.Affy.GBM.ISGS.csv")                                                     # Immunoscore
  rownames(immunoscore) <- immunoscore$X
  immunoscore$X <- NULL
  mutation.freq <- read.csv ("./3 ANALISYS/Mutations/GBM/Mutations.TCGA.GBM.MA.Affy.Gene.by.Cluster.csv")                               # mutation frequencies
  rownames(mutation.freq) <- mutation.freq$Patient_ID
  mutation.freq <-  mutation.freq[,c("Freq.All","Freq.Missense")]
  TP53.patients <- read.csv ("./3 ANALISYS/Mutations/GBM/TP53mutations.MA.Affy.bypatient.csv")                                                        # TP53 
  rownames(TP53.patients) <- TP53.patients$Patient_ID
  TP53.patients$Patient_ID <- NULL
  TP53.patients$X <- NULL

# merge data
  Master.file <- merge (ClinicalData.subset,CC.MA,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,MA.subset,by="row.names",all.x=TRUE)
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
write.csv (Master.file, file = "./3 ANALISYS/MASTER FILES/TCGA.GBM.MA.Affy_subset_ISGS.Master.csv",row.names = TRUE);
write.xlsx (Master.file, file = "./3 ANALISYS/MASTER FILES/TCGA.GBM.MA.Affy_subset_ISGS.Master.xlsx", sheetName ="MA ISGS.Master data", row.names=TRUE)
