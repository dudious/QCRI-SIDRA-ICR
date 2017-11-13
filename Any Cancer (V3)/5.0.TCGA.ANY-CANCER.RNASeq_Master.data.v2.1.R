#################################################################
###
### This Script Creates a Master File with all Data by Patient
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

# Set Parameters
  Cancerset <- "BRCA.BSF2"
  Parent.Cancerset <- substring(Cancerset,1,4)
  Geneset   <- "DBGS3.FLTR"       # SET GENESET HERE
  Parent.Geneset <- substring(Geneset,1,5)
  K <- 4          
  
# Load data files 
  #clinical data
  ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
  rownames(ClinicalData.subset) <- ClinicalData.subset$X 
  ClinicalData.subset$X <-NULL
  #RNASeq data
  load (paste0("./2 DATA/SUBSETS/",Parent.Cancerset,"/TCGA.",Parent.Cancerset,".RNASeq.subset.",Parent.Geneset,".RData"))
  #Consensus clustering data
  CC.RNASeq <- read.csv (paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",
                                Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",
                                Geneset,".reps5000.k=4.consensusClass.ICR.csv"))                                                # Cluster assignment
  rownames(CC.RNASeq) <- CC.RNASeq$X 
  CC.RNASeq$X <- NULL
  colnames(CC.RNASeq) <- c("PatientID",paste0("Cluster.",Geneset,".RNSeq"))
  CC.RNASeq$PatientID <-NULL
  #immunoscore data
  immunoscore <- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.",Parent.Cancerset,".",Parent.Geneset,".csv"))                            # Immunoscore
  rownames(immunoscore) <- immunoscore$X
  immunoscore$X <- NULL
  #mutation frequnecy data
  mutation.freq <- read.csv (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutations.TCGA.",Cancerset,".",Geneset,".Patient.by.Cluster.csv"))      # mutation frequencies
  rownames(mutation.freq) <- mutation.freq$Patient_ID
  mutation.freq <-  mutation.freq[,c("Freq.All","Freq.Missense")]
  #TP53.mutaion status
  TP53.patients <- read.csv (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".TP53.mutations.",Geneset,".bypatient.csv"))      # TP53 
  rownames(TP53.patients) <- TP53.patients$Patient_ID
  TP53.patients$Patient_ID <- NULL
  TP53.patients$X <- NULL
  #MAPK STATUS
  Mutstat <- read.csv ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
  

# merge data
  Master.file <- merge (ClinicalData.subset,CC.RNASeq,by="row.names",all.x=TRUE)
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
write.csv (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.csv"),row.names = TRUE);
write.xlsx (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.xlsx"), sheetName ="RNASeq ISGS.Master data", row.names=TRUE)

#multivariate analysis
Master.file$Cluster.DBGS3.FLTR.RNSeq
MV.data <- Master.file[,c("TCGA.PAM50.RMethod.RNASeq","GOF_mutation","Cluster.DBGS3.FLTR.RNSeq")]
MV.data$MAPKMUT <- Mutstat$MAP2K4.MAP3K1[match(rownames(MV.data),Mutstat$Sample)]
MV.data <- MV.data[complete.cases(MV.data),]
MV.data$Basal <- "FALSE"
MV.data$ICR4 <- "FALSE"
MV.data[MV.data$TCGA.PAM50.RMethod.RNASeq=="Basal-like",c("Basal")] <- "TRUE"
MV.data[MV.data$Cluster.DBGS3.FLTR.RNSeq=="ICR4",c("ICR4")] <- "TRUE"
write.csv (MV.data,"./3 ANALISYS/Multivariate analysis/ICRCluster_TP53_MAPK_IMS.csv")
