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
  required.packages <- c("Hmisc","HGNChelper")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  #library (xlsx) #xlsx needs java installed
  library (Hmisc)
  #setwd("~/Dropbox/BREAST_QATAR/")
  setwd("/mnt3/wouter/BREAST-QATAR/")

# Set Parameters
  Cancersets        = "ALL"
  Geneset           = "DBGS3.FLTR"             # SET GENESET HERE
  DL.Method         = "BIOLINKS"               #Choose "ASSEMBLER" or "BIOLINKS"
  sample.types      = "Selected"               #Alternatives TP , TP_TM , Selected   
  
# DO ALL
  TCGA.cancersets <- read.csv ("./2 DATA/TCGA.datasets.csv")
  if (Cancersets == "ALL") { 
    Cancersets = gsub("\\]","",gsub(".*\\[","",TCGA.cancersets$Cancername))
  }
  N.sets = length(Cancersets)
  for (i in 1:N.sets) {
    Cancerset = Cancersets[i]
    if (Cancerset %in% c("LAML","FPPP")) {next}
 
# Load data files
  Parent.Cancerset  = substring(Cancerset,1,4)
  Parent.Geneset    = substring(Geneset,1,5)
  #clinical data
  ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_",DL.Method,"_subset_clinicaldata.csv"))                       # Clinical data including IMS
  ClinicalData.subset <- ClinicalData.subset[-which(is.na(ClinicalData.subset$X)),]
  rownames(ClinicalData.subset) <- ClinicalData.subset$X 
  ClinicalData.subset$X <-NULL
  #RNASeq data
  load (paste0("./2 DATA/SUBSETS/",DL.Method,"/",Parent.Cancerset,"/TCGA.",Parent.Cancerset,".RNASeq.",sample.types,".subset.",Parent.Geneset,".RData"))
  #Consensus clustering data
  CC.RNASeq <- read.csv (paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000/",
                                Cancerset,".TCGA.",DL.Method,".EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"))                                                # Cluster assignment
  rownames(CC.RNASeq) <- CC.RNASeq$X 
  CC.RNASeq$X <- NULL
  colnames(CC.RNASeq) <- c("PatientID",paste0("Cluster.",Geneset,".RNSeq"))
  CC.RNASeq$PatientID <-NULL
  #immunoscore data
  immunoscore <- read.csv (paste0("./3 ANALISYS/IMMUNOSCORE/immunoscore.TCGA.",DL.Method,".",Parent.Cancerset,".",Parent.Geneset,".csv"))                            # Immunoscore
  rownames(immunoscore) <- immunoscore$X
  immunoscore$X <- NULL
  #mutation frequnecy data
  load(paste0("./3 ANALISYS/Mutations/",DL.Method,"/",Parent.Cancerset,"/Mutation.Data.TCGA.",Parent.Cancerset,".All.DBGS3.FLTR.Frequencies.RDATA"))
  mutation.freq <- Mutation.Frequency.Patient[,c("Patient_ID","Freq.NonSilent","Freq.All")]
  row.names(mutation.freq)<-mutation.freq$Patient_ID
  mutation.freq$Patient_ID <- NULL
  
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
  

# export data to Rdata , csv and excel
save (Master.file,file = paste0("./3 ANALISYS/MASTER FILES/",DL.Method,"/TCGA.",DL.Method,".",Cancerset,".RNASeq_subset_",Geneset,".Master.Summary.Rdata"))
write.csv (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/",DL.Method,"/TCGA.",DL.Method,".",Cancerset,".RNASeq_subset_",Geneset,".Master.Summary.csv"),row.names = TRUE);
#write.xlsx (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.xlsx"), sheetName ="RNASeq ISGS.Master data", row.names=TRUE)

#multivariate analysis
#Master.file$Cluster.DBGS3.FLTR.RNSeq
#MV.data <- Master.file[,c("TCGA.PAM50.RMethod.RNASeq","GOF_mutation","Cluster.DBGS3.FLTR.RNSeq")]
#MV.data$MAPKMUT <- Mutstat$MAP2K4.MAP3K1[match(rownames(MV.data),Mutstat$Sample)]
#MV.data <- MV.data[complete.cases(MV.data),]
#V.data$Basal <- "FALSE"
#MV.data$ICR4 <- "FALSE"
#MV.data[MV.data$TCGA.PAM50.RMethod.RNASeq=="Basal-like",c("Basal")] <- "TRUE"
#MV.data[MV.data$Cluster.DBGS3.FLTR.RNSeq=="ICR4",c("ICR4")] <- "TRUE"
#write.csv (MV.data,"./3 ANALISYS/Multivariate analysis/ICRCluster_TP53_MAPK_IMS.csv")
print (paste0("Masterfile for ",Cancerset," saved..."))
}
