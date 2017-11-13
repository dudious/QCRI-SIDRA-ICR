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
  load (paste0("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData"))
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
  #mutation GOF status
  load("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.All.DBGS3.FLTR.Mutation.Matrixes.NonSilent.oncoplot.Rdata")
  load("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.All.DBGS3.FLTR.Mutation.Matrixes.NonSilent.Rdata")
  genes.mutations.selected[is.na(genes.mutations.selected)] <- 0
  merged.matrix[is.na(merged.matrix)] <- "NONE"
  # merge data
  Master.file <- merge (ClinicalData.subset,CC.RNASeq,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,immunoscore,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,mutation.freq,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,genes.mutations.selected,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,merged.matrix[,c("CTCF","FCGBP","MAP2K4","MAP3K1","TP53")],by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,t(RNASeq.NORM_Log2),by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
 
  

# export data to txt and excel
write.csv (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.AllG.csv"),row.names = TRUE);
#write.xlsx (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.AllG.xlsx"), sheetName ="RNASeq ISGS.Master data", row.names=TRUE)
