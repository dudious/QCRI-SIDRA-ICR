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
  #library (xlsx) #xlsx needs java installed
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
  load("./3 ANALISYS/Mutations/BRCA.BSF2/Mutation.Data.TCGA.BRCA.BSF2.DBGS3.FLTR.Frequencies.RDATA")
  mutation.freq <- Mutation.Frequency.Patient[,c("Patient_ID","Freq.NonSilent","Freq.All")]
  row.names(mutation.freq)<-mutation.freq$Patient_ID
  mutation.freq$Patient_ID <- NULL
  #TP53.mutaion status
  matrix.type    = "NonSilent"    # Alterantives "Any" , "Missense", "NonSilent"
  IMS.filter     = "All"          # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
  load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
  TP53.NS.mutations <- as.data.frame(merged.matrix[,"TP53",drop=FALSE],stringsAsFactors = FALSE)
  TP53.NS.mutations$TP53[TP53.NS.mutations$TP53 == "HOMDEL"] <- "WT"
  TP53.NS.mutations$TP53[TP53.NS.mutations$TP53 == "AMP"] <- "WT"
  TP53.NS.mutations$TP53[TP53.NS.mutations$TP53 == "MUT;HOMDEL"] <- "MUT"
  TP53.NS.mutations$TP53[TP53.NS.mutations$TP53 == "MUT;AMP"] <- "MUT"
  TP53.NS.mutations$TP53[is.na(TP53.NS.mutations$TP53)] <- "WT"
  colnames(TP53.NS.mutations)<-"TP53.NS.mutations"
  matrix.type    = "Any"    # Alterantives "Any" , "Missense", "NonSilent"
  IMS.filter     = "All"    # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
  load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
  TP53.ALL.mutations <- as.data.frame(merged.matrix[,"TP53",drop=FALSE],stringsAsFactors = FALSE)
  TP53.ALL.mutations$TP53[TP53.ALL.mutations$TP53 == "HOMDEL"] <- "WT"
  TP53.ALL.mutations$TP53[TP53.ALL.mutations$TP53 == "AMP"] <- "WT"
  TP53.ALL.mutations$TP53[TP53.ALL.mutations$TP53 == "MUT;HOMDEL"] <- "MUT"
  TP53.ALL.mutations$TP53[TP53.ALL.mutations$TP53 == "MUT;AMP"] <- "MUT"
  TP53.ALL.mutations$TP53[is.na(TP53.ALL.mutations$TP53)] <- "WT"
  colnames(TP53.ALL.mutations)<-"TP53.ALL.mutations"
  
  #MAP2K4.mutaion status
  matrix.type    = "NonSilent"    # Alterantives "Any" , "Missense", "NonSilent"
  IMS.filter     = "All"          # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
  load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
  MAP2K4.NS.mutations <- as.data.frame(merged.matrix[,"MAP2K4",drop=FALSE],stringsAsFactors = FALSE)
  MAP2K4.NS.mutations$MAP2K4[MAP2K4.NS.mutations$MAP2K4 == "HOMDEL"] <- "WT"
  MAP2K4.NS.mutations$MAP2K4[MAP2K4.NS.mutations$MAP2K4 == "AMP"] <- "WT"
  MAP2K4.NS.mutations$MAP2K4[MAP2K4.NS.mutations$MAP2K4 == "MUT;HOMDEL"] <- "MUT"
  MAP2K4.NS.mutations$MAP2K4[MAP2K4.NS.mutations$MAP2K4 == "MUT;AMP"] <- "MUT"
  MAP2K4.NS.mutations$MAP2K4[is.na(MAP2K4.NS.mutations$MAP2K4)] <- "WT"
  colnames(MAP2K4.NS.mutations)<-"MAP2K4.NS.mutations"
  matrix.type    = "Any"    # Alterantives "Any" , "Missense", "NonSilent"
  IMS.filter     = "All"    # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
  load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
  MAP2K4.ALL.mutations <- as.data.frame(merged.matrix[,"MAP2K4",drop=FALSE],stringsAsFactors = FALSE)
  MAP2K4.ALL.mutations$MAP2K4[MAP2K4.ALL.mutations$MAP2K4 == "HOMDEL"] <- "WT"
  MAP2K4.ALL.mutations$MAP2K4[MAP2K4.ALL.mutations$MAP2K4 == "AMP"] <- "WT"
  MAP2K4.ALL.mutations$MAP2K4[MAP2K4.ALL.mutations$MAP2K4 == "MUT;HOMDEL"] <- "MUT"
  MAP2K4.ALL.mutations$MAP2K4[MAP2K4.ALL.mutations$MAP2K4 == "MUT;AMP"] <- "MUT"
  MAP2K4.ALL.mutations$MAP2K4[is.na(MAP2K4.ALL.mutations$MAP2K4)] <- "WT"
  colnames(MAP2K4.ALL.mutations)<-"MAP2K4.ALL.mutations"
  
  #MAP3K1.mutaion status
  matrix.type    = "NonSilent"    # Alterantives "Any" , "Missense", "NonSilent"
  IMS.filter     = "All"          # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
  load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
  MAP3K1.NS.mutations <- as.data.frame(merged.matrix[,"MAP3K1",drop=FALSE],stringsAsFactors = FALSE)
  MAP3K1.NS.mutations$MAP3K1[MAP3K1.NS.mutations$MAP3K1 == "HOMDEL"] <- "WT"
  MAP3K1.NS.mutations$MAP3K1[MAP3K1.NS.mutations$MAP3K1 == "AMP"] <- "WT"
  MAP3K1.NS.mutations$MAP3K1[MAP3K1.NS.mutations$MAP3K1 == "MUT;HOMDEL"] <- "MUT"
  MAP3K1.NS.mutations$MAP3K1[MAP3K1.NS.mutations$MAP3K1 == "MUT;AMP"] <- "MUT"
  MAP3K1.NS.mutations$MAP3K1[is.na(MAP3K1.NS.mutations$MAP3K1)] <- "WT"
  colnames(MAP3K1.NS.mutations)<-"MAP3K1.NS.mutations"
  matrix.type    = "Any"    # Alterantives "Any" , "Missense", "NonSilent"
  IMS.filter     = "All"    # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
  load(paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".oncoplot.Rdata"))
  MAP3K1.ALL.mutations <- as.data.frame(merged.matrix[,"MAP3K1",drop=FALSE],stringsAsFactors = FALSE)
  MAP3K1.ALL.mutations$MAP3K1[MAP3K1.ALL.mutations$MAP3K1 == "HOMDEL"] <- "WT"
  MAP3K1.ALL.mutations$MAP3K1[MAP3K1.ALL.mutations$MAP3K1 == "AMP"] <- "WT"
  MAP3K1.ALL.mutations$MAP3K1[MAP3K1.ALL.mutations$MAP3K1 == "MUT;HOMDEL"] <- "MUT"
  MAP3K1.ALL.mutations$MAP3K1[MAP3K1.ALL.mutations$MAP3K1 == "MUT;AMP"] <- "MUT"
  MAP3K1.ALL.mutations$MAP3K1[is.na(MAP3K1.ALL.mutations$MAP3K1)] <- "WT"
  colnames(MAP3K1.ALL.mutations)<-"MAP3K1.ALL.mutations"
  
  #MAPK STATUS from Ines
  Mutstat <- read.csv ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.csv")
  rownames(Mutstat)<-Mutstat$X
  Mutstat$X <- NULL
  Mutstat$Sample <- NULL

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
  Master.file <- merge (Master.file,TP53.NS.mutations,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,TP53.ALL.mutations,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,MAP2K4.NS.mutations,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,MAP2K4.ALL.mutations,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,MAP3K1.NS.mutations,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,MAP3K1.ALL.mutations,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  Master.file <- merge (Master.file,Mutstat,by="row.names",all.x=TRUE)
  rownames(Master.file) <- Master.file$Row.names
  Master.file$Row.names <- NULL
  

# export data to Rdata , csv and excel
save (Master.file,file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.Summary.Rdata"))
write.csv (Master.file, file = paste0("./3 ANALISYS/MASTER FILES/TCGA.",Cancerset,".RNASeq_subset_",Geneset,".Master.Summary.csv"),row.names = TRUE);
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
