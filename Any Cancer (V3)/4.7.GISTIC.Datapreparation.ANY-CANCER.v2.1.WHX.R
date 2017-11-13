################################################################
###
### This scripts prepares the TCGA CNV data for GISTIC analysis
###
###############################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

# Set Parameters
Cancerset   <- "BRCA"         #use combined cancerset not -GA or -hiseq
Geneset     <- "DBGS3.FLTR"
BRCA.Filter <- "PCF"          # "PCF" or "BSF" Pancer or Breast specific

# check for split dataset (GA/hiseq)
split = "FALSE"
Download.path <- paste0("./2 DATA/TCGA RNAseq/RNASeq_",Cancerset,"_ASSEMBLER/")
file.hiseq    <- list.files(Download.path,full.names = TRUE,pattern = "hiseq.DATA.txt" )
hiseq         <- length(file.hiseq)
file.GA       <- list.files(Download.path,full.names = TRUE,pattern = "GA.DATA.txt" )
GA            <- length(file.GA)
if (hiseq+GA == 2) {split = "TRUE"}

# Load Data
CNV.Data <- read.table(paste0("./2 DATA/TCGA CNV/TCGA CNV_",Cancerset,"_ASSEMBLER/",Cancerset,".CN.TCGA.ASSEMBLER.DATA.nocnv_hg19.txt"),as.is=TRUE,header = TRUE)
if (Cancerset == "BRCA"){
  if (substring(Geneset,7,10)=="FLTR"){
    Cancerset <- paste0(Cancerset,".",BRCA.Filter)
  }
}
if (split == "TRUE"){
  Geneset.GA = paste0(Geneset,"-GA")
  Consensuss.class.GA<-read.delim(file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"-GA/",
                                           Cancerset,"-GA.TCGA.EDASeq.k7.",Geneset,".reps5000/",
                                           Cancerset,"-GA.TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),sep=",")
  Geneset.hiseq = paste0(Geneset,"-GA")
  Consensuss.class.hiseq<-read.delim(file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"-hiseq/",
                                              Cancerset,"-hiseq.TCGA.EDASeq.k7.",Geneset,".reps5000/",
                                              Cancerset,"-hiseq.TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),sep=",")
  Consensuss.class <- unique(rbind(Consensuss.class.hiseq,Consensuss.class.GA))
}
 if (split == "FALSE"){
  Consensuss.class<-read.delim(file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",
                               Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",
                               Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),sep=",")
}

# merge CNV data with Cluster assignment
CNV.Data$PatientID <- substr(CNV.Data$Sample,1,12)
CNV.Data$Cluster <- Consensuss.class$Group[match(CNV.Data$PatientID,Consensuss.class$PatientID)]

# Split and save
dir.create(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset), showWarnings = FALSE)
dir.create(paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC"), showWarnings = FALSE)

write.table(CNV.Data[CNV.Data$Cluster=="ICR1",1:6],
            file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_ICR1_seg.txt"),
            row.names=FALSE,sep="\t",quote=FALSE)
write.table(CNV.Data[CNV.Data$Cluster=="ICR2",1:6],
            file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_ICR2_seg.txt"),
            row.names=FALSE,sep="\t",quote=FALSE)
write.table(CNV.Data[CNV.Data$Cluster=="ICR3",1:6],
            file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_ICR3_seg.txt"),
            row.names=FALSE,sep="\t",quote=FALSE)
write.table(CNV.Data[CNV.Data$Cluster=="ICR4",1:6],
            file=paste0("./3 ANALISYS/GISTIC/GISTIC.",Cancerset,"/",Geneset,".FILE.FOR.GISTIC/",Geneset,".Cluster_ICR4_seg.txt"),
            row.names=FALSE,sep="\t",quote=FALSE)

