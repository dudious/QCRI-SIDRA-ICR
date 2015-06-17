################################################################
###
### This scripts prepares the TCGA CNV data for GISTIC analysis
###
###############################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR")

# Set Parameters
Cancerset   <- "BRCA"
Geneset     <- "DBGS1"

# Load Data
CNV.Data <- read.table(paste0("./2 DATA/TCGA CNV/TCGA CNV_",Cancerset,"_ASSEMBLER/",Cancerset,".CN.TCGA.ASSEMBLER.DATA.nocnv.hg19.txt"),as.is=TRUE,header = TRUE)
Consensuss.class<-read.delim(file=paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",
                               Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",
                               Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),sep=",")

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

