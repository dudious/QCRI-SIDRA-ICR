# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

#install missing packages
required.packages <- c("gplots")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("edgeR")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
#source("http://bioconductor.org/biocLite.R")
if(length(missing.packages)) biocLite(missing.packages)

#Load the packages
library(gplots)
library(edgeR)

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

#load the functions

TimeUse <-function(func){ 
  time1 <- proc.time() 
  result <-func
  time2 <- proc.time()
  show(timeelapsed <- time2-time1)
}

edgeRDEA4GSEAOff <- function(geneCounts, groups, baseLine){
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  DEA <- DGEList(geneCounts, group = groups)
  DEA <- calcNormFactors(DEA)
  keep <- rowSums(cpm(DEA)>1)>=10
  DEA <- DEA[keep,]
  DEA$samples$lib.size <- colSums(DEA$counts)
  message("Estimating GLM Common Dispersion")
  TimeUse(DEA <- estimateGLMCommonDisp(DEA, design = design))
  DEA <- estimateGLMTrendedDisp(DEA, design = design)
  message("Estimating GLM Tagwise Dispersion")
  TimeUse(DEA <- estimateGLMTagwiseDisp(DEA, design = design))
  message("Fitting the GLM Model")
  TimeUse(fit <- glmFit(DEA, design = design))
  contrast <- rep(1, 2)
  blId <- which(colnames(design) == baseLine)
  contrast[blId] <- -1
  message("Maximum Likelihood Estimate")
  TimeUse(test <- glmLRT(fit, contrast = contrast))
  ans <- topTags(test, n = nrow(DEA$counts))
  ans <- ans$table
  diffGenes <- abs(ans[, "logFC"]) >= 1 & ans[, "FDR"] <= 0.01
  detags <- rownames(ans)[diffGenes]
  pdf("smearPlot.pdf")
  plotSmear(test,de.tags=detags)
  abline(h=c(-1,1),col="blue")
  dev.off()
  return(ans)
}

##############################################################################

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
Signif.tresh = 0.05
fch.tresh = 1
method = "edgeR"
IMS_filter = "Luminal"
RNASeq_filter = 0.25


#load MAPKMUT data
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.rdata",verbose = TRUE)

#load and subset RNAseq data
load("./2 DATA/TCGA RNAseq/RNASeq_BRCA_EDASeq/BRCA.RNASeq.TCGA.ASSEMBLER.NORMALIZED.RData")
selectedGenes <- which(rowMeans(RNASeq.NORM) >= quantile(rowMeans(RNASeq.NORM), RNASeq_filter))
RNASeq.NORM <- RNASeq.NORM[selectedGenes,]
dim(RNASeq.NORM)
RNASeq.subset <- as.data.frame(t(RNASeq.NORM[,])) #select all / no subset used

#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#Add MAPK status (MUT/WT) (non-silent)
RNASeq.subset$MAPK <- FinalClassification$MAP2K4.MAP3K1[match(rownames(RNASeq.subset),FinalClassification$Sample)]
RNASeq.subset <- RNASeq.subset[-which(is.na(RNASeq.subset$MAPK)),]

#Filter luminal only
if (IMS_filter=="Luminal"){
  subtype <- ClinicalData.subset[rownames(RNASeq.subset),"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
  luminal <- subtype[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A" | subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B",,drop=FALSE]
  RNASeq.subset <- RNASeq.subset[which(rownames(RNASeq.subset) %in% rownames(luminal)),]
}

##edgeR
if (method == "edgeR"){
  groups <- rep("",nrow(RNASeq.subset))
  names(groups)<- rownames(RNASeq.subset)

  #MAPK MUT vs WT
  G1<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]=="MUT",])
  G2<- rownames(RNASeq.subset[RNASeq.subset[,ncol(RNASeq.subset)]=="WT",])
  groups[G1]<-"G1"
  groups[G2]<-"G2"
  groups<-as.factor(groups)
  
  RNASeq.subset <- t(RNASeq.subset[,-ncol(RNASeq.subset)])  #remove grouping column
  
  system.time(DEGs <- edgeRDEA4GSEAOff(RNASeq.subset[,names(groups)], groups = groups,baseLine = "G2"))
  DEGs<-DEGs[order((DEGs$logFC),decreasing = T),]
  
  write.csv(DEGs,file=paste0("./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEG_BRCA_BSF2_MAPX_",IMS_filter,".csv"))
  signif.genes <- DEGs[DEGs$FDR<=Signif.tresh & abs(DEGs$logFC)>=fch.tresh,,drop=FALSE]
  save(DEGs,signif.genes,RNASeq.subset,file=paste0("./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEG_MAPX_",IMS_filter,".RDATA"))
  
  #reconstitute data for heatmap matrix 
  #RE-load and subset RNAseq data
  load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/FinalClassification.rdata")
  load(file=paste0("./3 ANALISYS/DIFFERENTIAL EXPRESSION/wout/DEG_MAPX_",IMS_filter,".RDATA"))
  RNASeq.subset <- as.data.frame(RNASeq.subset[rownames(signif.genes),])
  RNASeq.subset <- as.data.frame(t(RNASeq.subset))
  RNASeq.subset$MAPX_state <- FinalClassification$MAP2K4.MAP3K1[match(rownames(RNASeq.subset),FinalClassification$Sample)]
  
  #ordering of the clusters
  RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$MAPX_state,levels = c("MUT","WT"))),]     
  MAPX.state <- as.data.frame(RNASeq.subset$MAPX_state) 
  colnames(MAPX.state) <- "MAPX_state"
  RNASeq.subset$MAPX_state <- NULL
  
  # Heatmap 2 (simple no extra annotations)
  RNASeq.subset <- log(as.matrix(RNASeq.subset)+1,2)
  mode(RNASeq.subset) <- "numeric"
  patientcolors <- MAPX.state
  levels (patientcolors$MAPX_state) <- c(levels (patientcolors$MAPX_state),c("red","yellow"))  #Aply color scheme to patients
  patientcolors$MAPX_state[patientcolors$MAPX_state=="MUT"] <- "red"
  patientcolors$MAPX_state[patientcolors$MAPX_state=="WT"] <- "yellow"
  patientcolors <- as.character(patientcolors$MAPX_state)
  my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
  my.colors = unique(c(seq(-2,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))
  png(paste0("./4 FIGURES/Heatmaps/Heatmap.MAPX-MUT.RNASeq_Log2.TCGA.BRCA-BSF2.",IMS_filter,".DBGS3.DENDO.png"),res=600,height=6,width=6,unit="in")     # set filename
  heatmap.2(t(RNASeq.subset),
            main = paste0("RNASeq.log2 - MAPX MUT"),
            col=my.palette,                   #set color sheme RED High, GREEN low
            #breaks=my.colors,
            ColSideColors=patientcolors,      #set goup colors                 
            key=TRUE,
            symm=FALSE,
            symkey=FALSE,
            symbreaks=TRUE,             
            scale="row",
            density.info="none",
            trace="none",
            labCol=FALSE,
            cexRow=0.8,cexCol=0.1,margins=c(2,7),
            Colv=TRUE)
  par(lend = 1)
  legend("topright",legend = c("Mutated","Wild-type"),
         col = c("red","yellow"),lty= 1,lwd = 5,cex = 0.7)
  dev.off()
  
  #heatmap without dendogram
  
  png(paste0("./4 FIGURES/Heatmaps/Heatmap.MAPX-MUT.RNASeq_log2.TCGA.BRCA-BSF2.",IMS_filter,".DBGS3.png"),res=600,height=6,width=6,unit="in")     # set filename
  heatmap.2(t(RNASeq.subset),
            main = paste0("RNASeq.log2 - MAPX MUT"),
            col=my.palette,                   #set color sheme RED High, BLUE low
            #breaks=my.colors,
            ColSideColors=patientcolors,      #set goup colors                 
            key=TRUE,
            symm=FALSE,
            symkey=FALSE,
            symbreaks=TRUE,             
            scale="row",
            density.info="none",
            trace="none",
            labCol=FALSE,
            cexRow=0.8,cexCol=0.1,margins=c(2,7),
            Colv=FALSE)
  par(lend = 1)
  legend("topright",legend = c("Mutated","Wild-type"),
         col = c("red","yellow"),lty= 1,lwd = 5,cex = 0.7)
  dev.off()
  
  
  
  #lookup subtype
  subtype <- ClinicalData.subset[rownames(RNASeq.subset),"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
  
  levels (subtype$TCGA.PAM50.RMethod.RNASeq) <- c(levels (subtype$TCGA.PAM50.RMethod.RNASeq),c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3"))
  subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal A"]     <- "#eaff00"    
  subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Luminal B"]     <- "#00c0ff"     
  subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Basal-like"]    <- "#da70d6"   
  subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="HER2-enriched"] <- "#daa520"            
  subtype$TCGA.PAM50.RMethod.RNASeq[subtype$TCGA.PAM50.RMethod.RNASeq=="Normal-like"]   <- "#d3d3d3"
  subtypecolors <- as.character(subtype$TCGA.PAM50.RMethod.RNASeq)
  #lookup cluster
  cluster <- Consensus.class[rownames(RNASeq.subset),"Cluster",drop=FALSE]
  levels (cluster$Cluster) <- c(levels (cluster$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))        # Apply color scheme to patients
  cluster$Cluster[cluster$Cluster=="ICR4"] <- "#FF0000"
  cluster$Cluster[cluster$Cluster=="ICR3"] <- "#FFA500"
  cluster$Cluster[cluster$Cluster=="ICR2"] <- "#00FF00"
  cluster$Cluster[cluster$Cluster=="ICR1"] <- "#0000FF"
  clustercolors <- as.character(cluster$Cluster)
  
  color.matrix <- as.matrix(rbind (patientcolors,subtypecolors,clustercolors))
  rownames(color.matrix) <- c("MAP state","Subtype","Cluster")

  png(paste0("./4 FIGURES/Heatmaps/Heatmap.MAPX-MUT.RNASeq_log2.TCGA.BRCA-BSF2.",IMS_filter,".DBGS3.anotated.png"),res=600,height=6,width=6,unit="in")     # set filename
  heatmap.3(t(RNASeq.subset),
            main = "RNASeq.log2 - MAPX MUT LUM",
            col=my.palette,                                     # set color scheme RED High, GREEN low
            ColSideColors=t(color.matrix),                         # set goup colors
            key=TRUE,
            symm=FALSE,
            symkey=FALSE,
            symbreaks=TRUE,             
            scale="row",
            density.info="none",
            trace="none",
            labCol=FALSE,
            cexRow=0.8,cexCol=0.1,margins=c(9,9),
            Colv=FALSE,Rowv=TRUE                           # reorder row/columns by dendogram
  )
  par(lend = 1)
  #legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
  #       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
  legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1","","Mutated","Wild Type"),
         col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue","white","red","yellow"),lty= 1,lwd = 5,cex = 0.3)
  dev.off()
  
  png(paste0("./4 FIGURES/Heatmaps/Heatmap.MAPX-MUT.RNASeq_log2.TCGA.BRCA-BSF2.",IMS_filter,".DBGS3.anotated.DENDO.png"),res=600,height=6,width=6,unit="in")     # set filename
  heatmap.3(t(RNASeq.subset),
            main = "RNASeq.log2 - MAPX MUT LUM",
            col=my.palette,                                     # set color scheme RED High, GREEN low
            ColSideColors=t(color.matrix),                         # set goup colors
            key=TRUE,
            symm=FALSE,
            symkey=FALSE,
            symbreaks=TRUE,             
            scale="row",
            density.info="none",
            trace="none",
            labCol=FALSE,
            cexRow=0.8,cexCol=0.1,margins=c(9,9),
            Colv=TRUE,Rowv=TRUE                           # reorder row/columns by dendogram
  )
  par(lend = 1)
  #legend("left",legend = c("ICR4","ICR3","ICR2","ICR1"),
  #       col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 1.3)
  legend("topright",legend = c("Luminal A","Luminal B","Basal-like","HER2-enriched","Normal-like","","ICR4","ICR3","ICR2","ICR1","","Mutated","Wild Type"),
         col = c("#eaff00","#00c0ff","#da70d6","#daa520","#d3d3d3","white","red","orange","green","blue","white","red","yellow"),lty= 1,lwd = 5,cex = 0.3)
  dev.off()
}

if (method == "kruskal"){
  
  # filter for significant genes
  RNASeq.subset <- cbind(RNASeq.subset,cluster.assignment$Cluster)
  p.value.matrix <- as.data.frame(matrix(vector(),(ncol(RNASeq.subset)-1),2))
  colnames(p.value.matrix) <- c("Gene","p.value")
  p.value.matrix$Gene <- colnames(RNASeq.subset)[-length(colnames(RNASeq.subset))]
  rownames(p.value.matrix) <- p.value.matrix$Gene
  p.value.matrix$Gene <- NULL
  for (i in 1:(ncol(RNASeq.subset)-1)) {
    Gene = colnames(RNASeq.subset)[i]
    test = kruskal.test(RNASeq.subset[,i]~RNASeq.subset[,ncol(RNASeq.subset)])
    p.value.matrix[Gene,1] = test$p.value
    
  }
  signif.genes <- p.value.matrix[p.value.matrix$p.value<Signif.tresh,,drop=FALSE]
  RNASeq.subset <- as.data.frame(RNASeq.subset[,rownames(signif.genes)])
  
  #ordering of the clusters
  RNASeq.subset$Cluster <- cluster.assignment$Cluster[match(rownames(RNASeq.subset),cluster.assignment$PatientID)]
  RNASeq.subset <- RNASeq.subset[order(factor(RNASeq.subset$Cluster,levels = c("ICR4","ICR3","ICR2","ICR1"))),]     
  RNASeq.subset$Cluster <- NULL
  
  #re-order the labels
  cluster.assignment <- cluster.assignment[rownames(RNASeq.subset),]
  
  # Heatmap 2 (simple no extra annotations)
  RNASeq.subset <- as.matrix(RNASeq.subset)
  mode(RNASeq.subset) <- "numeric"
  patientcolors <- cluster.assignment
  levels (patientcolors$Cluster) <- c(levels (patientcolors$Cluster),c("#FF0000","#FFA500","#00FF00","#0000FF"))  #Aply color scheme to patients
  patientcolors$Cluster[patientcolors$Cluster=="ICR4"] <- "#FF0000"
  patientcolors$Cluster[patientcolors$Cluster=="ICR3"] <- "#FFA500"
  patientcolors$Cluster[patientcolors$Cluster=="ICR2"] <- "#00FF00"
  patientcolors$Cluster[patientcolors$Cluster=="ICR1"] <- "#0000FF"
  #patientcolors$Cluster <- droplevels(patientcolors$Cluster)
  patientcolors <- as.character(patientcolors$Cluster)
  my.palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 297)
  my.colors = unique(c(seq(-2,-0.5,length=100),seq(-0.5,0.5,length=100),seq(0.5,4,length=100)))
  png(paste0("./4 FIGURES/Heatmaps/Heatmap.Neo-antigens.kruskal.RNASeq.TCGA.BRCA-BSF2.",IMS_filter,".DBGS3.png"),res=600,height=6,width=6,unit="in")     # set filename
  heatmap.2(t(RNASeq.subset),
            main = paste0("Heatmap RNASeq - neo-antigens"),
            col=my.palette,                   #set color sheme RED High, GREEN low
            breaks=my.colors,                                 
            ColSideColors=patientcolors,      #set goup colors                 
            key=TRUE,
            symm=FALSE,
            symkey=FALSE,
            symbreaks=TRUE,             
            scale="row",
            density.info="none",
            trace="none",
            labCol=FALSE,
            cexRow=0.8,cexCol=0.1,margins=c(2,7),
            Colv=TRUE)
  par(lend = 1)
  legend("topright",legend = c("ICR4","ICR3","ICR2","ICR1"),
         col = c("red","orange","green","blue"),lty= 1,lwd = 5,cex = 0.7)
  dev.off()
  
}