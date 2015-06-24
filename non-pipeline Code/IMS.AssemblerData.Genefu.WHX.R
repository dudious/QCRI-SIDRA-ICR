#################################################################
###
### This Script uses genefu to predicit IMS for the TCGA Micro Array 
### patient population. (Hu,Sorlie and Parker(PAM50)).
### It uses as input 
### Micro array Data retrieved using TCGA Assembbler
### source data :
### "./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.rda"
### Results are saved in
### ./3 ANALISYS/IMS/Genefu IMS/
### File to use :
### "PredictionTable.ASSEMBLER.csv"
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR")
  ## dependencies
     required.packages <- c("xlsx")
     missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
     if(length(missing.packages)) install.packages(missing.packages)
     required.packages.BioC <- c("genefu","org.Hs.eg.db")
     missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
     source("http://bioconductor.org/biocLite.R")
     if(length(missing.packages)) biocLite(missing.packages)
  library(genefu);
  library(org.Hs.eg.db);
  library (xlsx);

# Load Data 
  load ("./2 DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.Rdata")
  data(pam50)
# build ID matrix from Gene symbols (HGNC)
  HGNC <- colnames(agilentData);
  ID.selection <- c("ENTREZID");
## exrta ID.selection <- c("ENTREZID","ENSEMBL","UNIGENE","REFSEQ");
   ID.table <- select(org.Hs.eg.db, keys=HGNC, columns=ID.selection, keytype="SYMBOL");
   colnames(ID.table)[1] <- "probe"
   colnames(ID.table)[2] <- "EntrezGene.ID"

# cluster samples by genelist (new model)
  #whx.model <- intrinsic.cluster(data=agilentData, annot=ID.table,
  #do.mapping = TRUE,
  #std = "robust",
  #intrinsicg=pam50$centroids.map[ ,c("probe", "EntrezGene.ID")],
  #number.cluster = 5, mins = 5,
  #method.centroids = c("mean", "median", "tukey"),
  #filen="./DATA/whx.model", verbose = TRUE);

# predict subclass based on PAM50
  pam50.prediction <- intrinsic.cluster.predict(sbt.model=pam50,
  data=agilentData, annot=ID.table,
  do.mapping = TRUE,
  do.prediction.strength = FALSE, verbose = TRUE);

# predict subclass based on Sorlie
  Sorlie.prediction <- intrinsic.cluster.predict(sbt.model=ssp2003,
  data=agilentData, annot=ID.table,
  do.mapping = TRUE,
  do.prediction.strength = FALSE, verbose = TRUE);

# predict subclass based on Hu
  Hu.prediction <- intrinsic.cluster.predict(sbt.model=ssp2006,
  data=agilentData, annot=ID.table,
  do.mapping = TRUE,
  do.prediction.strength = FALSE, verbose = TRUE);

# predict subclass based on whx.model
## whx.prediction <- intrinsic.cluster.predict(sbt.model=whx.model,
##  data=agilentData, annot=ID.table,
##  do.mapping = TRUE,
##  do.prediction.strength = FALSE, verbose = TRUE);

# subtype distribution
  print ("pam50 Distribution");
  print (table(pam50.prediction$subtype));
  print ("Sorlie Distribution");
  print(table(Sorlie.prediction$subtype));
  print ("Hu Distribution");
  print(table(Hu.prediction$subtype));
## print(table(whx.prediction$subtype));

# SAVE
  PredictionTable <- data.frame (pam50.prediction$subtype, Hu.prediction$subtype, Sorlie.prediction$subtype);
  PredictionTable <- cbind(Sample.ID = rownames(PredictionTable), PredictionTable); 
  rownames(PredictionTable) <- NULL;
  Hu.prediction.proba.table <- Hu.prediction$subtype.proba;
  Hu.prediction.proba.table <- cbind(Sample.ID = rownames(Hu.prediction.proba.table), Hu.prediction.proba.table);
  rownames(Hu.prediction.proba.table) <- NULL;
  Sorlie.prediction.proba.table <- Sorlie.prediction$subtype.proba;
  Sorlie.prediction.proba.table <- cbind(Sample.ID = rownames(Sorlie.prediction.proba.table), Sorlie.prediction.proba.table);
  rownames(Sorlie.prediction.proba.table) <- NULL;
  pam50.prediction.proba.table <- pam50.prediction$subtype.proba;
  pam50.prediction.proba.table <- cbind(Sample.ID = rownames(pam50.prediction.proba.table), pam50.prediction.proba.table);
  rownames(pam50.prediction.proba.table) <- NULL;
  ResultFilename <- "./3 ANALISYS/IMS/Genefu IMS/PredictionTable.ASSEMBLER";
  write.csv(PredictionTable, file = paste0 (ResultFilename, ".csv",sep = ""), sep="\t", quote=FALSE, row.names=FALSE);
  write.xlsx  (PredictionTable, file = paste0 (ResultFilename, ".xlsx",sep = ""), sheetName ="IMS Prediction overview", row.names=FALSE);
  write.xlsx  (Hu.prediction.proba.table, file = paste0 (ResultFilename, ".xlsx",sep = ""), sheetName ="Hu IMS probability ", row.names=FALSE,append =TRUE);
  write.xlsx  (Sorlie.prediction.proba.table, file = paste0 (ResultFilename, ".xlsx",sep = ""), sheetName ="Sorlie IMS probability ", row.names=FALSE,append =TRUE);
  write.xlsx  (pam50.prediction.proba.table, file = paste0 (ResultFilename, ".xlsx",sep = ""), sheetName ="pam50 IMS probability ", row.names=FALSE,append =TRUE);
  print (paste0 ("Results are saved in " , ResultFilename , ".xlsx and " , ResultFilename , ".csv.",sep = ""));

