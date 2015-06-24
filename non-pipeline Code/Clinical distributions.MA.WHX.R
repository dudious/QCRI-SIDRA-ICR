#################################################################
###
### This Script Annalyses the distribution of clinical parameters
### in the TCGA Micro Array patient population. It generates figures
### relating to these annalyses. 
### Data is saved :
### ./3 ANALYSIS/Clinical Information/BRCA/Distribution/
### Figures are saved :
### ./FIGURES
###
#################################################################

# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR/")
  #Dependencies
    required.packages <- c("ggplot2","reshape")
    missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
    if(length(missing.packages)) install.packages(missing.packages)
  library (ggplot2)
  library (reshape)
  
  
# Load data files
  ClinicalData.table <- read.csv ("./3 ANALISYS/CLINICAL DATA/Agilent_subset_clinicaldata.csv", header = TRUE)
  row.names(ClinicalData.table) <- ClinicalData.table$bcr_patient_barcode  
  ClinicalData.table$bcr_patient_barcode  <- NULL 

#analise Age distribution  
  Age <- as.data.frame(table( ClinicalData.table$age_at_diagnosis))
  row.names(Age) <- Age[,1] 
  Age[,1] <- NULL
  Age$Group <- cut(as.numeric(rownames(Age)), 
                         breaks = c(-Inf, 25, 35, 45, 65, Inf), 
                         labels = c("25", "26-35", "36-45", "46-65", "65+"), 
                         right = TRUE)
  Age.AG<- aggregate(Age["Freq"],Age["Group"],FUN=sum)
  rownames(Age.AG) <- Age.AG$Group
  Age.AG[,1] <- NULL
  write.csv (Age.AG,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_Age.csv")
  png("./FIGURES/clin_distr_MA_Age.png")
    par(mar=c(6,4,4,2))
    barplot( Age.AG$Freq, main="Age distribution", 
        xlab="Age Groups",ylab="Frequency",
        names.arg=row.names(Age.AG))
    Age.AG.graph <-  barplot( Age.AG$Freq, main="Age distribution", 
        xlab="Age Groups",ylab="Frequency",
        names.arg=row.names(Age.AG)) 
    text(x= Age.AG.graph, y= Age.AG$Freq+12, labels=as.character(Age.AG$Freq), xpd=TRUE)
  dev.off()
  
#analise Histology distribution
  Histotype <- as.data.frame(table( ClinicalData.table$histological_type ))
  row.names(Histotype) <- Histotype[,1] 
  Histotype[,1] <- NULL
  print (Histotype)
  write.csv (Histotype,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_Histotype.csv")
  png("./FIGURES/clin_distr_MA_Histotype.png")
    par(mar=c(12,4,4,2))
    barplot( Histotype$Freq, main="Histology distribution", 
           xlab=NULL,ylab="Frequency",
           names.arg=row.names(Histotype),
           cex.names=0.8,las=2)  
    Histotype.graph <- barplot( Histotype$Freq, main="Histology distribution", 
         xlab=NULL,ylab="Frequency",
         names.arg=row.names(Histotype),
         cex.names=0.8,las=2)
  text(x= Histotype.graph, y= Histotype$Freq+20, labels=as.character(Histotype$Freq), xpd=TRUE)
  dev.off()
  
#analise Vital Status distribution
  VS <- as.data.frame(table( ClinicalData.table$vital_status))
  row.names(VS) <- VS[,1] 
  VS[,1] <- NULL
  print (VS)  
  write.csv (VS,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_VS.csv")
  png("./FIGURES/clin_distr_MA_VS.png")
    par(mar=c(6,4,4,2))
    barplot( VS$Freq, main="Vital Status distribution", 
           xlab="Dead or Alive",ylab="Frequency",
           names.arg=row.names(VS))  
    VS.graph <- barplot( VS$Freq, main="Vital Status distribution", 
                xlab="Dead or Alive",ylab="Frequency",
                names.arg=row.names(VS))
    text(x= VS.graph, y= VS$Freq+12, labels=as.character(VS$Freq), xpd=TRUE)
  dev.off() 
  
  #analise Gender distribution
  Gender <- as.data.frame(table( ClinicalData.table$gender))
  row.names(Gender) <- Gender[,1] 
  Gender[,1] <- NULL
  print (Gender)  
  write.csv (Gender,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_Gender.csv")
  png("./FIGURES/clin_distr_MA_Gender.png")
  par(mar=c(6,4,4,2))
  barplot( Gender$Freq, main="Gender distribution", 
           xlab="Gender",ylab="Frequency",
           names.arg=row.names(Gender))  
  Gender.graph <- barplot( Gender$Freq, main="Gender distribution", 
                  xlab="Gender",ylab="Frequency",
                  names.arg=row.names(Gender))
  text(x= Gender.graph, y= Gender$Freq+12, labels=as.character(Gender$Freq), xpd=TRUE)
  dev.off() 
  
  #analise Neo-Adjuvant Treatment distribution
  NATreatment <- as.data.frame(table( ClinicalData.table$history_neoadjuvant_treatment))
  row.names(NATreatment) <- NATreatment[,1] 
  NATreatment[,1] <- NULL
  print (NATreatment)  
  write.csv (NATreatment,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_NATreatment.csv")
  png("./FIGURES/clin_distr_MA_NATreatment.png")
  par(mar=c(6,4,4,2))
  barplot( NATreatment$Freq, main="Neo-Adjuvant Treatment distribution", 
           xlab=NULL,ylab="Frequency",
           names.arg=row.names(NATreatment))  
  NATreatment.graph <- barplot( NATreatment$Freq, main="Neo-Adjuvant Treatment distribution", 
                       xlab=NULL,ylab="Frequency",
                       names.arg=row.names(NATreatment))
  text(x= NATreatment.graph, y= NATreatment$Freq+12, labels=as.character(NATreatment$Freq), xpd=TRUE)
  dev.off()

  #analise Patholigical tumor stage distribution
  pT <- as.data.frame(table( ClinicalData.table$ajcc_pathologic_tumor_stage))
  row.names(pT) <- pT[,1] 
  pT[,1] <- NULL
  print (pT)  
  write.csv (pT,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_pT.csv")
  png("./FIGURES/clin_distr_MA_pT.png")
  par(mar=c(6,4,4,2))
  barplot( pT$Freq, main="Patholigical tumor stage distribution", 
           xlab=NULL,ylab="Frequency",
           names.arg=row.names(pT))  
  pT.graph <- barplot( pT$Freq, main="Patholigical tumor stage distribution", 
                       xlab=NULL,ylab="Frequency",
                       names.arg=row.names(pT),
                       cex.names=0.8,las=2)
  text(x= pT.graph, y= pT$Freq+8, labels=as.character(pT$Freq), xpd=TRUE)
  dev.off()

  #analise IMS distribution prediction vs TCGA call
  IMS1 <- as.data.frame(table( ClinicalData.table$pam50.prediction.subtype))
  IMS2 <- as.data.frame(table( ClinicalData.table$Sorlie.prediction.subtype))
  IMS3 <- as.data.frame(table( ClinicalData.table$Hu.prediction.subtype))
  IMS4 <- as.data.frame(table( ClinicalData.table$TCGA.PAM50.RMethod))
  IMS5 <- as.data.frame(table( ClinicalData.table$PAM50))
  row.names(IMS1) <- IMS1[,1] 
  IMS1[,1] <- NULL
  row.names(IMS2) <- IMS2[,1] 
  IMS2[,1] <- NULL
  row.names(IMS3) <- IMS3[,1] 
  IMS3[,1] <- NULL
  row.names(IMS4) <- IMS4[,1] 
  IMS4[,1] <- NULL
  row.names(IMS5) <- IMS5[,1] 
  IMS5[,1] <- NULL
  IMS <- cbind (IMS1,IMS2,IMS3,IMS4,IMS5)
  colnames (IMS) <- c("Genefu.PAM50","Genefu.Sorlie","Genefu.Hu","TCGA.Pam50.RMethod","TCGA.Pam50.Call")
  print (IMS)  
  write.csv (IMS,"./DATA/Clinical Information/BRCA/Distribution/clin_distr_MA_IMS.csv")
  
  #combined IMS plot
  IMS.m <- melt(as.matrix(IMS),
                condition.vars=c("PAM50.prediction","Sorlie.prediction","Hu.prediction",,"PAM50.TCGA.Method","PAM50.TCGA.Call"),
                variable.name=c("Subtype","IMS.source"),
                value.name="Frequency");
  colnames(IMS.m) <- c("Subtype","IMS.source","Frequency")
  IMS.m.plot <- ggplot (data=IMS.m,aes(x=IMS.source,y=Frequency,fill=Subtype,))
  IMS.m.plot + geom_bar(stat="identity",position=position_dodge()) + ggtitle("Intrinsic subtype distribution")
  ggsave(file="./FIGURES/clin_distr_MA_IMS.png")
  
  #Single IMS prediction plot (TCGA)
  png("./FIGURES/clin_distr_MA_IMS.TCGA.png")
  par(mar=c(6,4,4,2))
  barplot( IMS$TCGA.Pam50.Call, main="Intrinsic Molecular subtype distribution", 
           xlab=NULL,ylab="Frequency",
           names.arg=row.names(IMS),
           cex.names=0.7,
           beside = TRUE)  
  IMS.graph <- barplot( IMS$TCGA.Pam50.Call, main="Intrinsic Molecular subtype distribution", 
                        xlab=NULL,ylab="Frequency",
                        names.arg=row.names(IMS),
                        cex.names=0.7) 
  text(x= IMS.graph, y= IMS$TCGA.Pam50.Call+12, labels=as.character(IMS$TCGA.Pam50.Call), xpd=TRUE)
  dev.off()
