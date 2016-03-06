# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

#load data
load ("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.Basal.DBGS3.FLTR.NonSilent.VariationTables.RData")
variation.table$X <- NULL
variation.table.basal <- variation.table
load ("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.LumB.DBGS3.FLTR.NonSilent.VariationTables.RData")
variation.table$X <- NULL
variation.table.LumB <- variation.table
load ("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.LumA.DBGS3.FLTR.NonSilent.VariationTables.RData")
variation.table$X <- NULL
variation.table.LumA <- variation.table
load ("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA
      .BSF2.Her2.DBGS3.FLTR.NonSilent.VariationTables.RData")
variation.table$X <- NULL
variation.table.Her2 <- variation.table
load ("./3 ANALISYS/Mutations/BRCA.BSF2/BRCA.BSF2.All.DBGS3.FLTR.NonSilent.VariationTables.RData")
variation.table$X <- NULL
variation.table.All <- variation.table

rm(list=ls()[1:7])

variation.table.LumA$ICR1min4 <- as.numeric(lapply(strsplit(as.character(variation.table.LumA$Cluster_Percentages),split= ":"),'[',1))-
  as.numeric(lapply(strsplit(as.character(variation.table.LumA$Cluster_Percentages),split= ":"),'[',4))
variation.table.LumB$ICR1min4 <- as.numeric(lapply(strsplit(as.character(variation.table.LumB$Cluster_Percentages),split= ":"),'[',1))-
  as.numeric(lapply(strsplit(as.character(variation.table.LumB$Cluster_Percentages),split= ":"),'[',4))
variation.table.basal$ICR1min4 <- as.numeric(lapply(strsplit(as.character(variation.table.basal$Cluster_Percentages),split= ":"),'[',1))-
  as.numeric(lapply(strsplit(as.character(variation.table.basal$Cluster_Percentages),split= ":"),'[',4))
variation.table.Her2$ICR1min4 <- as.numeric(lapply(strsplit(as.character(variation.table.Her2$Cluster_Percentages),split= ":"),'[',1))-
  as.numeric(lapply(strsplit(as.character(variation.table.Her2$Cluster_Percentages),split= ":"),'[',4))
variation.table.All$ICR1min4 <- as.numeric(lapply(strsplit(as.character(variation.table.All$Cluster_Percentages),split= ":"),'[',1))-
  as.numeric(lapply(strsplit(as.character(variation.table.All$Cluster_Percentages),split= ":"),'[',4))


nrow(variation.table.LumA[variation.table.LumA$ChiSquare_ICR1vs4<0.05 & variation.table.LumA$ICR1min4>0,])
nrow(variation.table.LumA[variation.table.LumA$ChiSquare_ICR1vs4<0.05 & variation.table.LumA$ICR1min4<0,])
nrow(variation.table.LumB[variation.table.LumB$ChiSquare_ICR1vs4<0.05 & variation.table.LumB$ICR1min4>0,])
nrow(variation.table.LumB[variation.table.LumB$ChiSquare_ICR1vs4<0.05 & variation.table.LumB$ICR1min4<0,])
nrow(variation.table.basal[variation.table.basal$ChiSquare_ICR1vs4 <0.05 & variation.table.basal$ICR1min4>0,])
nrow(variation.table.basal[variation.table.basal$ChiSquare_ICR1vs4 <0.05 & variation.table.basal$ICR1min4<0,])
nrow(variation.table.Her2[variation.table.Her2$ChiSquare_ICR1vs4<0.05 & variation.table.Her2$ICR1min4>0,])
nrow(variation.table.Her2[variation.table.Her2$ChiSquare_ICR1vs4<0.05 & variation.table.Her2$ICR1min4<0,])
nrow(variation.table.All[variation.table.All$ChiSquare_ICR1vs4<0.05 & variation.table.All$ICR1min4>0,])
nrow(variation.table.All[variation.table.All$ChiSquare_ICR1vs4<0.05 & variation.table.All$ICR1min4<0,])

nrow(variation.table.LumA[variation.table.LumA$db.test == "TRUE" & variation.table.LumA$ICR1min4>0,])
nrow(variation.table.LumA[variation.table.LumA$db.test == "TRUE" & variation.table.LumA$ICR1min4<0,])
nrow(variation.table.LumB[variation.table.LumB$db.test == "TRUE" & variation.table.LumB$ICR1min4>0,])
nrow(variation.table.LumB[variation.table.LumB$db.test == "TRUE" & variation.table.LumB$ICR1min4<0,])
nrow(variation.table.basal[variation.table.basal$db.test == "TRUE" & variation.table.basal$ICR1min4>0,])
nrow(variation.table.basal[variation.table.basal$db.test == "TRUE" & variation.table.basal$ICR1min4<0,])
nrow(variation.table.Her2[variation.table.Her2$db.test == "TRUE" & variation.table.Her2$ICR1min4>0,])
nrow(variation.table.Her2[variation.table.Her2$db.test == "TRUE" & variation.table.Her2$ICR1min4<0,])
nrow(variation.table.All[variation.table.All$db.test == "TRUE" & variation.table.All$ICR1min4>0,])
nrow(variation.table.All[variation.table.All$db.test == "TRUE" & variation.table.All$ICR1min4<0,])

nrow(variation.table.LumA[(variation.table.LumA$db.test == "TRUE" | variation.table.LumA$ChiSquare_ICR1vs4<0.05) & variation.table.LumA$ICR1min4>0,])
nrow(variation.table.LumA[(variation.table.LumA$db.test == "TRUE" | variation.table.LumA$ChiSquare_ICR1vs4<0.05) & variation.table.LumA$ICR1min4<0,])
nrow(variation.table.LumB[(variation.table.LumB$db.test == "TRUE" | variation.table.LumB$ChiSquare_ICR1vs4<0.05) & variation.table.LumB$ICR1min4>0,])
nrow(variation.table.LumB[(variation.table.LumB$db.test == "TRUE" | variation.table.LumB$ChiSquare_ICR1vs4<0.05) & variation.table.LumB$ICR1min4<0,])
nrow(variation.table.basal[(variation.table.basal$db.test == "TRUE" | variation.table.basal$ChiSquare_ICR1vs4<0.05) & variation.table.basal$ICR1min4>0,])
nrow(variation.table.basal[(variation.table.basal$db.test == "TRUE" | variation.table.basal$ChiSquare_ICR1vs4<0.05) & variation.table.basal$ICR1min4<0,])
nrow(variation.table.Her2[(variation.table.Her2$db.test == "TRUE" | variation.table.Her2$ChiSquare_ICR1vs4<0.05) & variation.table.Her2$ICR1min4>0,])
nrow(variation.table.Her2[(variation.table.Her2$db.test == "TRUE" | variation.table.Her2$ChiSquare_ICR1vs4<0.05) & variation.table.Her2$ICR1min4<0,])
nrow(variation.table.All[(variation.table.All$db.test == "TRUE" | variation.table.All$ChiSquare_ICR1vs4<0.05) & variation.table.All$ICR1min4>0,])
nrow(variation.table.All[(variation.table.All$db.test == "TRUE" | variation.table.All$ChiSquare_ICR1vs4<0.05) & variation.table.All$ICR1min4<0,])
