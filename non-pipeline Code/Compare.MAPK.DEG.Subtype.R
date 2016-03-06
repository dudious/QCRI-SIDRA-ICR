# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

source ("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/heatmap.3.R")

#load data
MAPK.PW.genes <- read.csv ("./3 ANALISYS/MAPK/MAPK.pathway.genes.csv",header = FALSE)

load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGinICR4vs1.RDATA")
DEGinICR1vs4.All <- DEGinICR4vs1
DEGinICR1vs4.All$logFC <- -(DEGinICR1vs4.All$logFC) #reversed for ICR1vs4 instaed of 1vs4
rm(DEGinICR4vs1)
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.rdata")
DEGinMUTvsWT.Lum <-DEGs
rm(DEGs)
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/BRCA.BSF2.subtypes separate/DEGinICR1vs4_Basal-like.RDATA")
DEGinICR1vs4.Bas <- DEGinICR1vs4
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/BRCA.BSF2.subtypes separate/DEGinICR1vs4_HER2-enriched.RDATA")
DEGinICR1vs4.Her <- DEGinICR1vs4
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/BRCA.BSF2.subtypes separate/DEGinICR1vs4_Luminal A.RDATA")
DEGinICR1vs4.LmA <- DEGinICR1vs4
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/BRCA.BSF2.subtypes separate/DEGinICR1vs4_Luminal B.RDATA")
DEGinICR1vs4.LmB <- DEGinICR1vs4
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/BRCA.BSF2.subtypes separate/DEGinICR1vs4_Luminal AB.RDATA")
DEGinICR1vs4.Lum <- DEGinICR1vs4
rm(DEGinICR1vs4)
#load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGsLM_ICR4vs1.rdata")


#select significant MAPK.PW genes in TCGA ICR DEG  or MAPK DEG

#### ICR1
##### All subtypes
###### UPREGULATED
DEGinICR1vs4.All.MAPK.PW.UP<-DEGinICR1vs4.All[rownames(DEGinICR1vs4.All)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.All.MAPK.PW.UP<-DEGinICR1vs4.All.MAPK.PW.UP[DEGinICR1vs4.All.MAPK.PW.UP$logFC>0,]   
DEGinICR1vs4.All.MAPK.PW.UP<-DEGinICR1vs4.All.MAPK.PW.UP[order(DEGinICR1vs4.All.MAPK.PW.UP$FDR),]
DEGinICR1vs4.All.MAPK.PW.UP<-DEGinICR1vs4.All.MAPK.PW.UP[DEGinICR1vs4.All.MAPK.PW.UP$FDR<0.05,]
###### DOWNREGULATED
DEGinICR1vs4.All.MAPK.PW.DOWN<-DEGinICR1vs4.All[rownames(DEGinICR1vs4.All)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.All.MAPK.PW.DOWN<-DEGinICR1vs4.All.MAPK.PW.DOWN[DEGinICR1vs4.All.MAPK.PW.DOWN$logFC<0,]   
DEGinICR1vs4.All.MAPK.PW.DOWN<-DEGinICR1vs4.All.MAPK.PW.DOWN[order(DEGinICR1vs4.All.MAPK.PW.DOWN$FDR),]
DEGinICR1vs4.All.MAPK.PW.DOWN<-DEGinICR1vs4.All.MAPK.PW.DOWN[DEGinICR1vs4.All.MAPK.PW.DOWN$FDR<0.05,]

##### Luminal A
###### UPREGULATED
DEGinICR1vs4.LmA.MAPK.PW.UP<-DEGinICR1vs4.LmA[rownames(DEGinICR1vs4.LmA)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.LmA.MAPK.PW.UP<-DEGinICR1vs4.LmA.MAPK.PW.UP[DEGinICR1vs4.LmA.MAPK.PW.UP$logFC>0,]   
DEGinICR1vs4.LmA.MAPK.PW.UP<-DEGinICR1vs4.LmA.MAPK.PW.UP[order(DEGinICR1vs4.LmA.MAPK.PW.UP$FDR),]
DEGinICR1vs4.LmA.MAPK.PW.UP<-DEGinICR1vs4.LmA.MAPK.PW.UP[DEGinICR1vs4.LmA.MAPK.PW.UP$FDR<0.05,]
###### DOWNREGULATED
DEGinICR1vs4.LmA.MAPK.PW.DOWN<-DEGinICR1vs4.LmA[rownames(DEGinICR1vs4.LmA)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.LmA.MAPK.PW.DOWN<-DEGinICR1vs4.LmA.MAPK.PW.DOWN[DEGinICR1vs4.LmA.MAPK.PW.DOWN$logFC<0,]   
DEGinICR1vs4.LmA.MAPK.PW.DOWN<-DEGinICR1vs4.LmA.MAPK.PW.DOWN[order(DEGinICR1vs4.LmA.MAPK.PW.DOWN$FDR),]
DEGinICR1vs4.LmA.MAPK.PW.DOWN<-DEGinICR1vs4.LmA.MAPK.PW.DOWN[DEGinICR1vs4.LmA.MAPK.PW.DOWN$FDR<0.05,]

##### Luminal B
###### UPREGULATED
DEGinICR1vs4.LmB.MAPK.PW.UP<-DEGinICR1vs4.LmB[rownames(DEGinICR1vs4.LmB)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.LmB.MAPK.PW.UP<-DEGinICR1vs4.LmB.MAPK.PW.UP[DEGinICR1vs4.LmB.MAPK.PW.UP$logFC>0,]   
DEGinICR1vs4.LmB.MAPK.PW.UP<-DEGinICR1vs4.LmB.MAPK.PW.UP[order(DEGinICR1vs4.LmB.MAPK.PW.UP$FDR),]
DEGinICR1vs4.LmB.MAPK.PW.UP<-DEGinICR1vs4.LmB.MAPK.PW.UP[DEGinICR1vs4.LmB.MAPK.PW.UP$FDR<0.05,]
###### DOWNREGULATED
DEGinICR1vs4.LmB.MAPK.PW.DOWN<-DEGinICR1vs4.LmB[rownames(DEGinICR1vs4.LmB)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.LmB.MAPK.PW.DOWN<-DEGinICR1vs4.LmB.MAPK.PW.DOWN[DEGinICR1vs4.LmB.MAPK.PW.DOWN$logFC<0,]   
DEGinICR1vs4.LmB.MAPK.PW.DOWN<-DEGinICR1vs4.LmB.MAPK.PW.DOWN[order(DEGinICR1vs4.LmB.MAPK.PW.DOWN$FDR),]
DEGinICR1vs4.LmB.MAPK.PW.DOWN<-DEGinICR1vs4.LmB.MAPK.PW.DOWN[DEGinICR1vs4.LmB.MAPK.PW.DOWN$FDR<0.05,]

##### Basal-like
###### UPREGULATED
DEGinICR1vs4.Bas.MAPK.PW.UP<-DEGinICR1vs4.Bas[rownames(DEGinICR1vs4.Bas)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.Bas.MAPK.PW.UP<-DEGinICR1vs4.Bas.MAPK.PW.UP[DEGinICR1vs4.Bas.MAPK.PW.UP$logFC>0,]   
DEGinICR1vs4.Bas.MAPK.PW.UP<-DEGinICR1vs4.Bas.MAPK.PW.UP[order(DEGinICR1vs4.Bas.MAPK.PW.UP$FDR),]
DEGinICR1vs4.Bas.MAPK.PW.UP<-DEGinICR1vs4.Bas.MAPK.PW.UP[DEGinICR1vs4.Bas.MAPK.PW.UP$FDR<0.05,]
###### DOWNREGULATED
DEGinICR1vs4.Bas.MAPK.PW.DOWN<-DEGinICR1vs4.Bas[rownames(DEGinICR1vs4.Bas)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.Bas.MAPK.PW.DOWN<-DEGinICR1vs4.Bas.MAPK.PW.DOWN[DEGinICR1vs4.Bas.MAPK.PW.DOWN$logFC<0,]   
DEGinICR1vs4.Bas.MAPK.PW.DOWN<-DEGinICR1vs4.Bas.MAPK.PW.DOWN[order(DEGinICR1vs4.Bas.MAPK.PW.DOWN$FDR),]
DEGinICR1vs4.Bas.MAPK.PW.DOWN<-DEGinICR1vs4.Bas.MAPK.PW.DOWN[DEGinICR1vs4.Bas.MAPK.PW.DOWN$FDR<0.05,]

##### HER2-Enriched
###### UPREGULATED
DEGinICR1vs4.Her.MAPK.PW.UP<-DEGinICR1vs4.Her[rownames(DEGinICR1vs4.Her)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.Her.MAPK.PW.UP<-DEGinICR1vs4.Her.MAPK.PW.UP[DEGinICR1vs4.Her.MAPK.PW.UP$logFC>0,]   
DEGinICR1vs4.Her.MAPK.PW.UP<-DEGinICR1vs4.Her.MAPK.PW.UP[order(DEGinICR1vs4.Her.MAPK.PW.UP$FDR),]
DEGinICR1vs4.Her.MAPK.PW.UP<-DEGinICR1vs4.Her.MAPK.PW.UP[DEGinICR1vs4.Her.MAPK.PW.UP$FDR<0.05,]
###### DOWNREGULATED
DEGinICR1vs4.Her.MAPK.PW.DOWN<-DEGinICR1vs4.Her[rownames(DEGinICR1vs4.Her)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.Her.MAPK.PW.DOWN<-DEGinICR1vs4.Her.MAPK.PW.DOWN[DEGinICR1vs4.Her.MAPK.PW.DOWN$logFC<0,]   
DEGinICR1vs4.Her.MAPK.PW.DOWN<-DEGinICR1vs4.Her.MAPK.PW.DOWN[order(DEGinICR1vs4.Her.MAPK.PW.DOWN$FDR),]
DEGinICR1vs4.Her.MAPK.PW.DOWN<-DEGinICR1vs4.Her.MAPK.PW.DOWN[DEGinICR1vs4.Her.MAPK.PW.DOWN$FDR<0.05,]

##### Luminal
###### UPREGULATED
DEGinICR1vs4.Lum.MAPK.PW.UP<-DEGinICR1vs4.Lum[rownames(DEGinICR1vs4.Lum)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.Lum.MAPK.PW.UP<-DEGinICR1vs4.Lum.MAPK.PW.UP[DEGinICR1vs4.Lum.MAPK.PW.UP$logFC>0,]   
DEGinICR1vs4.Lum.MAPK.PW.UP<-DEGinICR1vs4.Lum.MAPK.PW.UP[order(DEGinICR1vs4.Lum.MAPK.PW.UP$FDR),]
DEGinICR1vs4.Lum.MAPK.PW.UP<-DEGinICR1vs4.Lum.MAPK.PW.UP[DEGinICR1vs4.Lum.MAPK.PW.UP$FDR<0.05,]
###### DOWNREGULATED
DEGinICR1vs4.Lum.MAPK.PW.DOWN<-DEGinICR1vs4.Lum[rownames(DEGinICR1vs4.Lum)%in%MAPK.PW.genes$V2,]
DEGinICR1vs4.Lum.MAPK.PW.DOWN<-DEGinICR1vs4.Lum.MAPK.PW.DOWN[DEGinICR1vs4.Lum.MAPK.PW.DOWN$logFC<0,]   
DEGinICR1vs4.Lum.MAPK.PW.DOWN<-DEGinICR1vs4.Lum.MAPK.PW.DOWN[order(DEGinICR1vs4.Lum.MAPK.PW.DOWN$FDR),]
DEGinICR1vs4.Lum.MAPK.PW.DOWN<-DEGinICR1vs4.Lum.MAPK.PW.DOWN[DEGinICR1vs4.Lum.MAPK.PW.DOWN$FDR<0.05,]

#### Mutated
##### Luminal
###### UPREGULATED Luminal MUT
DEGinMUTvsWT.Lum.MAPK.PW.UP<-DEGinMUTvsWT.Lum[rownames(DEGinMUTvsWT.Lum)%in%MAPK.PW.genes$V2,]
DEGinMUTvsWT.Lum.MAPK.PW.UP<-DEGinMUTvsWT.Lum.MAPK.PW.UP[DEGinMUTvsWT.Lum.MAPK.PW.UP$logFC>0,]
DEGinMUTvsWT.Lum.MAPK.PW.UP<-DEGinMUTvsWT.Lum.MAPK.PW.UP[order(DEGinMUTvsWT.Lum.MAPK.PW.UP$FDR),]
DEGinMUTvsWT.Lum.MAPK.PW.UP<-DEGinMUTvsWT.Lum.MAPK.PW.UP[DEGinMUTvsWT.Lum.MAPK.PW.UP$FDR<0.05,]

###### DOWNREGULATED Luminal MUT
DEGinMUTvsWT.Lum.MAPK.PW.DOWN<-DEGinMUTvsWT.Lum[rownames(DEGinMUTvsWT.Lum)%in%MAPK.PW.genes$V2,]
DEGinMUTvsWT.Lum.MAPK.PW.DOWN<-DEGinMUTvsWT.Lum.MAPK.PW.DOWN[DEGinMUTvsWT.Lum.MAPK.PW.DOWN$logFC<0,]
DEGinMUTvsWT.Lum.MAPK.PW.DOWN<-DEGinMUTvsWT.Lum.MAPK.PW.DOWN[order(DEGinMUTvsWT.Lum.MAPK.PW.DOWN$FDR),]
DEGinMUTvsWT.Lum.MAPK.PW.DOWN<-DEGinMUTvsWT.Lum.MAPK.PW.DOWN[DEGinMUTvsWT.Lum.MAPK.PW.DOWN$FDR<0.05,]

## Overview table
MAPK.DEG.comparison <- data.frame (gene.symbol=MAPK.PW.genes$V2)
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.All.MAPK.PW.UP),"ICR1vs4.all"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.All.MAPK.PW.DOWN),"ICR1vs4.all"] <- "DOWN"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.LmA.MAPK.PW.UP),"ICR1vs4.LmA"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.LmA.MAPK.PW.DOWN),"ICR1vs4.LmA"] <- "DOWN"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.LmB.MAPK.PW.UP),"ICR1vs4.LmB"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.LmB.MAPK.PW.DOWN),"ICR1vs4.LmB"] <- "DOWN"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.Bas.MAPK.PW.UP),"ICR1vs4.Bas"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.Bas.MAPK.PW.DOWN),"ICR1vs4.Bas"] <- "DOWN"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.Her.MAPK.PW.UP),"ICR1vs4.Her"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.Her.MAPK.PW.DOWN),"ICR1vs4.Her"] <- "DOWN"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.Lum.MAPK.PW.UP),"ICR1vs4.Lum"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinICR1vs4.Lum.MAPK.PW.DOWN),"ICR1vs4.Lum"] <- "DOWN"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinMUTvsWT.Lum.MAPK.PW.UP),"MUTvsWT.Lum"] <- "UP"
MAPK.DEG.comparison[MAPK.DEG.comparison$gene.symbol %in% rownames(DEGinMUTvsWT.Lum.MAPK.PW.DOWN),"MUTvsWT.Lum"] <- "DOWN"

write.csv (MAPK.DEG.comparison,"./3 ANALISYS/MAPK/DEG.comparison.csv")
