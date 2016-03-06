# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
# Dependencies
required.packages <- c("mygene")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)

library(mygene)

# Load data
load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/BRCA.BSF2.subtypes separate/DEGinICR1vs4_Luminal AB.RDATA")
DEGinICR1vs4.Luminal <- DEGinICR1vs4
rm(DEGinICR1vs4)
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGinICR4vs1.RDATA" )
DEGinICR4vs1$logFC <- -(DEGinICR4vs1$logFC)
DEGinICR1vs4 <- DEGinICR4vs1 
rm(DEGinICR4vs1)
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/Cluster/DEGsLM_ICR4vs1.rdata")
rm(DEmerge)
DEGsLM_ICR4vs1$dm <- -(DEGsLM_ICR4vs1$dm)
DEGsLM_ICR1vs4 <- DEGsLM_ICR4vs1 
rm(DEGsLM_ICR4vs1)
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsall.rdata")
DEGinMAPKMUT<-DEGs
rm(DEGs)
load("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.rdata")
DEGinMAPKMUT.Luminal<-DEGs
rm(DEGs)

#master gene.list
master.gene.list<-unique(c(rownames(DEGinMAPKMUT),rownames(DEGinMAPKMUT.Luminal),rownames(DEGinICR1vs4),rownames(DEGinICR1vs4.Luminal),DEGsLM_ICR1vs4$Symbol))
gene.query <- queryMany(master.gene.list, scopes="symbol", fields=c("entrezgene", "go"), species="human")
gene.table <- as.data.frame(gene.query)

#append entrezID to datasets
DEGinICR1vs4$Entrez         <- gene.table$entrezgene[match(rownames(DEGinICR1vs4),gene.table$query)]
DEGinICR1vs4.Luminal$Entrez <- gene.table$entrezgene[match(rownames(DEGinICR1vs4.Luminal),gene.table$query)]
DEGinMAPKMUT$Entrez         <- gene.table$entrezgene[match(rownames(DEGinMAPKMUT),gene.table$query)]
DEGinMAPKMUT.Luminal$Entrez <- gene.table$entrezgene[match(rownames(DEGinMAPKMUT.Luminal),gene.table$query)]
DEGsLM_ICR1vs4$Entrez       <- gene.table$entrezgene[match(DEGsLM_ICR1vs4$Symbol,gene.table$query)]

#Create Entrez lists
CutFDR = 0.05
CutLogFC = 0.5

DEGinICR1vs4.UP <- DEGinICR1vs4[DEGinICR1vs4$FDR<CutFDR & DEGinICR1vs4$logFC>CutLogFC & !is.na(DEGinICR1vs4$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinICR1vs4.UP <- DEGinICR1vs4.UP[order(DEGinICR1vs4.UP$FDR),]
if (nrow(DEGinICR1vs4.UP) > 3000) {DEGinICR1vs4.UP <- DEGinICR1vs4.UP[1:3000,1,drop=FALSE]}
DEGinICR1vs4.DOWN <- DEGinICR1vs4[DEGinICR1vs4$FDR<CutFDR & DEGinICR1vs4$logFC<(-CutLogFC) & !is.na(DEGinICR1vs4$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinICR1vs4.DOWN <- DEGinICR1vs4.DOWN[order(DEGinICR1vs4.DOWN$FDR),]
if (nrow(DEGinICR1vs4.DOWN) > 3000) {DEGinICR1vs4.DOWN <- DEGinICR1vs4.DOWN[1:3000,1,drop=FALSE]}
DEGinICR1vs4.Luminal.UP <- DEGinICR1vs4.Luminal[DEGinICR1vs4.Luminal$FDR<CutFDR & DEGinICR1vs4.Luminal$logFC>CutLogFC & !is.na(DEGinICR1vs4.Luminal$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinICR1vs4.Luminal.UP <- DEGinICR1vs4.Luminal.UP[order(DEGinICR1vs4.Luminal.UP$FDR),]
if (nrow(DEGinICR1vs4.Luminal.UP) > 3000) {DEGinICR1vs4.Luminal.UP <- DEGinICR1vs4.Luminal.UP[1:3000,1,drop=FALSE]}
DEGinICR1vs4.Luminal.DOWN <- DEGinICR1vs4.Luminal[DEGinICR1vs4.Luminal$FDR<CutFDR & DEGinICR1vs4.Luminal$logFC<(-CutLogFC) & !is.na(DEGinICR1vs4.Luminal$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinICR1vs4.Luminal.DOWN <- DEGinICR1vs4.Luminal.DOWN[order(DEGinICR1vs4.Luminal.DOWN$FDR),]
if (nrow(DEGinICR1vs4.Luminal.DOWN) > 3000) {DEGinICR1vs4.Luminal.DOWN <- DEGinICR1vs4.Luminal.DOWN[1:3000,1,drop=FALSE]}
DEGinMAPKMUT.UP <- DEGinMAPKMUT[DEGinMAPKMUT$FDR<CutFDR & DEGinMAPKMUT$logFC>CutLogFC & !is.na(DEGinMAPKMUT$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinMAPKMUT.UP <- DEGinMAPKMUT.UP[order(DEGinMAPKMUT.UP$FDR),]
if (nrow(DEGinMAPKMUT.UP) > 3000) {DEGinMAPKMUT.UP <- DEGinMAPKMUT.UP[1:3000,1,drop=FALSE]}
DEGinMAPKMUT.DOWN <- DEGinMAPKMUT[DEGinMAPKMUT$FDR<CutFDR & DEGinMAPKMUT$logFC<(-CutLogFC) & !is.na(DEGinMAPKMUT$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinMAPKMUT.DOWN <- DEGinMAPKMUT.DOWN[order(DEGinMAPKMUT.DOWN$FDR),]
if (nrow(DEGinMAPKMUT.DOWN) > 3000) {DEGinMAPKMUT.DOWN <- DEGinMAPKMUT.DOWN[1:3000,1,drop=FALSE]}
DEGinMAPKMUT.Luminal.UP <- DEGinMAPKMUT.Luminal[DEGinMAPKMUT.Luminal$FDR<CutFDR & DEGinMAPKMUT.Luminal$logFC>CutLogFC & !is.na(DEGinMAPKMUT.Luminal$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinMAPKMUT.Luminal.UP <- DEGinMAPKMUT.Luminal.UP[order(DEGinMAPKMUT.Luminal.UP$FDR),]
if (nrow(DEGinMAPKMUT.Luminal.UP) > 3000) {DEGinMAPKMUT.Luminal.UP <- DEGinMAPKMUT.Luminal.UP[1:3000,1,drop=FALSE]}
DEGinMAPKMUT.Luminal.DOWN <- DEGinMAPKMUT.Luminal[DEGinMAPKMUT.Luminal$FDR<CutFDR & DEGinMAPKMUT.Luminal$logFC<(-CutLogFC) & !is.na(DEGinMAPKMUT.Luminal$Entrez) ,c("Entrez","FDR"),drop=FALSE]
DEGinMAPKMUT.Luminal.DOWN <- DEGinMAPKMUT.Luminal.DOWN[order(DEGinMAPKMUT.Luminal.DOWN$FDR),]
if (nrow(DEGinMAPKMUT.Luminal.DOWN) > 3000) {DEGinMAPKMUT.Luminal.DOWN <- DEGinMAPKMUT.Luminal.DOWN[1:3000,1,drop=FALSE]}
DEGsLM_ICR1vs4.UP <- DEGsLM_ICR1vs4[DEGsLM_ICR1vs4$FDR<CutFDR & DEGsLM_ICR1vs4$dm>CutLogFC & !is.na(DEGsLM_ICR1vs4$Entrez) ,c("Entrez","FDR","Symbol"),drop=FALSE]
DEGsLM_ICR1vs4.DOWN <- DEGsLM_ICR1vs4[DEGsLM_ICR1vs4$FDR<CutFDR & DEGsLM_ICR1vs4$dm<(-CutLogFC) & !is.na(DEGsLM_ICR1vs4$Entrez) ,c("Entrez","FDR","Symbol"),drop=FALSE]
rownames(DEGsLM_ICR1vs4.UP)<-NULL
rownames(DEGsLM_ICR1vs4.DOWN)<-NULL
DEGsLM_ICR1vs4.UP <- DEGsLM_ICR1vs4.UP[order(DEGsLM_ICR1vs4.UP$FDR),]
DEGsLM_ICR1vs4.DOWN <- DEGsLM_ICR1vs4.DOWN[order(DEGsLM_ICR1vs4.DOWN$FDR),]
DEGsLM_ICR1vs4.UP$FDR <- NULL
DEGsLM_ICR1vs4.DOWN$FDR <- NULL
DEGsLM_ICR1vs4.UP<-unique(DEGsLM_ICR1vs4.UP)
DEGsLM_ICR1vs4.DOWN<-unique(DEGsLM_ICR1vs4.DOWN)
rownames(DEGsLM_ICR1vs4.UP)<-DEGsLM_ICR1vs4.UP$Symbol
rownames(DEGsLM_ICR1vs4.DOWN)<-DEGsLM_ICR1vs4.DOWN$Symbol
DEGsLM_ICR1vs4.UP$Symbol <- NULL
DEGsLM_ICR1vs4.DOWN$Symbol <- NULL
if (nrow(DEGsLM_ICR1vs4.UP) > 3000) {DEGsLM_ICR1vs4.UP <- DEGsLM_ICR1vs4.UP[1:3000,1,drop=FALSE]}
if (nrow(DEGsLM_ICR1vs4.DOWN) > 3000) {DEGsLM_ICR1vs4.DOWN <- DEGsLM_ICR1vs4.DOWN[1:3000,1,drop=FALSE]}

DEGinICR1vs4.UP$FDR <- NULL
DEGinICR1vs4.DOWN$FDR <- NULL
DEGinICR1vs4.Luminal.UP$FDR <- NULL
DEGinICR1vs4.Luminal.DOWN$FDR <- NULL
DEGinMAPKMUT.UP$FDR <- NULL
DEGinMAPKMUT.DOWN$FDR <- NULL
DEGinMAPKMUT.Luminal.UP$FDR <- NULL
DEGinMAPKMUT.Luminal.DOWN$FDR <- NULL
DEGsLM_ICR1vs4.UP$FDR <- NULL
DEGsLM_ICR1vs4.DOWN$FDR <- NULL

write.csv (DEGinICR1vs4.UP,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/FORDAVID/DEGinICR1vs4.UP.csv")

DEGinICR1vs4.BOTH <- DEGinICR1vs4[DEGinICR1vs4$FDR<CutFDR & abs(DEGinICR1vs4$logFC)>CutLogFC & !is.na(DEGinICR1vs4$Entrez) ,c("Entrez","FDR"),drop=FALSE]
class(DEGinICR1vs4.BOTH$FDR) = "numeric"
DEGinICR1vs4.BOTH <- DEGinICR1vs4.BOTH[order(DEGinICR1vs4.BOTH$FDR),]
if (nrow(DEGinICR1vs4.BOTH) > 3000) {DEGinICR1vs4.BOTH <- DEGinICR1vs4.BOTH[1:3000,1,drop=FALSE]}

DEGinICR1vs4.Luminal.BOTH <- DEGinICR1vs4.Luminal[DEGinICR1vs4.Luminal$FDR<CutFDR & abs(DEGinICR1vs4.Luminal$logFC)>CutLogFC & !is.na(DEGinICR1vs4.Luminal$Entrez) ,c("Entrez","FDR"),drop=FALSE]
class(DEGinICR1vs4.Luminal.BOTH$FDR) = "numeric"
DEGinICR1vs4.Luminal.BOTH <- DEGinICR1vs4.Luminal.BOTH[order(DEGinICR1vs4.Luminal.BOTH$FDR),]
if (nrow(DEGinICR1vs4.Luminal.BOTH) > 3000) {DEGinICR1vs4.Luminal.BOTH <- DEGinICR1vs4.Luminal.BOTH[1:3000,1,drop=FALSE]}

DEGinMAPKMUT.BOTH <-DEGinMAPKMUT[DEGinMAPKMUT$FDR<CutFDR & abs(DEGinMAPKMUT$logFC)>CutLogFC & !is.na(DEGinMAPKMUT$Entrez) ,c("Entrez","FDR"),drop=FALSE]
class(DEGinMAPKMUT.BOTH$FDR) = "numeric"
DEGinMAPKMUT.BOTH <- DEGinMAPKMUT.BOTH[order(DEGinMAPKMUT.BOTH$FDR),]
if (nrow(DEGinMAPKMUT.BOTH) > 3000) {DEGinMAPKMUT.BOTH <- DEGinMAPKMUT.BOTH[1:3000,1,drop=FALSE]}

DEGinMAPKMUT.Luminal.BOTH <-DEGinMAPKMUT.Luminal[DEGinMAPKMUT.Luminal$FDR<CutFDR & abs(DEGinMAPKMUT.Luminal$logFC)>CutLogFC & !is.na(DEGinMAPKMUT.Luminal$Entrez) ,c("Entrez","FDR"),drop=FALSE]
class(DEGinMAPKMUT.Luminal.BOTH$FDR) = "numeric"
DEGinMAPKMUT.Luminal.BOTH <- DEGinMAPKMUT.Luminal.BOTH[order(DEGinMAPKMUT.Luminal.BOTH$FDR),]
if (nrow(DEGinMAPKMUT.Luminal.BOTH) > 3000) {DEGinMAPKMUT.Luminal.BOTH <- DEGinMAPKMUT.Luminal.BOTH[1:3000,1,drop=FALSE]}

DEGsLM_ICR1vs4.BOTH <-DEGsLM_ICR1vs4[DEGsLM_ICR1vs4$FDR<CutFDR & abs(DEGsLM_ICR1vs4$dm)>CutLogFC & !is.na(DEGsLM_ICR1vs4$Entrez) ,c("Entrez","FDR","Symbol"),drop=FALSE]
class(DEGsLM_ICR1vs4.BOTH$FDR) = "numeric"
DEGsLM_ICR1vs4.BOTH <- DEGsLM_ICR1vs4.BOTH[order(DEGsLM_ICR1vs4.BOTH$FDR),]
rownames(DEGsLM_ICR1vs4.BOTH) <- NULL
DEGsLM_ICR1vs4.BOTH$FDR <- NULL
DEGsLM_ICR1vs4.BOTH <- unique(DEGsLM_ICR1vs4.BOTH)
rownames(DEGsLM_ICR1vs4.BOTH)<-DEGsLM_ICR1vs4.BOTH$Symbol
DEGsLM_ICR1vs4.BOTH$Symbol <-NULL
if (nrow(DEGsLM_ICR1vs4.BOTH) > 3000) {DEGsLM_ICR1vs4.BOTH <- DEGsLM_ICR1vs4.BOTH[1:3000,1,drop=FALSE]}


#DAVID
#source("https://bioconductor.org/biocLite.R")
#biocLite("RDAVIDWebService")
library("RDAVIDWebService")
#install.packages('rJava', .libPaths()[1], 'http://www.rforge.net/')
library(rJava)
david<-DAVIDWebService$new(email="whendrickx@sidra.org",url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
#TCGA DATA
#ICR 1 vs ICR4 
result<-addList(david, DEGinICR1vs4.UP,idType="ENTREZ_GENE_ID",listName="DEGinICR1vs4.UP", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinICR1vs4.UP.KEGG <- getFunctionalAnnotationChart(david)
result<-addList(david, DEGinICR1vs4.DOWN,idType="ENTREZ_GENE_ID",listName="DEGinICR1vs4.DOWN", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinICR1vs4.DOWN.KEGG <- getFunctionalAnnotationChart(david)
#ICR 1 vs ICR4 LUMINAL
result<-addList(david, DEGinICR1vs4.Luminal.UP,idType="ENTREZ_GENE_ID",listName="DEGinICR1vs4.Luminal.UP", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinICR1vs4.Luminal.UP.KEGG <- getFunctionalAnnotationChart(david)
result<-addList(david, DEGinICR1vs4.Luminal.DOWN,idType="ENTREZ_GENE_ID",listName="DEGinICR1vs4.Luminal.DOWN", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinICR1vs4.Luminal.DOWN.KEGG <- getFunctionalAnnotationChart(david)
#MUT vs WT
result<-addList(david, DEGinMAPKMUT.UP,idType="ENTREZ_GENE_ID",listName="DEGinMAPKMUT.UP", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinMAPKMUT.UP.KEGG <- getFunctionalAnnotationChart(david)
result<-addList(david, DEGinMAPKMUT.DOWN,idType="ENTREZ_GENE_ID",listName="DEGinMAPKMUT.DOWN", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinMAPKMUT.DOWN.KEGG <- getFunctionalAnnotationChart(david)
#MUT vs WT LUMINAL
result<-addList(david, DEGinMAPKMUT.Luminal.UP,idType="ENTREZ_GENE_ID",listName="DEGinMAPKMUT.Luminal.UP", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinMAPKMUT.Luminal.UP.KEGG <- getFunctionalAnnotationChart(david)
result<-addList(david, DEGinMAPKMUT.Luminal.DOWN,idType="ENTREZ_GENE_ID",listName="DEGinMAPKMUT.Luminal.DOWN", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinMAPKMUT.Luminal.DOWN.KEGG <- getFunctionalAnnotationChart(david)
#LANCE MILLER DATA
#ICR 1 vs ICR4
result<-addList(david, DEGsLM_ICR1vs4.UP,idType="ENTREZ_GENE_ID",listName="DEGsLM_ICR1vs4.UP", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGsLM_ICR1vs4.UP.KEGG <- getFunctionalAnnotationChart(david)
result<-addList(david, DEGsLM_ICR1vs4.DOWN,idType="ENTREZ_GENE_ID",listName="DEGsLM_ICR1vs4.DOWN", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGsLM_ICR1vs4.DOWN.KEGG <- getFunctionalAnnotationChart(david)

#combine UP and DOWN for TCGA ICR1 vs ICR4
result<-addList(david,DEGinICR1vs4.BOTH,idType="ENTREZ_GENE_ID",listName="DEGinICR1vs4.BOTH", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinICR1vs4.BOTH.KEGG <- getFunctionalAnnotationChart(david)
#combine UP and DOWN for TCGA ICR1 vs ICR4 Luminal
result<-addList(david,DEGinICR1vs4.Luminal.BOTH,idType="ENTREZ_GENE_ID",listName="DEGinICR1vs4.Luminal.BOTH", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinICR1vs4.Luminal.BOTH.KEGG <- getFunctionalAnnotationChart(david)
#combine UP and DOWN for TCGA MUT vs WT
result<-addList(david,DEGinMAPKMUT.BOTH,idType="ENTREZ_GENE_ID",listName="DEGinMAPKMUT.BOTH", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinMAPKMUT.BOTH.KEGG <- getFunctionalAnnotationChart(david)
#combine UP and DOWN for TCGA MUT vs WT Luminal
result<-addList(david,DEGinMAPKMUT.Luminal.BOTH,idType="ENTREZ_GENE_ID",listName="DEGinMAPKMUT.Luminal.BOTH", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGinMAPKMUT.Luminal.BOTH.KEGG <- getFunctionalAnnotationChart(david)
#combine UP and DOWN for LM dataset ICR1 vs ICR4
result<-addList(david,DEGsLM_ICR1vs4.BOTH,idType="ENTREZ_GENE_ID",listName="DEGsLM_ICR1vs4.BOTH", listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY"))
DEGsLM_ICR1vs4.BOTH.KEGG <- getFunctionalAnnotationChart(david)


#Save KEGG
save (DEGinICR1vs4.UP.KEGG,DEGinICR1vs4.UP.KEGG,
      DEGinICR1vs4.Luminal.UP.KEGG,DEGinICR1vs4.Luminal.DOWN.KEGG,
      DEGinMAPKMUT.UP.KEGG,DEGinMAPKMUT.UP.KEGG,
      DEGinMAPKMUT.Luminal.UP.KEGG,DEGinMAPKMUT.Luminal.DOWN.KEGG,
      DEGsLM_ICR1vs4.UP.KEGG,DEGsLM_ICR1vs4.DOWN.KEGG,
      DEGinICR1vs4.BOTH.KEGG,
      DEGinICR1vs4.Luminal.BOTH.KEGG,
      DEGinMAPKMUT.BOTH.KEGG,
      DEGinMAPKMUT.Luminal.BOTH.KEGG,
      DEGsLM_ICR1vs4.BOTH.KEGG,
      file= "./3 ANALISYS/DIFFERENTIAL EXPRESSION/FORDAVID/R-david-results.FC_cut_0.0.Rdata")
