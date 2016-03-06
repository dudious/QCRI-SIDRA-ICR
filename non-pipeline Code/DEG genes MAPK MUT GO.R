# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

# Dependencies
source("https://bioconductor.org/biocLite.R")
required.packages <- c("ReactomePA","mygene","biomaRt","tmod","org.Hs.eg.db")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) biocLite(missing.packages)
#library (gage)
#library(gageData)
library(ReactomePA)
library(mygene)
library(org.Hs.eg.db)
library(biomaRt)
library(reshape)
library(limma)
library(tmod)

load ("./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.rdata")

gene.list.name <- rownames(DEGs)

gene.IDs <- queryMany(gene.list.name, scopes="symbol", fields=c("entrezgene", "go"), species="human")
gene.table <- as.data.frame(gene.IDs)
gene.list.entrez <- unique(as.character(gene.table$entrezgene[-which(is.na(gene.table$entrezgene))]))

#entrez
DEGs$entrez <- gene.table$entrezgene[match(rownames(DEGs),gene.table$query)]
DEGs<-DEGs[-which(is.na(DEGs$entrez)),]

write.csv (DEGs,file="./3 ANALISYS/DIFFERENTIAL EXPRESSION/15Oct/MAPKs/DEGsMAPKsLUMINAL.entrez.csv")


## Reactome pathways
yy = enrichPathway(gene.list.entrez, pvalueCutoff=0.1)
head(summary(yy))


#ID                         Description GeneRatio BgRatio       pvalue   p.adjust     qvalue
#380108 380108 Chemokine receptors bind chemokines   46/3815 56/6750 4.388787e-05 0.05586926 0.05446715
#geneID
#380108 6387/3579/7852/3577/3627/6373/4283/5196/5473/6363/6366/10850/2919/2920/2921/6374/6372/6364/10563/6361/6367/6376/2826/6356/6347/6354/1236/6348/56477/
#       58191/1233/1230/729230/1232/1234/9034/10803/10663/2829/1524/6352/6349/414062/6351/6360/1237

## Gene ontology
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

results <- getBM(attributes = c("hgnc_symbol", "go_id","name_1006"),
                 filters = "hgnc_symbol",
                 values = gene.list.name, mart = mart)

GO.data<-results[-which(results$go_id==""),]

#overall

GO.counts <- as.data.frame(table (GO.data$go_id))
GO.ID_Name <- unique(GO.data[,-1])
GO.counts$name_1006 <- GO.ID_Name$name_1006[match(GO.counts$Var1,GO.ID_Name$go_id)] 
GO.counts <- GO.counts[order(-GO.counts$Freq),]

write.csv (GO.counts,file="./3 ANALISYS/GISTIC/GISTIC.BRCA.BSF2/GO.counts.WHX.csv")

#AMP4DEL1 & DEL4AMP1 (using OR)

GO.data.AMP4DEL1 <- GO.data[GO.data$hgnc_symbol %in% AMP4DEL1.genes$genes,]
GO.counts.AMP4DEL1 <- as.data.frame(table (GO.data.AMP4DEL1$go_id))
GO.counts.AMP4DEL1$name_1006 <- GO.ID_Name$name_1006[match(GO.counts.AMP4DEL1$Var1,GO.ID_Name$go_id)] 
GO.counts.AMP4DEL1 <- GO.counts.AMP4DEL1[order(-GO.counts.AMP4DEL1$Freq),]

GO.data.DEL4AMP1 <- GO.data[GO.data$hgnc_symbol %in% DEL4AMP1.genes$genes,]
GO.counts.DEL4AMP1 <- as.data.frame(table (GO.data.DEL4AMP1$go_id))
GO.counts.DEL4AMP1$name_1006 <- GO.ID_Name$name_1006[match(GO.counts.DEL4AMP1$Var1,GO.ID_Name$go_id)] 
GO.counts.DEL4AMP1 <- GO.counts.DEL4AMP1[order(-GO.counts.DEL4AMP1$Freq),]

write.csv (GO.counts.AMP4DEL1,file="./3 ANALISYS/GISTIC/GISTIC.BRCA.BSF2/GO.counts.AMP4DEL1.WHX.csv")
write.csv (GO.counts.DEL4AMP1,file="./3 ANALISYS/GISTIC/GISTIC.BRCA.BSF2/GO.counts.DEL4AMP1.WHX.csv")

#AMP4DEL1 & DEL4AMP1 (using AND)



GO.data.AMP4DEL1 <- GO.data[GO.data$hgnc_symbol %in% AMP4DEL1.genes.2way$genes,]
GO.counts.AMP4DEL1 <- as.data.frame(table (GO.data.AMP4DEL1$go_id))
GO.counts.AMP4DEL1$name_1006 <- GO.ID_Name$name_1006[match(GO.counts.AMP4DEL1$Var1,GO.ID_Name$go_id)] 
GO.counts.AMP4DEL1 <- GO.counts.AMP4DEL1[order(-GO.counts.AMP4DEL1$Freq),]

GO.data.DEL4AMP1 <- GO.data[GO.data$hgnc_symbol %in% DEL4AMP1.genes.2way$genes,]
GO.counts.DEL4AMP1 <- as.data.frame(table (GO.data.DEL4AMP1$go_id))
GO.counts.DEL4AMP1$name_1006 <- GO.ID_Name$name_1006[match(GO.counts.DEL4AMP1$Var1,GO.ID_Name$go_id)] 
GO.counts.DEL4AMP1 <- GO.counts.DEL4AMP1[order(-GO.counts.DEL4AMP1$Freq),]

write.csv (GO.counts.AMP4DEL1,file="./3 ANALISYS/GISTIC/GISTIC.BRCA.BSF2/GO.counts.AMP4DEL1.2way.WHX.csv")
write.csv (GO.counts.DEL4AMP1,file="./3 ANALISYS/GISTIC/GISTIC.BRCA.BSF2/GO.counts.DEL4AMP1.2way.WHX.csv")


