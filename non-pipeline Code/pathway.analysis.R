
# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

# Dependencies
source("https://bioconductor.org/biocLite.R")
required.packages <- c("ReactomePA","mygene","biomaRt","tmod")
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

## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "db.test.strict"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict"
IMS.filter     = "All"     # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "1vs4"

# Load Data
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))
diff.mut.genes <- read.csv ("./3 ANALISYS/Mutations/BRCA.BSF2/diferentially.mutated.genes.dbtest.loose.csv")
diff.mut.genes.melted <- melt(diff.mut.genes[,c("Luminal.A","Luminal.B","Her2","basal","All")],measure.vars=c("Luminal.A","Luminal.B","Her2","basal","All"))
overlap <- as.data.frame(table (diff.mut.genes.melted$value))
table(overlap$Freq)

gene.list.name  <- unique(diff.mut.genes.melted$value)
gene.list.name  <- gene.list.name[-which(gene.list.name %in% "")]
length (gene.list.name) #5156

gene.IDs <- queryMany(gene.list.name, scopes="symbol", fields=c("entrezgene", "go"), species="human")
gene.table <- as.data.frame(gene.IDs)
gene.list.entrez <- as.character(gene.table$entrezgene[-which(is.na(gene.table$entrezgene))])

## Reactome pathways
yy <- enrichPathway(gene.list.entrez,organism = "human", pvalueCutoff = 0.05,
              pAdjustMethod = "none", qvalueCutoff = 0.2,
              readable = TRUE)
barplot(yy, showCategory = 15)
Enrichment.result <- data.frame(ID= yy@result$ID ,class = yy@result$Description,genes=yy@result$geneID,p.value=yy@result$p.adjust)
write.csv ( Enrichment.result , file = "./3 ANALISYS/Mutations/BRCA.BSF2/Enriched.pathways.",IMS.filter ,".csv")

## Gene ontology
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("hgnc_symbol", "go_id","name_1006"),
                 filters = "hgnc_symbol",
                 values = gene.list.name, mart = mart)
GO.data<-results[-which(results$go_id==""),]
GO.counts <- as.data.frame(table (GO.data$go_id))
GO.filtered <- GO.counts[GO.counts$Freq >= 100 ,] #101 GO's left
GO.data.filtered <- GO.data[GO.data$go_id %in% GO.filtered$Var1, ]
#Go.data.sign <- GO.data[GO.data$go_id %in% GO.filtered$Var1, ]
write.csv ( GO.data.filtered , file = "./3 ANALISYS/Mutations/BRCA.BSF2/GO.Datat.bygene.",IMS.filter ,".csv")

# methode chausebell
load ("./2 DATA/TCGA Mutations/BRCA/Somatic_Mutations/BRCA.TCGA.combined.Mutation.Data.maf.Rdata")

fg = gene.list.name
bg = unique (maf.merged.table$Hugo_Symbol)
tmod.result = tmodHGtest(fg=fg, bg=bg)
tmodCERNOtest(fg)


save(GO.filtered,GO.data.filtered,Enrichment.result,file="./3 ANALISYS/Mutations/BRCA.BSF2/diferentially.mutated.dbtest.loose.GO.Rdata")
