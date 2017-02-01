rm(list=ls())
setwd(("/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/2 DATA/favors/"))
IS <- read.csv ("./file1.csv")
length(which(duplicated(substr(IS$TCGA_ID,1,12))))
IS.filtered <- IS[-which(duplicated(substr(IS$TCGA_ID,1,12))),]
rownames(IS.filtered) <- substr(IS.filtered$TCGA_ID,1,12)
ESTIMATE <- read.csv("./file2.csv")
ESTIMATE$X <- NULL
length(which(duplicated(substr(ESTIMATE$Sample.ID,1,12))))
ESTIMATE.filtered <- ESTIMATE[-which(duplicated(substr(ESTIMATE$Sample.ID,1,12))),]
rownames(ESTIMATE.filtered) <- substr(ESTIMATE.filtered$Sample.ID,1,12)
IS.merged <- merge(IS.filtered,ESTIMATE.filtered,by="row.names")
row.names(IS.merged) <- IS.merged$Row.names
IS.merged$Row.names <- NULL

IS.merged.estimate <-IS.merged[-which(is.na(IS.merged$ESTIMATE)),]
estimate.cor <- cor(x = IS.merged.estimate$Immune_signature_score,
                    y = IS.merged.estimate$ESTIMATE,
                    method = "pearson")
IS.merged.absolute <-IS.merged[-which(is.na(IS.merged$ABSOLUTE)),]
absolute.cor <- cor(x = IS.merged.absolute$Immune_signature_score,
                    y = IS.merged.absolute$ABSOLUTE,
                    method = "pearson")
IS.merged.lump <-IS.merged[-which(is.na(IS.merged$LUMP)),]
lump.cor <- cor(x = IS.merged.lump$Immune_signature_score,
                    y = IS.merged.lump$LUMP,
                    method = "pearson")
IS.merged.ihc <-IS.merged[-which(is.na(IS.merged$IHC)),]
ihc.cor <- cor(x = IS.merged.ihc$Immune_signature_score,
                y = IS.merged.ihc$IHC,
                method = "pearson")
IS.merged.cpe <-IS.merged[-which(is.na(IS.merged$CPE)),]
cpe.cor <- cor(x = IS.merged.cpe$Immune_signature_score,
               y = IS.merged.cpe$CPE,
               method = "pearson")


plot (x = IS.merged.estimate$Immune_signature_score,
      y = IS.merged.estimate$ESTIMATE,
      main = paste0("spearman : ", round(estimate.cor,3)))

plot (x = IS.merged.absolute$Immune_signature_score,
      y = IS.merged.absolute$ABSOLUTE,
      main = paste0("spearman : ", round(absolute.cor,3)))

plot (x = IS.merged.lump$Immune_signature_score,
      y = IS.merged.lump$LUMP,
      main = paste0("spearman : ", round(lump.cor,3)))

plot (x = IS.merged.ihc$Immune_signature_score,
      y = IS.merged.ihc$IHC,
      main = paste0("spearman : ", round(ihc.cor,3)))

plot (x = IS.merged.cpe$Immune_signature_score,
      y = IS.merged.cpe$CPE,
      main = paste0("spearman : ", round(cpe.cor,3)))


load("/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF2.RNASeq_subset_DBGS3.FLTR.Master.Summary.Rdata")
BRCA.data <- Master.file [,c("Freq.NonSilent","unscaled.IS")]


IS.merged.BRCA <- merge (IS.merged,BRCA.data,by="row.names")
rownames(IS.merged.BRCA) <- IS.merged.BRCA$Row.names
IS.merged.BRCA$Row.names <- NULL

BRCA.IS.cor <- cor(x = IS.merged.BRCA$Immune_signature_score,
                   y = IS.merged.BRCA$unscaled.IS,
                   method = "pearson")

IS.merged.BRCA.SM <- IS.merged.BRCA[-which(is.na(IS.merged.BRCA$Freq.NonSilent)),]
BRCA.SM.cor <- cor(x = IS.merged.BRCA.SM$Immune_signature_score,
                   y = IS.merged.BRCA.SM$Freq.NonSilent,
                   method = "pearson")

plot (x = IS.merged.BRCA$Immune_signature_score,
      y = IS.merged.BRCA$unscaled.IS,
      main = paste0("spearman : ", round(BRCA.IS.cor,3)))

plot (x = IS.merged.BRCA.SM$Immune_signature_score,
      y = log(IS.merged.BRCA.SM$Freq.NonSilent),
      main = paste0("spearman : ", round(BRCA.SM.cor,3)))


