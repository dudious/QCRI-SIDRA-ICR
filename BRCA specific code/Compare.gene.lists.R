PAM50.genes = c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
                "NUF2","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4",
                "FOXA1","FOXC1","GPR160","GRB7","KIF2C","NDC80","KRT14","KRT17","KRT5","MAPT",
                "MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","ORC6L",
                "PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")

MAPKMUTSIG.UP <- c("TAOK2","TP53","MAPK3","MAP3K1","MAPT","HSPA1A","FLNB","TAOK3","CRK","RPS6KA2",
                    "MAP2K4","DUSP5","CACNA1D","MAPK8","RASGRP1","CACNA1G")
MAPKMUTSIG.DOWN <- c("CACNG6","CACNA1B","CACNA2D3","FASLG","RASGRF1","JUN","JUND","DUSP16","PPM1B",
                      "SOS1","FGF12","RASGRP2","PRKCB","MAP4K1","PTPN7","GADD45G","DDIT3","DUSP8",
                      "DUSP10","FGFR4","FGF14","FGF13","MAP2K6","DUSP2")

MAPKMUTSIG.UP[which(MAPKMUTSIG.UP %in% PAM50.genes)]
MAPKMUTSIG.DOWN[which(MAPKMUTSIG.DOWN %in% PAM50.genes)]

Master.gene.list.BRCA <- read.csv("./2 DATA/Master_gene_list_breast_cancer.txt",header = FALSE)[,-1]
colnames(Master.gene.list.BRCA) <- c("Source","Symbol","DAVID_ID","Gene_name","Entrez_ID","Unigene_ID","Ensemble_ID","Refseq_ID")

UP.overlap <- MAPKMUTSIG.UP[which(MAPKMUTSIG.UP %in% Master.gene.list.BRCA$Symbol)]
DOWN.overlap <- MAPKMUTSIG.DOWN[which(MAPKMUTSIG.DOWN %in% Master.gene.list.BRCA$Symbol)]

Master.gene.list.BRCA[Master.gene.list.BRCA$Symbol %in% UP.overlap,]
