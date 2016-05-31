VT=prop.trend.test(c(27,76,105,74),c(189,295,292,127)) # non-silent count ANY
gof=prop.trend.test(c(27,76,106,75),c(189,295,292,127))
VT$p.value
gof$p.value


MAPX.nonsilent <- prop.trend.test(c(34,32,26,2),c(189,295,292,127))
MAP3K1.nonsilent <- prop.trend.test(c(31,29,20,2),c(189,295,292,127))
MAP2K4.nonsilent <- prop.trend.test(c(10,9,8,0),c(189,295,292,127))
TP53.nonsilent <- prop.trend.test(c(27,76,105,74),c(189,295,292,127))

MAPX.nonsilent$p.value
MAP3K1.nonsilent$p.value
MAP2K4.nonsilent$p.value
TP53.nonsilent$p.value

MAPX.nonsilent <- fisher.test(matrix(c(34,2,189-34,127-2),ncol=2))
MAP3K1.nonsilent <- fisher.test(matrix(c(31,2,189-31,127-2),ncol=2))
MAP2K4.nonsilent <- fisher.test(matrix(c(10,0,189-10,127-0),ncol=2))
TP53.nonsilent <- fisher.test(matrix(c(27,74,189-27,127-74),ncol=2))

MAPX.nonsilent$p.value
MAP3K1.nonsilent$p.value
MAP2K4.nonsilent$p.value
TP53.nonsilent$p.value

MAPX.nonsilent <- chisq.test(matrix(c(34,2,189-34,127-2),ncol=2))
MAP3K1.nonsilent <- chisq.test(matrix(c(31,2,189-31,127-2),ncol=2))
MAP2K4.nonsilent <- chisq.test(matrix(c(10,0,189-10,127-0),ncol=2))
TP53.nonsilent <- chisq.test(matrix(c(27,74,189-27,127-74),ncol=2))

MAPX.nonsilent$p.value
MAP3K1.nonsilent$p.value
MAP2K4.nonsilent$p.value
TP53.nonsilent$p.value

Gof.mut <- read.csv ("./6 Project details/Paper/tables/BRCA.BSF2.All.DBGS3.FLTR.NonSilent.GOF.Mutcounts.csv",stringsAsFactors = FALSE)
colnames(Gof.mut) <- c("variant","Subtype","Cluster","TP53.muts","TP53.percent","MAP2K4.muts","MAP2K4.percent","MAP3K1.muts","MAP3K1.percent",
                       "MAPX.muts","MAPX.percent","N")
Gof.mut <- Gof.mut[-1,]

#TP53 - trend test
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("TP53.muts","N")])
mode(trend.test.matrix) <- "numeric"
TP53.nonsilent.basal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("TP53.muts","N")])
mode(trend.test.matrix) <- "numeric"
TP53.nonsilent.her2.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("TP53.muts","N")])
mode(trend.test.matrix) <- "numeric"
TP53.nonsilent.lumA.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("TP53.muts","N")])
mode(trend.test.matrix) <- "numeric"
TP53.nonsilent.lumB.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("TP53.muts","N")])
mode(trend.test.matrix) <- "numeric"
TP53.nonsilent.normal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("TP53.muts","N")])
mode(trend.test.matrix) <- "numeric"
TP53.nonsilent.luminal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
TP53.nonsilent.basal.trend$p.value
TP53.nonsilent.her2.trend$p.value
TP53.nonsilent.lumA.trend$p.value
TP53.nonsilent.lumB.trend$p.value
TP53.nonsilent.normal.trend$p.value
TP53.nonsilent.luminal.trend$p.value

#TP53 - chisq test
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("TP53.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
TP53.nonsilent.basal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("TP53.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
TP53.nonsilent.her2.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("TP53.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
TP53.nonsilent.lumA.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("TP53.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
TP53.nonsilent.lumB.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("TP53.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
TP53.nonsilent.normal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("TP53.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
TP53.nonsilent.luminal.chisq <- chisq.test(chisq.test.matrix)
TP53.nonsilent.basal.chisq$p.value
TP53.nonsilent.her2.chisq$p.value
TP53.nonsilent.lumA.chisq$p.value
TP53.nonsilent.lumB.chisq$p.value
TP53.nonsilent.normal.chisq$p.value
TP53.nonsilent.luminal.chisq$p.value

#TP53 - chisq test
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("TP53.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
TP53.nonsilent.basal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("TP53.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
TP53.nonsilent.her2.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("TP53.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
TP53.nonsilent.lumA.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("TP53.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
TP53.nonsilent.lumB.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("TP53.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
TP53.nonsilent.normal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("TP53.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
TP53.nonsilent.luminal.fisher <- fisher.test(fisher.test.matrix)
TP53.nonsilent.basal.fisher$p.value
TP53.nonsilent.her2.fisher$p.value
TP53.nonsilent.lumA.fisher$p.value
TP53.nonsilent.lumB.fisher$p.value
TP53.nonsilent.normal.fisher$p.value
TP53.nonsilent.luminal.fisher$p.value

#MAP3K1 - trend test
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAP3K1.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP3K1.nonsilent.basal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAP3K1.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP3K1.nonsilent.her2.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAP3K1.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP3K1.nonsilent.lumA.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAP3K1.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP3K1.nonsilent.lumB.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAP3K1.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP3K1.nonsilent.normal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAP3K1.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP3K1.nonsilent.luminal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
MAP3K1.nonsilent.basal.trend$p.value
MAP3K1.nonsilent.her2.trend$p.value
MAP3K1.nonsilent.lumA.trend$p.value
MAP3K1.nonsilent.lumB.trend$p.value
MAP3K1.nonsilent.normal.trend$p.value
MAP3K1.nonsilent.luminal.trend$p.value

#MAP3K1 - chisq test
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAP3K1.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP3K1.nonsilent.basal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAP3K1.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP3K1.nonsilent.her2.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAP3K1.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP3K1.nonsilent.lumA.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAP3K1.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP3K1.nonsilent.lumB.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAP3K1.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP3K1.nonsilent.normal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAP3K1.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP3K1.nonsilent.luminal.chisq <- chisq.test(chisq.test.matrix)
MAP3K1.nonsilent.basal.chisq$p.value
MAP3K1.nonsilent.her2.chisq$p.value
MAP3K1.nonsilent.lumA.chisq$p.value
MAP3K1.nonsilent.lumB.chisq$p.value
MAP3K1.nonsilent.normal.chisq$p.value
MAP3K1.nonsilent.luminal.chisq$p.value

#MAP3K1 - chisq test
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAP3K1.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP3K1.nonsilent.basal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAP3K1.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP3K1.nonsilent.her2.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAP3K1.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP3K1.nonsilent.lumA.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAP3K1.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP3K1.nonsilent.lumB.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAP3K1.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP3K1.nonsilent.normal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAP3K1.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP3K1.nonsilent.luminal.fisher <- fisher.test(fisher.test.matrix)
MAP3K1.nonsilent.basal.fisher$p.value
MAP3K1.nonsilent.her2.fisher$p.value
MAP3K1.nonsilent.lumA.fisher$p.value
MAP3K1.nonsilent.lumB.fisher$p.value
MAP3K1.nonsilent.normal.fisher$p.value
MAP3K1.nonsilent.luminal.fisher$p.value


#MAP2K4 - trend test
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAP2K4.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP2K4.nonsilent.basal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAP2K4.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP2K4.nonsilent.her2.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAP2K4.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP2K4.nonsilent.lumA.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAP2K4.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP2K4.nonsilent.lumB.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAP2K4.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP2K4.nonsilent.normal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAP2K4.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAP2K4.nonsilent.luminal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
MAP2K4.nonsilent.basal.trend$p.value
MAP2K4.nonsilent.her2.trend$p.value
MAP2K4.nonsilent.lumA.trend$p.value
MAP2K4.nonsilent.lumB.trend$p.value
MAP2K4.nonsilent.normal.trend$p.value
MAP2K4.nonsilent.luminal.trend$p.value

#MAP2K4 - chisq test
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAP2K4.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP2K4.nonsilent.basal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAP2K4.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP2K4.nonsilent.her2.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAP2K4.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP2K4.nonsilent.lumA.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAP2K4.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP2K4.nonsilent.lumB.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAP2K4.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP2K4.nonsilent.normal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAP2K4.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAP2K4.nonsilent.luminal.chisq <- chisq.test(chisq.test.matrix)
MAP2K4.nonsilent.basal.chisq$p.value
MAP2K4.nonsilent.her2.chisq$p.value
MAP2K4.nonsilent.lumA.chisq$p.value
MAP2K4.nonsilent.lumB.chisq$p.value
MAP2K4.nonsilent.normal.chisq$p.value
MAP2K4.nonsilent.luminal.chisq$p.value

#MAP2K4 - chisq test
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAP2K4.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP2K4.nonsilent.basal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAP2K4.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP2K4.nonsilent.her2.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAP2K4.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP2K4.nonsilent.lumA.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAP2K4.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP2K4.nonsilent.lumB.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAP2K4.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP2K4.nonsilent.normal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAP2K4.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAP2K4.nonsilent.luminal.fisher <- fisher.test(fisher.test.matrix)
MAP2K4.nonsilent.basal.fisher$p.value
MAP2K4.nonsilent.her2.fisher$p.value
MAP2K4.nonsilent.lumA.fisher$p.value
MAP2K4.nonsilent.lumB.fisher$p.value
MAP2K4.nonsilent.normal.fisher$p.value
MAP2K4.nonsilent.luminal.fisher$p.value


#MAPX - trend test
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAPX.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAPX.nonsilent.basal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAPX.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAPX.nonsilent.her2.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAPX.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAPX.nonsilent.lumA.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAPX.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAPX.nonsilent.lumB.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAPX.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAPX.nonsilent.normal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
trend.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAPX.muts","N")])
mode(trend.test.matrix) <- "numeric"
MAPX.nonsilent.luminal.trend <- prop.trend.test(trend.test.matrix[,1],trend.test.matrix[,2])
MAPX.nonsilent.basal.trend$p.value
MAPX.nonsilent.her2.trend$p.value
MAPX.nonsilent.lumA.trend$p.value
MAPX.nonsilent.lumB.trend$p.value
MAPX.nonsilent.normal.trend$p.value
MAPX.nonsilent.luminal.trend$p.value

#MAPX - chisq test
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAPX.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAPX.nonsilent.basal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAPX.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAPX.nonsilent.her2.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAPX.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAPX.nonsilent.lumA.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAPX.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAPX.nonsilent.lumB.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAPX.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAPX.nonsilent.normal.chisq <- chisq.test(chisq.test.matrix)
chisq.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAPX.muts","N")])
mode(chisq.test.matrix) <- "numeric"
chisq.test.matrix[,2]<- chisq.test.matrix[,2] - chisq.test.matrix[,1]
chisq.test.matrix<-chisq.test.matrix[-c(2,3),]
MAPX.nonsilent.luminal.chisq <- chisq.test(chisq.test.matrix)
MAPX.nonsilent.basal.chisq$p.value
MAPX.nonsilent.her2.chisq$p.value
MAPX.nonsilent.lumA.chisq$p.value
MAPX.nonsilent.lumB.chisq$p.value
MAPX.nonsilent.normal.chisq$p.value
MAPX.nonsilent.luminal.chisq$p.value

#MAPX - chisq test
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Basal-like",c("MAPX.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAPX.nonsilent.basal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "HER2-enriched",c("MAPX.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAPX.nonsilent.her2.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal A",c("MAPX.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAPX.nonsilent.lumA.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal B",c("MAPX.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAPX.nonsilent.lumB.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Normal-like",c("MAPX.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAPX.nonsilent.normal.fisher <- fisher.test(fisher.test.matrix)
fisher.test.matrix <- as.matrix(Gof.mut[Gof.mut$Subtype == "Luminal",c("MAPX.muts","N")])
mode(fisher.test.matrix) <- "numeric"
fisher.test.matrix[,2]<- fisher.test.matrix[,2] - fisher.test.matrix[,1]
fisher.test.matrix<-fisher.test.matrix[-c(2,3),]
MAPX.nonsilent.luminal.fisher <- fisher.test(fisher.test.matrix)
MAPX.nonsilent.basal.fisher$p.value
MAPX.nonsilent.her2.fisher$p.value
MAPX.nonsilent.lumA.fisher$p.value
MAPX.nonsilent.lumB.fisher$p.value
MAPX.nonsilent.normal.fisher$p.value
MAPX.nonsilent.luminal.fisher$p.value
