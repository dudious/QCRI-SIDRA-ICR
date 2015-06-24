
# Setup environment
rm(list=ls())
## dependencies
## install java for xlsx export
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
required.packages <- c("xlsx","plyr","ggplot2","reshape")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library (xlsx) #xlsx needs java installed
library (plyr)
library (reshape)
library (ggplot2)
setwd("~/Dropbox/BREAST_QATAR/")

# Load data files
load ("./2 DATA/SUBSETS/RNASeq.subset.TP53andCo.RData")

TP53.patients <- read.csv ("./3 ANALISYS/Mutations/TP53.patients.csv")                                                      
rownames(TP53.patients) <- TP53.patients$X
TP53.patients$X <- NULL

ClinicalData.subset <- read.csv ("./3 ANALISYS/CLINICAL DATA/TCGA.BRCA.RNASeq_subset_clinicaldata.csv")                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL

RNASeq.subset<-merge (RNASeq.subset,TP53.patients[,c("Mutated","Missense")],by="row.names",all.x=FALSE)
row.names(RNASeq.subset) <- RNASeq.subset$Row.names
RNASeq.subset$Row.names <- NULL
RNASeq.subset<-merge (RNASeq.subset,ClinicalData.subset[,c("TCGA.PAM50.RMethod.MA"),drop=FALSE],by="row.names",all.x=TRUE)
row.names(RNASeq.subset) <- RNASeq.subset$Row.names
RNASeq.subset$Row.names <- NULL

means <- ddply(RNASeq.subset,c("TCGA.PAM50.RMethod.MA","Mutated"), numcolwise(mean))
sds<- ddply(RNASeq.subset,c("TCGA.PAM50.RMethod.MA","Mutated"), numcolwise(sd))

melted <- melt(as.matrix(RNASeq.subset[,1:6]), id.vars = "row.names")
colnames(melted) <- c("Patient_ID","Gene","Expr")
melted.mut <- merge(melted,RNASeq.subset[,c("Mutated","TCGA.PAM50.RMethod.MA"),drop=FALSE],by.y="row.names",by.x="Patient_ID",all.X=TRUE)
melted.mut <- melted.mut[-which(is.na(melted.mut$TCGA.PAM50.RMethod.MA)),]

ggplot(melted.mut, aes(x=TCGA.PAM50.RMethod.MA, y=Expr, fill=Mutated)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean,
               geom="point",
               shape=5,
               size=4) + 
  facet_grid(. ~Gene) +
  scale_fill_manual(values = c("blue", "red"), .3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("WTvsMUT Expression of selected genes by subtype")
