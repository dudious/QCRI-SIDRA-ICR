load ("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/2 DATA/SUBSETS/BIOLINKS/OV/TCGA.OV.RNASeq.Selected.subset.DBGS3.RData")
load ("f:/DropBox Wouter/Dropbox (TBI-Lab)/BREAST_QATAR/2 DATA/TCGA RNAseq/RNASeq_OV_EDASeq/OV.RNASeq.TCGA.ASSEMBLER.NORMALIZED.LOG2.RData")
ISvsGENE <- data.frame(IS=rowMeans(RNASeq.subset),FLG=RNASeq.NORM_Log2["FLG",rownames(RNASeq.subset)])
corr <- cor.test(x = ISvsGENE$IS, y=ISvsGENE$FLG)
plot (x = ISvsGENE$IS, y=ISvsGENE$FLG,main = paste0("FLGvsIS in OV R:",signif(corr$estimate,2)," p:",signif(corr$p.value,2)))
abline(lm(ISvsGENE$FLG~ISvsGENE$IS), col="red") # regression line (y~x) 
#lines(lowess(ISvsGENE$IS,ISvsGENE$FLG), col="blue") # lowess line (x,y)


