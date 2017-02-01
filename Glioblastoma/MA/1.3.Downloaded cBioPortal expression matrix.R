load ("./2 DATA/SUBSETS/ASSEMBLER/GBM/TCGA.GBM.MA.subset.DBGS1.RData")
cPB.table <- read.csv ("./2 DATA/SUBSETS/ASSEMBLER/GBM/GBM.cBP.txt", sep="\t",header = FALSE)
cPB.table <- cPB.table[-c(1,5,12),]
cPB.table <- t(cPB.table[,-1])
colnames (cPB.table) <- cPB.table[1,]
cPB.table <- cPB.table[-1,]
rownames(cPB.table) <- substr(cPB.table[,1],1,12)
cPB.table <- cPB.table[,-1]
cPB.table.MA <- cPB.table[-which(cPB.table[,2]=="NaN"),]
cPB.table.MA <- cPB.table.MA[complete.cases(cPB.table.MA),]
cPB.table.test <- as.data.frame(cPB.table.MA[,c("IFNG","CXCL10")])
MA.subset <- as.data.frame(MA.subset)
cPB.table.test$CXCL10.ASMLBR <- MA.subset$CXCL10[match(rownames(cPB.table.test),rownames(MA.subset))]
cPB.table.test <- cPB.table.test[order(cPB.table.test$CXCL10),]
