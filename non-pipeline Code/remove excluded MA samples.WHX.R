load ("./DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.Rdata")
Clinical.data <- read.csv ("./DATA/Agilent_subset_clinicaldata.csv")
rownames(Clinical.data) <- Clinical.data[,1]
Clinical.data[,1] <-NULL
Clinical.data.subset <- subset (Clinical.data,Clinical.data$exclude == "No")
agilentData.subset <- subset (agilentData,is.element (row.names(agilentData),rownames(Clinical.data.subset)))

