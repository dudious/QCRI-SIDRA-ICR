TCGA.PAM50.Genes <- read.csv ("./CODE/PAM50.TCGA.method/bioclassifier_R/pam50_annotation.txt", sep="\t")
load ("./DATA/TCGA BC MA/BRCA.MA.TCGA.ASSEMBLER.CLEANED.rda")
agilent.data.Genes <- colnames(agilentData)
missing.genes <- TCGA.PAM50.Genes[which(!is.element (TCGA.PAM50.Genes$pcrID,agilent.data.Genes)),]
