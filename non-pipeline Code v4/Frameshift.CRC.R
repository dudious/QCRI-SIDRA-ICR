# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")                                                                   # Setwd to location were output files have to be saved.

load ("./3_DataProcessing/TCGA_Assembler/COAD/SomaticMutationData/COAD_Pairs_Aggregated_Capture_Processed_mutationLevel.Rdata")

summary(as.factor(Des$Variant_Classification))
#15843 FS_del & 8666 FS_ins

Des.FS <- Des[Des$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins"),]
Data.FS <- Data[Des$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins"),]

Des.FS$count <- rowSums(Data.FS)
Des.FS$pct <- (Des.FS$count / ncol(Data.FS))*100

Des.FS.Filtered <- Des.FS[Des.FS$pct > 3,]


