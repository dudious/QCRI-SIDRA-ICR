#########################
## Script to perform specific gene barplots
## Input: Mutation MAF file, and the cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
##         Number of clusters (num.clusters)
##         Paths to mutation file, cluster assignment file, and output filename
## 
#########################


## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("ggplot2", "plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")

## Parameters
Cancerset <- "BRCA"           # do not use -GA or -hiseq (data is merged)
BRCA.Filter <- "BSF2"          # "PCF" or "BSF" Pancer or Breast specific
Geneset <- "DBGS3.FLTR"       # SET GENESET HERE !!!!!!!!!!!!!!
GOF = "MAPX"

## Load Mutaion data
load (paste0("./2 DATA/TCGA Mutations/",Cancerset,"/Somatic_Mutations/",Cancerset,".TCGA.combined.Mutation.Data.maf.Rdata"))
Mutation.selected.data <- data.frame(Hugo_Symbol = maf.merged.table$Hugo_Symbol, Variant_Classification = maf.merged.table$Variant_Classification, Patient_ID = substr(maf.merged.table$Tumor_Sample_Barcode,1,12)) #add Mutation.data$Variant_Type for del,snp,ins

#clinical data
Cancerset = paste0(Cancerset,".",BRCA.Filter)
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL
#add subtype to Mutation.Frequency.Patient Table
Mutation.selected.data$Subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(Mutation.selected.data$Patient_ID,rownames(ClinicalData.subset))]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Subtype)),]
#Count the number of samples by IMS
Count.IMS.samples<-count(unique(Mutation.selected.data[,c("Patient_ID","Subtype")]),vars=c("Subtype"))

#cluster assignment
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]
# Add Class to mutation data
Mutation.selected.data$Cluster <- Consensus.class$Cluster [match(Mutation.selected.data$Patient_ID,Consensus.class$Patient_ID)]
Mutation.selected.data <- Mutation.selected.data [-which(is.na(Mutation.selected.data$Cluster)),]

#collapse IMS and mutation type
# Frame_Shift_Ins + Frame_Shift_Del = Frame_Shift
levels (Mutation.selected.data$Variant_Classification) <- c(levels (Mutation.selected.data$Variant_Classification),"Frame_Shift")
Mutation.selected.data[Mutation.selected.data$Variant_Classification %in% c("Frame_Shift_Ins","Frame_Shift_Del"),"Variant_Classification"] <- "Frame_Shift"
# Luminal A + Luminal B = Luminal
Mutation.selected.data.luminal <- Mutation.selected.data
levels (Mutation.selected.data.luminal$Subtype) <- c(levels (Mutation.selected.data.luminal$Subtype),"Luminal")
Mutation.selected.data.luminal[Mutation.selected.data.luminal$Subtype %in% c("Luminal A","Luminal B"),"Subtype"] <- "Luminal"
Mutation.selected.data.luminal <- Mutation.selected.data.luminal[Mutation.selected.data.luminal$Subtype=="Luminal",]
# No IMS 
Mutation.selected.data.merged_IMS <- Mutation.selected.data
Mutation.selected.data.merged_IMS$Subtype <- "All Subtypes"
# merge table
Mutation.selected.data <- rbind(Mutation.selected.data,Mutation.selected.data.luminal,Mutation.selected.data.merged_IMS)

#count patients / Cluster-IMS
Cluster_IMS.counts <- count(unique(Mutation.selected.data[,c("Patient_ID","Subtype","Cluster")]),vars = c("Subtype","Cluster"))
Cluster_IMS.counts$Group <- paste0(Cluster_IMS.counts$Cluster,".",Cluster_IMS.counts$Subtype)

# TEST : Number of IMS/cluster from Clindata.subset+Consensus.class vs mutationselected.data
#test.df <- ClinicalData.subset[,"TCGA.PAM50.RMethod.RNASeq",drop=FALSE]
#test.df$Cluster <- Consensus.class$Cluster[match(rownames(test.df),Consensus.class$Patient_ID)]
#test.df <- test.df[complete.cases(test.df),]
#test.count <- count(test.df)
#test.count$Group <- paste0(test.count$Cluster,".",test.count$TCGA.PAM50.RMethod.RNASeq)
#test.count$Matched <- Cluster_IMS.counts$freq[match(test.count$Group,Cluster_IMS.counts$Group)]
#NOT EQUAL ??? NOT ALL samples have mutation data

#Filter by GOF
Mutation.selected.data.gof <- Mutation.selected.data [Mutation.selected.data$Hugo_Symbol == GOF,]
if (GOF == "MAPX") {
  Mutation.selected.data.gof <- Mutation.selected.data [Mutation.selected.data$Hugo_Symbol == "MAP3K1" | Mutation.selected.data$Hugo_Symbol == "MAP2K4",]
  Mutation.selected.data.gof$Hugo_Symbol <- GOF
}
if (GOF == "COMBO") {
  Mutation.selected.data.gof <- Mutation.selected.data [Mutation.selected.data$Hugo_Symbol == "MAP3K1" | 
                                                        Mutation.selected.data$Hugo_Symbol == "MAP2K4" |
                                                        Mutation.selected.data$Hugo_Symbol == "CTCF"   |
                                                        Mutation.selected.data$Hugo_Symbol == "FCGBP" ,]
  Mutation.selected.data.gof$Hugo_Symbol <- GOF
}
#count 1 mutation per Mutation class allowed
#count.Any <- count (Mutation.selected.data.gof.Any[,-c(1,3)])
count.Separate.old <- count (Mutation.selected.data.gof,vars=c("Variant_Classification","Subtype","Cluster"))
count.Separate <- count (unique(Mutation.selected.data.gof),vars=c("Variant_Classification","Subtype","Cluster"))
if (GOF == "MAP2K4") {
  count.NonSilent <- data.frame(Variant_Classification = "NonSilent", 
                                count (unique(Mutation.selected.data.gof)[,c("Patient_ID","Subtype","Cluster")],#[-which(Mutation.selected.data.gof$Variant_Classification=="Silent"),]),
                                       vars=c("Subtype","Cluster")))
  count.NonSilent.byIMS.only <- data.frame(Variant_Classification = "NonSilent", 
                                           count (unique(Mutation.selected.data.gof[,c("Patient_ID","Subtype")]),
                                                  vars=c("Subtype")))
} 
if (GOF != "MAP2K4") {
  count.NonSilent <- data.frame(Variant_Classification = "NonSilent", 
                              count (unique(Mutation.selected.data.gof[-which(Mutation.selected.data.gof$Variant_Classification=="Silent"),c("Patient_ID","Subtype","Cluster")]),
                                     vars=c("Subtype","Cluster")))
  count.NonSilent.byIMS.only <- data.frame(Variant_Classification = "NonSilent", 
                                         count (unique(Mutation.selected.data.gof[-which(Mutation.selected.data.gof$Variant_Classification=="Silent"),c("Patient_ID","Subtype")]),
                                                vars=c("Subtype")))
}


write.csv (count.NonSilent,file=paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"count.NonSilent.By.IMS.and.Cluster",".csv", sep=""))
write.csv (count.NonSilent.byIMS.only,file=paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"count.NonSilent.By.IMS.",".csv", sep=""))

#count.All.Missense <- count.All[count.All$Variant_Classification=="Missense_Mutation",] 
blot.df <- rbind(count.Separate,count.NonSilent)
colnames(blot.df) <- c("Mutation_Type","Molecular_Subtype","Cluster_Assignment","Mutation_Count")
#blot.df <- blot.df[-which (blot.df$Variant_Classification %in% c("In_Frame_Del","Silent")),]
#blot.df <- rbind(blot.df,c("Nonsense_Mutation","HER2-enriched","ICR1",0))
#blot.df <- rbind(blot.df,c("Nonsense_Mutation","Luminal","ICR1",0))
#blot.df <- rbind(blot.df,c("Splice_Site","Luminal","ICR4",0))
#blot.df <- rbind(blot.df,c("Frame_Shift","HER2-enriched","ICR1",0))
#blot.df <- rbind(blot.df,c("Frame_Shift","HER2-enriched","ICR3",0))
#blot.df <- blot.df[-27,]

class (blot.df$Mutation_Count) <-"numeric"
blot.df$Group <- paste0(blot.df$Cluster_Assignment,".",blot.df$Molecular_Subtype)
blot.df$Group_Count <- Cluster_IMS.counts$freq[match(blot.df$Group,Cluster_IMS.counts$Group)]
#blot.df$Group_Count <- test.count$freq[match(blot.df$Group,test.count$Group)] #alternative group count
class (blot.df$Group_Count) <-"numeric"
blot.df$Mutation_Frequency <- round((blot.df$Mutation_Count / blot.df$Group_Count) *100,1)

#dud values for MAP2K4
if (GOF == "MAP2K4") {
blot.df <- rbind(blot.df,c("NonSilent","Basal-like","ICR4",0,"ICR4.Basal-like",127,0))  #add missing dud value
blot.df <- rbind(blot.df,c("NonSilent","All Subtypes","ICR4",0,"ICR4.Basal-like",127,0))  #add missing dud value
blot.df <- rbind(blot.df,c("Silent","Basal-like","ICR4",0,"ICR4.Basal-like",127,0))
blot.df <- rbind(blot.df,c("In_Frame_Del","Basal-like","ICR4",0,"ICR4.Basal-like",127,0))
blot.df <- rbind(blot.df,c("In_Frame_Del","All Subtypes","ICR4",0,"ICR4.Basal-like",127,0))
class (blot.df$Mutation_Frequency) <- "numeric"
class(blot.df$Mutation_Count) <- "numeric"
class(blot.df$Group_Count) <- "numeric"
}
#dud values for CTCF
if (GOF == "CTCF") {
  blot.df <- rbind(blot.df,c("In_Frame_Del","Basal-like","ICR4",0,"ICR4.Basal-like",127,0))  #add missing dud value
  blot.df <- rbind(blot.df,c("In_Frame_Del","All Subtypes","ICR4",0,"ICR4.All Subtypes",127,0))  #add missing dud value
  class (blot.df$Mutation_Frequency) <- "numeric"
  class(blot.df$Mutation_Count) <- "numeric"
  class(blot.df$Group_Count) <- "numeric"
}
#dud values for FCGBP
if (GOF == "FCGBP") {
  blot.df <- rbind(blot.df,c("Frame_Shift","Basal-like","ICR4",0,"ICR4.Basal-like",127,0))  #add missing dud value
  blot.df <- rbind(blot.df,c("Frame_Shift","All Subtypes","ICR4",0,"ICR4.Basal-like",127,0))  #add missing dud value
  blot.df <- rbind(blot.df,c("Splice_Site","HER2-enriched","ICR4",0,"ICR4.HER2-enriched",127,0))  #add missing dud value
  blot.df <- rbind(blot.df,c("Splice_Site","All Subtypes","ICR4",0,"ICR4.HER2-enriched",127,0))  #add missing dud value
  blot.df <- rbind(blot.df,c("NonSilent","HER2-enriched","ICR4",1,"ICR4.HER2-enriched",127,0))  #add missing dud value
  class (blot.df$Mutation_Frequency) <- "numeric"
  class(blot.df$Mutation_Count) <- "numeric"
  class(blot.df$Group_Count) <- "numeric"
}


blot.df$Mutation_Type <- gsub("_"," ",blot.df$Mutation_Type)
dir.create(paste0("./4 FIGURES/Mutation Plots/",GOF,"/"), showWarnings = FALSE)
write.csv (blot.df,file=paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,".By.IMS.Cluster.and.Mutationtype",".csv", sep=""))

#remove silent/nonsilent/normal-like
blot.df.stack <- blot.df[blot.df$Mutation_Type!="NonSilent",]
blot.df.stack <- blot.df.stack[blot.df.stack$Mutation_Type!="Silent",]
blot.df.stack <- blot.df.stack[blot.df.stack$Molecular_Subtype!="Normal-like",]
blot.df.stack$Mutation_Type <- as.factor(blot.df.stack$Mutation_Type)

#reorder subtypes
blot.df.stack$Mutation_Type <- droplevels(blot.df.stack$Mutation_Type)
blot.df.stack$Molecular_Subtype <- droplevels(blot.df.stack$Molecular_Subtype)
blot.df.stack <- blot.df.stack[order(factor(blot.df.stack$Molecular_Subtype,levels = c("All Subtypes","Basal-like","HER2-enriched","Luminal A","Luminal B","Luminal"))),]
blot.df.stack$Molecular_Subtype <- factor(blot.df.stack$Molecular_Subtype,levels = c("All Subtypes","Basal-like","HER2-enriched","Luminal A","Luminal B","Luminal"))
print(levels(blot.df.stack$Mutation_Type))
print(levels(blot.df.stack$Molecular_Subtype))

#rescale mutation percentage to fit silent mutation percentages
blot.df.NS<-blot.df[blot.df$Mutation_Type=="NonSilent",]
blot.df.stack$NS.percent.ICR.IMS <- blot.df.NS$Mutation_Frequency[match(blot.df.stack$Group,blot.df.NS$Group)]
NS.subclas.SUM.ICR.IMS <- aggregate(data = blot.df.stack, Mutation_Frequency~Group, sum)
blot.df.stack$NS.subclas.SUM <- NS.subclas.SUM.ICR.IMS$Mutation_Frequency[match(blot.df.stack$Group,NS.subclas.SUM.ICR.IMS$Group)]
blot.df.stack$scaled.freq <- round(blot.df.stack$Mutation_Frequency*blot.df.stack$NS.percent.ICR.IMS/blot.df.stack$NS.subclas.SUM,1)
blot.df.stack[blot.df.stack$scaled.freq == "NaN","scaled.freq"] <- 0

#stats
#create complete matrix for NonSilent
selected.data.NS <- blot.df[0,]
Empty.matrix <- matrix(rep(NA,(ncol(selected.data.NS)*nlevels(blot.df$Molecular_Subtype)*4)),ncol = ncol(selected.data.NS))
colnames(Empty.matrix) <- colnames(selected.data.NS)
selected.data.NS<-rbind (selected.data.NS,Empty.matrix)
selected.data.NS$Cluster_Assignment<-rep(c(paste0("ICR",c(1:4))),nlevels(blot.df$Molecular_Subtype))
selected.data.NS$Molecular_Subtype<-sort(rep(c(levels(blot.df$Molecular_Subtype)),4))
selected.data.NS$Mutation_Type<-"NonSilent"
selected.data.NS$Group <- paste0(selected.data.NS$Cluster_Assignment,".",selected.data.NS$Molecular_Subtype)
selected.data.NS$Group_Count <- Cluster_IMS.counts$freq[match(selected.data.NS$Group,Cluster_IMS.counts$Group)]
blot.df.NS <- blot.df[blot.df$Mutation_Type=="NonSilent",]
selected.data.NS$Mutation_Count <- blot.df.NS$Mutation_Count[match(selected.data.NS$Group,blot.df.NS$Group)]
selected.data.NS[is.na(selected.data.NS$Mutation_Count),"Mutation_Count"]<-0
selected.data.NS$Mutation_Frequency <- round(selected.data.NS$Mutation_Count/selected.data.NS$Group_Count*100,1)

selected.data <- selected.data.NS[selected.data.NS$Molecular_Subtype=="All Subtypes",]
test.trend.all <- prop.trend.test(selected.data$Mutation_Count,selected.data$Group_Count)
square.matrix <- matrix(c(selected.data$Mutation_Count[c(1,4)],(selected.data$Group_Count[c(1,4)]-selected.data$Mutation_Count[c(1,4)])),nrow=2)
chisq.all <- chisq.test(square.matrix)
fisher.all <- fisher.test(square.matrix)
selected.data <- selected.data.NS[selected.data.NS$Molecular_Subtype=="Basal-like",] 
if (all(selected.data$Mutation_Count == 0)== FALSE) {
  test.trend.BL <- prop.trend.test(selected.data$Mutation_Count,selected.data$Group_Count)
}
if (all(selected.data$Mutation_Count == 0)== TRUE) {
  test.trend.BL <- list (p.value = 1)
}
square.matrix <- matrix(c(selected.data$Mutation_Count[c(1,4)],(selected.data$Group_Count[c(1,4)]-selected.data$Mutation_Count[c(1,4)])),nrow=2)
chisq.BL <- chisq.test(square.matrix)
fisher.BL <- fisher.test(square.matrix)
selected.data <- selected.data.NS[selected.data.NS$Molecular_Subtype=="HER2-enriched",] 
test.trend.HE <- prop.trend.test(selected.data$Mutation_Count,selected.data$Group_Count)
square.matrix <- matrix(c(selected.data$Mutation_Count[c(1,4)],(selected.data$Group_Count[c(1,4)]-selected.data$Mutation_Count[c(1,4)])),nrow=2)
chisq.HE <- chisq.test(square.matrix)
fisher.HE <- fisher.test(square.matrix)
selected.data <- selected.data.NS[selected.data.NS$Molecular_Subtype=="Luminal A",] 
test.trend.LA <- prop.trend.test(selected.data$Mutation_Count,selected.data$Group_Count)
square.matrix <- matrix(c(selected.data$Mutation_Count[c(1,4)],(selected.data$Group_Count[c(1,4)]-selected.data$Mutation_Count[c(1,4)])),nrow=2)
chisq.LA <- chisq.test(square.matrix)
fisher.LA <- fisher.test(square.matrix)
selected.data <- selected.data.NS[selected.data.NS$Molecular_Subtype=="Luminal B",] 
test.trend.LB <- prop.trend.test(selected.data$Mutation_Count,selected.data$Group_Count)
square.matrix <- matrix(c(selected.data$Mutation_Count[c(1,4)],(selected.data$Group_Count[c(1,4)]-selected.data$Mutation_Count[c(1,4)])),nrow=2)
chisq.LB <- chisq.test(square.matrix)
fisher.LB <- fisher.test(square.matrix)
selected.data <- selected.data.NS[selected.data.NS$Molecular_Subtype=="Luminal",] 
test.trend.AB <- prop.trend.test(selected.data$Mutation_Count,selected.data$Group_Count)
square.matrix <- matrix(c(selected.data$Mutation_Count[c(1,4)],(selected.data$Group_Count[c(1,4)]-selected.data$Mutation_Count[c(1,4)])),nrow=2)
chisq.AB <- chisq.test(square.matrix)
fisher.AB <- fisher.test(square.matrix)

subtype = c("All Subtypes","Basal-like"    ,"HER2-enriched" ,"Luminal A"        ,"Luminal B"          ,"Luminal")
IMS_colors  = c("#696969", "#da70d6"   ,"#daa520"      ,"#eaff00"          ,"#00c0ff"            ,"#009999" )
mut_type = c("Frame Shift"  ,"In Frame Del" ,"Missense Mutation","Nonsense Mutation","Splice Site" )
MUT_colors  = c("#00c800"   ,"#aa14f0"      ,"#dd3768"          ,"#ec9b2b"          ,"#12e1d5") #green , purple , auquamarine , orange , redish

stats.p.values <- data.frame(Molecular_Subtype = subtype,
                             p.value.trend = paste0("p = ",c(signif(test.trend.all$p.value,digits=3),
                                                             signif(test.trend.BL$p.value,digits=3),
                                                             signif(test.trend.HE$p.value,digits=3),
                                                             signif(test.trend.LA$p.value,digits=3),
                                                             signif(test.trend.LB$p.value,digits=3),
                                                             signif(test.trend.AB$p.value,digits=3))),
                             p.value.chisq = paste0("p = ",c(signif(chisq.all$p.value,digits=3),
                                                             signif(chisq.BL$p.value,digits=3),
                                                             signif(chisq.HE$p.value,digits=3),
                                                             signif(chisq.LA$p.value,digits=3),
                                                             signif(chisq.LB$p.value,digits=3),
                                                             signif(chisq.AB$p.value,digits=3))),
                             p.value.fisher = paste0("p = ",c(signif(fisher.all$p.value,digits=3),
                                                             signif(fisher.BL$p.value,digits=3),
                                                             signif(fisher.HE$p.value,digits=3),
                                                             signif(fisher.LA$p.value,digits=3),
                                                             signif(fisher.LB$p.value,digits=3),
                                                             signif(fisher.AB$p.value,digits=3))),
                             Mutation_Type="NonSilent")

write.csv (stats.p.values,file=paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,".stats.NonSilent",".csv", sep=""))
stats.p.values$Mutation_Type = NA

y.value = 30 #TP53 110 , MAP3K1 35, MAP2K4 7, MAPx 00, CTCF 6.5,FCGBP 18, COMBO 42
if (GOF == "TP53") {y.value = 110}
if (GOF == "MAP3K1") {y.value = 20}
if (GOF == "MAP2K4") {y.value = 15}

#create grid plot of all mutation types en molecular subtypes (stacked bar chart)

png(paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"By.IMS.stacked.noIMScolor",".png", sep=""), height = 500, width= 2000)
 gg = ggplot(blot.df.stack, aes(x = Cluster_Assignment, y = Mutation_Frequency , fill = Mutation_Type )) + #, colour = Molecular_Subtype
              geom_bar(stat="identity",size=0,width=0.8) + #position="dodge",drop=FALSE,ylim=c(0,100)
              facet_grid(.~Molecular_Subtype, space="free") +
              xlab("ICR Cluster Assignment") + ylab("Mutation Frequency") + theme_bw() +
              scale_fill_manual(values = MUT_colors) +
              scale_colour_manual(values = IMS_colors) +
              theme(strip.text.x = element_text(size = 20),strip.text.y = element_text(size = 20)) +
              theme(text = element_text(size=20),axis.text.x = element_text(size=20,angle=0),axis.text.y = element_text(size=20)) +
              ylim(0,y.value) +  
              ggtitle((paste0(GOF,".",Cancerset,".",Geneset,"By.IMS.stacked"))) +
              theme(plot.title = element_text(vjust = 3))
 gg = gg + geom_text(data = stats.p.values,
                     aes(y = y.value, x = 1.3 ,label = p.value.trend),
                     size = 7, vjust = 1.2)
 print(gg)
dev.off()

#drop all Subtypes
blot.df.stack.IMSonly <- blot.df.stack[blot.df.stack$Molecular_Subtype!="All Subtypes",]
stats.p.values <- stats.p.values[stats.p.values$Molecular_Subtype!="All Subtypes",]

png(paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"By.IMS.stacked.noIMScolor.noAllIMS",".png", sep=""), height = 500, width= 2000)
gg = ggplot(blot.df.stack.IMSonly, aes(x = Cluster_Assignment, y = Mutation_Frequency , fill = Mutation_Type )) + #, colour = Molecular_Subtype
  geom_bar(stat="identity",size=0,width=0.8) + #position="dodge" drop=FALSE,ylim=c(0,100)
  facet_grid(.~Molecular_Subtype, space="free") +
  xlab("ICR Cluster Assignment") + ylab("Mutation Frequency") + theme_bw() +
  scale_fill_manual(values = MUT_colors) +
  scale_colour_manual(values = IMS_colors) +
  theme(strip.text.x = element_text(size = 20),strip.text.y = element_text(size = 20)) +
  theme(text = element_text(size=20),axis.text.x = element_text(size=20,angle=0),axis.text.y = element_text(size=20)) +
  ylim(0,y.value) +  #TP53 115 , 
  ggtitle((paste0(GOF,".",Cancerset,".",Geneset,"By.IMS.stacked"))) +
  theme(plot.title = element_text(vjust = 3))
gg = gg + geom_text(data = stats.p.values,
                    aes(y = y.value, x = 1.3 ,label = p.value.trend),
                    size = 7, vjust = 1.2)
print(gg)
dev.off()

#SCALED version

Labels <-selected.data.NS [selected.data.NS$Molecular_Subtype!="All Subtypes"&selected.data.NS$Molecular_Subtype!="Normal-like",
                           c("Molecular_Subtype","Cluster_Assignment","Mutation_Frequency")]
Labels <- Labels[order(match(Labels$Molecular_Subtype,subtype)),]
Labels$Mutation_Type <- NA
Labels$Molecular_Subtype <- factor(Labels$Molecular_Subtype,levels=subtype)

png(paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,"By.IMS.stacked.scaled.noIMScolor.noAllIMS",".png", sep=""), height = 500, width= 2000)
gg = ggplot(blot.df.stack.IMSonly, aes(x = Cluster_Assignment, y = scaled.freq , fill = Mutation_Type )) + #, colour = Molecular_Subtype
  geom_bar(stat="identity",size=0,width=0.8) + #position="dodge"
  facet_grid(.~Molecular_Subtype, space="free") +
  xlab("ICR Cluster Assignment") + ylab("Mutation Frequency") + theme_bw() +
  scale_fill_manual(values = MUT_colors) +
  scale_colour_manual(values = IMS_colors) +
  theme(strip.text.x = element_text(size = 20),strip.text.y = element_text(size = 20)) +
  theme(text = element_text(size=20),axis.text.x = element_text(size=20,angle=0),axis.text.y = element_text(size=20)) +
  ylim(0,y.value) +  
  ggtitle((paste0(GOF,".",Cancerset,".",Geneset,"Scaled.By.IMS.stacked"))) +
  theme(plot.title = element_text(vjust = 3))
gg = gg + geom_text(data = stats.p.values,
                    aes(y = y.value, x = 1.3 ,label = p.value.trend),
                    size = 7, vjust = 1.2)
gg = gg + geom_text(data = Labels,
                    aes(y = Mutation_Frequency+2, x = Cluster_Assignment ,label = Mutation_Frequency),
                    size = 7, vjust = 1.2)
print(gg)
dev.off()


#simplyfy to all subtypes

blot.df.stack.noIMS <- blot.df.stack[blot.df.stack$Molecular_Subtype=="All Subtypes",]
if (GOF == "TP53") {y.value = 70}
if (GOF == "MAP3K1") {y.value = 15}
if (GOF == "MAP2K4") {y.value = 7}
if (GOF == "MAPX") {y.value = 20}
#simplyfy to all mutation types 

blot.df <- blot.df[blot.df$Mutation_Type=="NonSilent",]
blot.df <- blot.df[blot.df$Molecular_Subtype=="All Subtypes",]
test.trend <- prop.trend.test(blot.df$Mutation_Count,blot.df$Group_Count)

png(paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,".stacked.scaled.png", sep=""), height = 500, width= 700)
gg = ggplot(blot.df.stack.noIMS, aes(x = Cluster_Assignment, y = scaled.freq , fill = Mutation_Type )) + #, colour = Molecular_Subtype
  geom_bar(stat="identity",size=0,width=0.8) + #position="dodge",drop=FALSE,ylim=c(0,100)
  xlab("ICR Cluster Assignment") + ylab("Mutation Frequency (%)") + theme_bw() +
  scale_fill_manual(values = MUT_colors) +
  scale_colour_manual(values = IMS_colors) +
  theme(strip.text.x = element_text(size = 20),strip.text.y = element_text(size = 20)) +
  theme(text = element_text(size=20),axis.text.x = element_text(size=20,angle=0),axis.text.y = element_text(size=20)) +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,y.value) +
  ggtitle((paste0(GOF,".",Cancerset,".",Geneset,"Scaled.By.IMS.stacked"))) + 
  theme(plot.title = element_text(vjust = 3))+
  annotate("text", label = paste0("trend : p = ",signif(test.trend.all$p.value,digits=3)), size = 6, x=1.5 , y = y.value) + #(nlevels(blot.df$Cluster_Assignment)-1)
  annotate("text", label = blot.df$Mutation_Frequency, size = 6, x=blot.df$Cluster_Assignment , y = blot.df$Mutation_Frequency+1)
print(gg)
dev.off()



png(paste0("./4 FIGURES/Mutation Plots/",GOF,"/",GOF,".",Cancerset,".",Geneset,".Overall.png", sep=""))
#cluster.order = rev(clusters)
gg = ggplot(blot.df, aes(x = Cluster_Assignment, y = Mutation_Frequency))  +
  geom_bar(stat = "identity", width = 0.8, position="dodge", fill = c("blue", "green", "orange", "red")) +
  annotate("text", label = paste0("p = ",signif(test.trend.all$p.value,digits=3)), size = 6, x=1 , y = y.value) + #(nlevels(blot.df$Cluster_Assignment)-1)
  geom_text(label=paste0(blot.df$Mutation_Frequency," %"), size = 6,y=y.value)
gg = gg + #scale_x_discrete(limits = cluster.order) +
  xlab("ICR Cluster Assignment") + ylab("Mutation Frequency") + theme_bw() +
  theme(text = element_text(size=18),axis.text.x = element_text(size=18),axis.text.y = element_text(size=18)) +
  ylim(0,y.value)
gg = gg+ theme(legend.position="none", strip.text.x = element_text(size = 18)) +
  #scale_fill_manual(values=c("gold")) +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle((paste0("Mutation Frequency - ",GOF))) #,".",Cancerset,".",Geneset
print(gg)
dev.off()

