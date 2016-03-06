#########################
## Script to perform mutation spectrum stacked barplot
## Input: Mutation .maf file, and the cluster assignment file (sample name, cluster assignment)
## Modify: Cancer Type (cancer)
##         Paths to mutation file, cluster assignment file, and output filename
## 
######


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
Cancerset  = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF ,Dont use -GA or -hiseq
Geneset    = "DBGS3.FLTR"  # SET GENESET HERE !!!!!!!!!!!!!!
IMS.filter = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" "Luminal A" "Luminal B"
stats      = ""            # Alterantives : "stats" ""
GOF        = "FALSE"

## Load Data
# Load MAF file
load ("./2 DATA/TCGA Mutations/BRCA/Somatic_Mutations/BRCA.TCGA.combined.Mutation.Data.maf.Rdata")
mutation.type.data <- maf.merged.table[,c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","Start_Position","End_Position","Matched_Norm_Sample_Barcode","Line_Number")]
mutation.type.data$Patient_ID <- substr(mutation.type.data$Tumor_Sample_Barcode,1,12)
dim (mutation.type.data) # 90490
mutation.type.data <- unique(mutation.type.data)
dim (mutation.type.data) # 90490
mutation.type.data <- mutation.type.data[,c("Patient_ID","Hugo_Symbol","Variant_Classification","Start_Position","Tumor_Sample_Barcode")]
dim (mutation.type.data) # 87352 
mutation.type.data$Tumor_Sample_Barcode <- NULL
mutation.type.data <- unique(mutation.type.data)
dim (mutation.type.data) # 86921 (start OR end position),87352 (tumor bacode), 89747 (Matched_Norm_Sample_Barcode), 90490( Line_Number)

# Frame_Shift_Ins + Frame_Shift_Del = Frame_Shift
levels (mutation.type.data$Variant_Classification) <- c(levels (mutation.type.data$Variant_Classification),"Frame_Shift")
mutation.type.data[mutation.type.data$Variant_Classification %in% c("Frame_Shift_Ins","Frame_Shift_Del"),"Variant_Classification"] <- "Frame_Shift"
mutation.type.data$Variant_Classification <- gsub("_"," ",mutation.type.data$Variant_Classification)

#cluster assignment
Consensus.class <- read.csv(paste0("./3 ANALISYS/CLUSTERING/RNAseq/",Cancerset,"/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000/",Cancerset,".TCGA.EDASeq.k7.",Geneset,".reps5000.k=4.consensusClass.ICR.csv"),header=TRUE) # select source data
Consensus.class <- Consensus.class[,-1]
colnames (Consensus.class) <- c("Patient_ID","Cluster")
rownames(Consensus.class) <- Consensus.class[,1]

#clinical data
ClinicalData.subset <- read.csv (paste0("./3 ANALISYS/CLINICAL DATA/TCGA.",Cancerset,".RNASeq_subset_clinicaldata.csv"))                       # Clinical data including IMS
rownames(ClinicalData.subset) <- ClinicalData.subset$X 
ClinicalData.subset$X <-NULL


mutation.type.count <- count (mutation.type.data,vars="Patient_ID")
mutation.type.count$Variant_Classification <- "All_types"
mutation.type.count <- mutation.type.count[order(-mutation.type.count$freq),c("Patient_ID","Variant_Classification","freq")]
mutation.type.count.All <- mutation.type.count

mutation.type.count <- count (mutation.type.data,vars=c("Patient_ID","Variant_Classification"))
mutation.type.count$log_freq <- log(mutation.type.count$freq)
#mutation.type.count <- rbind(mutation.type.count,mutation.type.count.type)
sumOFlogs<-aggregate(mutation.type.count$log_freq, by= list(mutation.type.count$Patient_ID), FUN=sum)
colnames(sumOFlogs)<- c("Patient_ID","sum_logs")
mutation.type.count$sum_logs <- sumOFlogs$sum_logs[match(mutation.type.count$Patient_ID,sumOFlogs$Patient_ID)]
mutation.type.count <- mutation.type.count[order(-mutation.type.count$sum_logs,mutation.type.count$Variant_Classification),c("Patient_ID","Variant_Classification","freq","sum_logs")]


## Prepare data for blotting
# Add Class to mutation data
mutation.type.count$Cluster <- Consensus.class$Cluster [match(mutation.type.count$Patient_ID,Consensus.class$Patient_ID)]
mutation.type.count <- mutation.type.count [-which(is.na(mutation.type.count$Cluster)),]
#add subtype to mutation data
mutation.type.count$Subtype <- ClinicalData.subset$TCGA.PAM50.RMethod.RNASeq[match(mutation.type.count$Patient_ID,rownames(ClinicalData.subset))]
mutation.type.count <- mutation.type.count [-which(is.na(mutation.type.count$Subtype)),]

#drop the All_types mutations
mutation.type.count <- mutation.type.count[mutation.type.count$Variant_Classification!="All_types",]

#put silent last
mutation.type.count$Variant_Classification <- factor(mutation.type.count$Variant_Classification)
mutation.type.count$Variant_Classification <- factor(mutation.type.count$Variant_Classification,levels=levels(mutation.type.count$Variant_Classification)[c(1:7,9,8)])
                                                       #c("Frame Shift","In Frame Del","In Frame Ins","Missense Mutation","Nonsense Mutation","Nonstop Mutation","RNA","Splice Site","Silent"))
# ICR 1, 4 only
mutation.type.count.blot <- mutation.type.count[mutation.type.count$Cluster== "ICR3",]



mut_type    = c("Frame Shift","In Frame Del" ,"In Frame Ins" ,"Missense Mutation","Nonsense Mutation","Nonstop Mutation","RNA"    ,"Splice Site","Silent" )
MUT_colors  = c("#00c800"    ,"#aa14f0"      ,"#c45cf5"      ,"#dd3768"          ,"#ec9b2b"          ,"#f0b15a"         ,"#0c9a92","#12e1d5"    ,"grey") #green , purple , auquamarine , orange , redish

# plot
dir.create("./4 FIGURES/mutationspectrum/",showWarnings = FALSE)
tittle = paste0(GOF,".mutation.type.spectrum.from.",IMS.filter,".LOG.bysample.ICR3.")
png (filename = paste0("./4 FIGURES/mutationspectrum/",tittle,".stacked.png") , height = 600, width= 2000)
gg = ggplot(mutation.type.count.blot, aes(x = reorder(Patient_ID,-sum_logs) ,y=log(freq) ,fill = Variant_Classification, order=Variant_Classification ))  +
  geom_bar(stat = "identity", width = 1) + scale_y_continuous(limits = c(0,40)) + # + scale_y_log10(limits = c(1,1e10))
  scale_fill_manual(values = MUT_colors)
gg = gg + ggtitle(tittle) +
  theme(plot.title = element_text(size = 15, lineheight=5, face="bold")) +
  theme_bw() +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme (axis.text.x=element_blank())
#gg = gg + geom_text(data = mutated.allele.blot.bysamp,
#                    aes(y = max(frequency)/10, label = paste0("Tot.Count = ",sample.Count)),
#                   size = 4, vjust = 1.2)
#gg= gg + facet_grid(.~Cluster) +
print(gg)
dev.off()

subtype = c("Basal-like"    ,"HER2-enriched" ,"Luminal A"        ,"Luminal B"          ,"Luminal")
IMS_colors  = c("#da70d6"   ,"#daa520"      ,"#eaff00"          ,"#00c0ff"            ,"#009999" )

subtype.df <- data.frame(Patient_ID=mutation.type.count.blot$Patient_ID,Subtype=mutation.type.count.blot$Subtype,Value=1)
gg= ggplot(subtype.df,aes(x=Patient_ID,y=Value,fill=Subtype))+
  geom_bar(stat = "identity", width = 1)+
  scale_fill_manual(values = IMS_colors)
print(gg)

