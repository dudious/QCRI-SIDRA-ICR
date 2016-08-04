# Setup environment
  rm(list=ls())
  setwd("~/Dropbox/BREAST_QATAR/")
  
# Dependencies
  required.packages <- c("ggplot2", "plyr")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) install.packages(missing.packages)
  library("ggplot2")
  library("plyr")
  
# Load data files
  pct.mut.strong.neoeps.byGene <- read.csv ("./3 ANALISYS/NEOEPITOPES/Fraction of mutations leading to strong neoepitopes per gene", header = TRUE,sep = "\t")
  strong.neoeps.byPatient <- read.csv ("./3 ANALISYS/NEOEPITOPES/Number of strong neoepitopes per sample per ICR", header = TRUE,sep = "\t")
  mutations.strong.neoeps.byPatient <- read.csv ("./3 ANALISYS/NEOEPITOPES/Number of mutations leading to strong neoepitopes per sample per ICR", header = TRUE,sep = "\t")
  pct.mutations.strong.neoeps.byPatient <- read.csv ("./3 ANALISYS/NEOEPITOPES/Fraction of mutations leading to strong neoepitopes per sample per ICR", header = TRUE,sep = "\t")
  load ("./3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF2.RNASeq_subset_DBGS3.FLTR.Master.Summary.Rdata")
  
# merge byPatient Data
  strong.neoeps.byPatient <- strong.neoeps.byPatient[,c("Sample","ICR.group","Val")]
  colnames(strong.neoeps.byPatient) <- c("PatientID","Cluster","Strong.NE.count")
  strong.neoeps.byPatient$SNE.mutation.count <- mutations.strong.neoeps.byPatient$Val[match(strong.neoeps.byPatient$PatientID,mutations.strong.neoeps.byPatient$Sample)]
  strong.neoeps.byPatient$All.mutation.count <- Master.file$Freq.All[match(strong.neoeps.byPatient$PatientID,rownames(Master.file))]
  strong.neoeps.byPatient$SNE.mutation.pct.Y <- pct.mutations.strong.neoeps.byPatient$Val[match(strong.neoeps.byPatient$PatientID,pct.mutations.strong.neoeps.byPatient$Sample)]
  strong.neoeps.byPatient$SNE.mutation.pct.W <- strong.neoeps.byPatient$SNE.mutation.count/strong.neoeps.byPatient$All.mutation.count*100
  strong.neoeps.byPatient$Subtype <- Master.file$TCGA.PAM50.RMethod.RNASeq[match(strong.neoeps.byPatient$PatientID,rownames(Master.file))]

# Master updata
  #Master.file$Strong.NE.count <- strong.neoeps.byPatient$Strong.NE.count[match(rownames(Master.file),strong.neoeps.byPatient$PatientID)]
  #Master.file$SNE.mutation.pct <- strong.neoeps.byPatient$SNE.mutation.pct.W[match(rownames(Master.file),strong.neoeps.byPatient$PatientID)]
  #save(Master.file, file="./3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF2.RNASeq_subset_DBGS3.FLTR.Master.Summary.Rdata")
  #write.csv(Master.file,file="./3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF2.RNASeq_subset_DBGS3.FLTR.Master.Summary.csv")
  


# Filter By IMS
  IMS.filter = "All" #default = "All"     
  if (IMS.filter == "Basal-like") {strong.neoeps.byPatient<-strong.neoeps.byPatient[strong.neoeps.byPatient$Subtype=="Basal-like",]}
  if (IMS.filter == "HER2-enriched") {strong.neoeps.byPatient<-strong.neoeps.byPatient[strong.neoeps.byPatient$Subtype=="HER2-enriched",]}
  if (IMS.filter == "Luminal A") {strong.neoeps.byPatient<-strong.neoeps.byPatient[strong.neoeps.byPatient$Subtype=="Luminal A",]}
  if (IMS.filter == "Luminal B") {strong.neoeps.byPatient<-strong.neoeps.byPatient[strong.neoeps.byPatient$Subtype=="Luminal B",]}
  if (IMS.filter == "Luminal") {strong.neoeps.byPatient<-strong.neoeps.byPatient[strong.neoeps.byPatient$Subtype=="Luminal A" | strong.neoeps.byPatient$Subtype=="Luminal B",]}
  
# Stats
  kruskal.Strong.NE.count <- kruskal.test(data = strong.neoeps.byPatient , Strong.NE.count~Cluster)$p.value
  kruskal.SNE.mutation.count <- kruskal.test(data = strong.neoeps.byPatient , SNE.mutation.count~Cluster)$p.value
  kruskal.SNE.mutation.pct <- kruskal.test(data = strong.neoeps.byPatient , SNE.mutation.pct.W~Cluster)$p.value
  
  
# p-value DMG vs all
  levels(pct.mut.strong.neoeps.byGene$DM2.status)<-c(levels(pct.mut.strong.neoeps.byGene$DM2.status),"REST")
  pct.mut.strong.neoeps.byGene[which(is.na(pct.mut.strong.neoeps.byGene$DM2.status)),]$DM2.status<- "REST"
  kruskal.DMGvsREST <- kruskal.test(data = pct.mut.strong.neoeps.byGene,DM2.status~Val)$p.value
  
# boxplots by cluster
dir.create("./4 FIGURES/neo-epitopes/",showWarnings = FALSE)
cluster.order = c("ICR4", "ICR3", "ICR2", "ICR1")
colors = c("blue", "green", "orange", "red")
#dev.new()
limit = 40
stats = data.frame(Cluster=cluster.order,p.value = paste0("p=",signif(kruskal.Strong.NE.count,3)),Strong.NE.count=NA)
meds = ddply(strong.neoeps.byPatient, .(Cluster), summarise, med = median(Strong.NE.count,na.rm = TRUE))
means = ddply(strong.neoeps.byPatient, .(Cluster), summarise, mean = mean(Strong.NE.count,na.rm = TRUE))
png(paste0("./4 FIGURES/neo-epitopes/StrongNeoEpitopes.TCGA.BRCA.",IMS.filter,".png"), height = 1000, width= 300)   #set filename
    gg = ggplot(strong.neoeps.byPatient, aes(Cluster,Strong.NE.count, fill=Cluster)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    geom_jitter(position=position_jitter(width=0.3,height=0.2))
  gg = gg + ylab("Number of Strong NeoEpitopes per sample") +
    scale_x_discrete(limits=cluster.order) +
    xlab("Clusters") + theme_bw()+ 
    scale_fill_manual(values = colors)
  gg = gg + ggtitle(paste0("StrongNeoEpitopes.TCGA.BRCA.",IMS.filter)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  gg = gg + geom_text(data = stats,
                      aes(y = limit, x = 1.3 ,label = p.value),
                      size = 7, vjust = 1.2)
  gg = gg + geom_text(data = meds,
                      aes(y = -0.25, label = round(med,2)),
                      size = 4, vjust = 1.2)
  gg = gg + geom_text(data = means,
                      aes(y = -1, label = round(mean,2)),
                      size = 4, vjust = 1.2)
  gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                  legend.position = "none",
                  axis.text.x = element_text(size = 12, vjust=1),
                  axis.title.x = element_text(size = 18, vjust = -1),
                  axis.text.y = element_text(size = 18, vjust=1),
                  axis.title.y = element_text(size = 18, vjust = 1))
  gg = gg + coord_cartesian(ylim=c(-2, limit))
print(gg)
dev.off()


#dev.new()
limit = 25
stats = data.frame(Cluster=cluster.order,p.value = paste0("p=",signif(kruskal.SNE.mutation.count,3)),SNE.mutation.count=NA)
meds = ddply(strong.neoeps.byPatient, .(Cluster), summarise, med = median(SNE.mutation.count,na.rm = TRUE))
means = ddply(strong.neoeps.byPatient, .(Cluster), summarise, mean = mean(SNE.mutation.count,na.rm = TRUE))
png(paste0("./4 FIGURES/neo-epitopes/SNE.mutations.TCGA.BRCA.",IMS.filter,".png"), height = 1000, width= 300)   #set filename
  gg = ggplot(strong.neoeps.byPatient, aes(Cluster,SNE.mutation.count, fill=Cluster)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    geom_jitter(position=position_jitter(width=0.3,height=0.2))
  gg = gg + ylab("Number of Mutation giving rise to SNE`s per sample") +
    scale_x_discrete(limits=cluster.order) +
    xlab("Clusters") + theme_bw()+ 
    scale_fill_manual(values = colors)
  gg = gg + ggtitle(paste0("SNE.mutations.TCGA.BRCA.",IMS.filter)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  gg = gg + geom_text(data = stats,
                      aes(y = limit, x = 1.3 ,label = p.value),
                      size = 7, vjust = 1.2)
  gg = gg + geom_text(data = meds,
                      aes(y = -0.25, label = round(med,2)),
                      size = 4, vjust = 1.2)
  gg = gg + geom_text(data = means,
                      aes(y = -1, label = round(mean,2)),
                      size = 4, vjust = 1.2)
  gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                  legend.position = "none",
                  axis.text.x = element_text(size = 12, vjust=1),
                  axis.title.x = element_text(size = 18, vjust = -1),
                  axis.text.y = element_text(size = 18, vjust=1),
                  axis.title.y = element_text(size = 18, vjust = 1))
  gg = gg + coord_cartesian(ylim=c(-2, limit))
print(gg)
dev.off()

#dev.new()
limit = 25
stats = data.frame(Cluster=cluster.order,p.value = paste0("p=",signif(kruskal.SNE.mutation.pct,3)),SNE.mutation.pct.W=NA)
meds = ddply(strong.neoeps.byPatient, .(Cluster), summarise, med = median(SNE.mutation.pct.W,na.rm = TRUE))
means = ddply(strong.neoeps.byPatient, .(Cluster), summarise, mean = mean(SNE.mutation.pct.W,na.rm = TRUE))
png(paste0("./4 FIGURES/neo-epitopes/Pct.SNE.Mutations.TCGA.BRCA.",IMS.filter,".png"), height = 1000, width= 300)   #set filename
  gg = ggplot(strong.neoeps.byPatient, aes(Cluster,SNE.mutation.pct.W, fill=Cluster)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE) +
    geom_jitter(position=position_jitter(width=0.3,height=0.2))
    #, aes(color="lightgray")) 
  gg = gg + ylab("Percent of Mutations giving rise to SNE's per samplee") +
    scale_x_discrete(limits=cluster.order) +
    xlab("Clusters") + theme_bw()+ 
    scale_fill_manual(values = colors)
  gg = gg + ggtitle(paste0("Pct.SNE.Mutations.TCGA.BRCA.",IMS.filter)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  gg = gg + geom_text(data = stats,
                      aes(y = limit, x = 1.3 ,label = p.value),
                      size = 7, vjust = 1.2)
  gg = gg + geom_text(data = meds,
                      aes(y = -0.25, label = round(med,2)),
                      size = 4, vjust = 1.2)
  gg = gg + geom_text(data = means,
                      aes(y = -1, label = round(mean,2)),
                      size = 4, vjust = 1.2)
  gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                  legend.position = "none",
                  axis.text.x = element_text(size = 12, vjust=1),
                  axis.title.x = element_text(size = 18, vjust = -1),
                  axis.text.y = element_text(size = 18, vjust=1),
                  axis.title.y = element_text(size = 18, vjust = 1))
  gg = gg + coord_cartesian(ylim=c(-2, limit))
print(gg)
dev.off()

# boxplots by IMS

if (IMS.filter == "All") {
  kruskal.Strong.NE.count.byIMS <- kruskal.test(data = strong.neoeps.byPatient , Strong.NE.count~Subtype)$p.value
  kruskal.SNE.mutation.count.byIMS <- kruskal.test(data = strong.neoeps.byPatient , SNE.mutation.count~Subtype)$p.value
  kruskal.SNE.mutation.pct.byIMS <- kruskal.test(data = strong.neoeps.byPatient , SNE.mutation.pct.W~Subtype)$p.value
subtype.order = c("HER2-enriched", "Basal-like" , "Luminal B" ,"Luminal A","Normal-like")
colors = c("#daa520","#da70d6","#eaff00","#00c0ff","#d3d3d3")
    
png(paste0("./4 FIGURES/neo-epitopes/StrongNeoEpitopes.TCGA.BRCA.byIMS.png"), height = 1000, width= 300)   #set filename
limit = 40
stats = data.frame(Subtype=subtype.order,p.value = paste0("p=",signif(kruskal.Strong.NE.count.byIMS,3)),Strong.NE.count=NA)
meds = ddply(strong.neoeps.byPatient, .(Subtype), summarise, med = median(Strong.NE.count,na.rm = TRUE))
means = ddply(strong.neoeps.byPatient, .(Subtype), summarise, mean = mean(Strong.NE.count,na.rm = TRUE))
gg = ggplot(strong.neoeps.byPatient, aes(Subtype,Strong.NE.count, fill=Subtype)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.3,height=0.2))
gg = gg + ylab("Number of Strong NeoEpitopes per sample") +
  scale_x_discrete(limits=subtype.order) +
  xlab("Subtypes") + theme_bw()+ 
  scale_fill_manual(values = colors)
gg = gg + ggtitle(paste0("StrongNeoEpitopes.TCGA.BRCA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg = gg + geom_text(data = stats,
                    aes(y = limit, x = 2 ,label = p.value),
                    size = 7, vjust = 1.2)
gg = gg + geom_text(data = meds,
                    aes(y = -0.25, label = round(med,2)),
                    size = 4, vjust = 1.2)
gg = gg + geom_text(data = means,
                    aes(y = -1, label = round(mean,2)),
                    size = 4, vjust = 1.2)
gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                legend.position = "none",
                axis.text.x = element_text(size = 12, vjust =0, angle = 90),
                axis.title.x = element_text(size = 18, vjust = -1),
                axis.text.y = element_text(size = 18, vjust = 1),
                axis.title.y = element_text(size = 18, vjust = 1))
gg = gg + coord_cartesian(ylim=c(-2, limit))
print(gg)
dev.off()


#dev.new()
png(paste0("./4 FIGURES/neo-epitopes/SNE.mutations.TCGA.BRCA.byIMS.png"), height = 1000, width= 300)   #set filename
limit = 25
stats = data.frame(Subtype=subtype.order,p.value = paste0("p=",signif(kruskal.SNE.mutation.count.byIMS,3)),SNE.mutation.count=NA)
meds = ddply(strong.neoeps.byPatient, .(Subtype), summarise, med = median(SNE.mutation.count,na.rm = TRUE))
means = ddply(strong.neoeps.byPatient, .(Subtype), summarise, mean = mean(SNE.mutation.count,na.rm = TRUE))
gg = ggplot(strong.neoeps.byPatient, aes(Subtype,SNE.mutation.count, fill=Subtype)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.3,height=0.2))
gg = gg + ylab("Number of Mutation giving rise to SNE`s per sample") +
  scale_x_discrete(limits=subtype.order) +
  xlab("Subtypes") + theme_bw()+ 
  scale_fill_manual(values = colors)
gg = gg + ggtitle(paste0("SNE.mutations.TCGA.BRCA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg = gg + geom_text(data = stats,
                    aes(y = limit, x = 2 ,label = p.value),
                    size = 7, vjust = 1.2)
gg = gg + geom_text(data = meds,
                    aes(y = -0.25, label = round(med,2)),
                    size = 4, vjust = 1.2)
gg = gg + geom_text(data = means,
                    aes(y = -1, label = round(mean,2)),
                    size = 4, vjust = 1.2)
gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                legend.position = "none",
                axis.text.x = element_text(size = 12, vjust = 0, angle = 90),
                axis.title.x = element_text(size = 18, vjust = -1),
                axis.text.y = element_text(size = 18, vjust = 1),
                axis.title.y = element_text(size = 18, vjust = 1))
gg = gg + coord_cartesian(ylim=c(-2, limit))
print(gg)
dev.off()

#dev.new()
png(paste0("./4 FIGURES/neo-epitopes/Pct.SNE.Mutations.TCGA.BRCA.byIMS.png"), height = 1000, width= 300)   #set filename
limit = 25
stats = data.frame(Subtype=subtype.order,p.value = paste0("p=",signif(kruskal.SNE.mutation.pct.byIMS,3)),SNE.mutation.pct.W=NA)
meds = ddply(strong.neoeps.byPatient, .(Subtype), summarise, med = median(SNE.mutation.pct.W,na.rm = TRUE))
means = ddply(strong.neoeps.byPatient, .(Subtype), summarise, mean = mean(SNE.mutation.pct.W,na.rm = TRUE))
gg = ggplot(strong.neoeps.byPatient, aes(Subtype,SNE.mutation.pct.W, fill=Subtype)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(notch=TRUE) +
  geom_jitter(position=position_jitter(width=0.3,height=0.2))
#, aes(color="lightgray")) 
gg = gg + ylab("Percent of Mutations giving rise to SNE's per samplee") +
  scale_x_discrete(limits=subtype.order) +
  xlab("Subtypes") + theme_bw()+ 
  scale_fill_manual(values = colors)
gg = gg + ggtitle(paste0("Pct.SNE.Mutations.TCGA.BRCA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg = gg + geom_text(data = stats,
                    aes(y = limit, x = 2 ,label = p.value),
                    size = 7, vjust = 1.2)
gg = gg + geom_text(data = meds,
                    aes(y = -0.25, label = round(med,2)),
                    size = 4, vjust = 1.2)
gg = gg + geom_text(data = means,
                    aes(y = -1, label = round(mean,2)),
                    size = 4, vjust = 1.2)
gg = gg + theme(strip.text.x = element_text(size = 20, colour = "black"),
                legend.position = "none",
                axis.text.x = element_text(size = 12, vjust = 0, angle = 90),
                axis.title.x = element_text(size = 18, vjust = -1),
                axis.text.y = element_text(size = 18, vjust = 1),
                axis.title.y = element_text(size = 18, vjust = 1))
gg = gg + coord_cartesian(ylim=c(-2, limit))
print(gg)
dev.off()
}

