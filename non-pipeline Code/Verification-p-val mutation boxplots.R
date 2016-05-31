#Verification of p-val statistics

master.file = read.csv("3 ANALISYS/MASTER FILES/TCGA.BRCA.BSF2.RNASeq_subset_DBGS3.FLTR.Master.csv", header =T, stringsAsFactors = F)
master.file.sub = data.frame(Sample = master.file$X, Cluster = master.file$Cluster.DBGS3.FLTR.RNSeq, Subtype = master.file$TCGA.PAM50.RMethod.RNASeq, Freq.NonSilent = master.file$Freq.NonSilent, Freq.All = master.file$Freq.All)

## remove samples with no cluster assignment
master.file.sub = master.file.sub[-which(is.na(master.file.sub$Cluster)),]
#master.file.sub = master.file.sub[-which(is.na(master.file.sub$Subtype)),]

## p-val of all non-silent by cluster
test.NonSilent.All = kruskal.test(master.file.sub$Freq.NonSilent, master.file.sub$Cluster)
p.value.NonSilent.All    = test.NonSilent.All$p.value #summary(test.NonSilent.Any)[[1]][["Pr(>F)"]][[1]]

test.All.All = kruskal.test(master.file.sub$Freq.All, master.file.sub$Cluster)
p.value.All.All    = test.All.All$p.value #summary(test.NonSilent.Any)[[1]][["Pr(>F)"]][[1]]

##luminalA
test.Basal.data = master.file.sub[which(master.file.sub$Subtype=="Basal-like"),]
test.NonSilent.Basal = kruskal.test(test.Basal.data$Freq.NonSilent, test.Basal.data$Cluster)
p.value.NonSilent.Basal    = test.NonSilent.Basal$p.value #summary(test.NonSilent.Any)[[1]][["Pr(>F)"]][[1]]

test.All.Basal = kruskal.test(test.Basal.data$Freq.All, test.Basal.data$Cluster)
p.value.All.Basal    = test.All.Basal$p.value #summary(test.NonSilent.Any)[[1]][["Pr(>F)"]][[1]]
