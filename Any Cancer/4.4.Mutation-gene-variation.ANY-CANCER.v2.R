#################################################################
###
### This Script calculates the statistics for
### diferentially mutated genes and trends.
###
#################################################################

# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

## Parameters
Cancerset <- "COAD"  # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset = "DBGS3.FLTR"   # SET GENESET HERE !!!!!!!!!!!!!!
K = 4                    # SET K here
Filter = 3               # at least one clutser has to have x% mutation frequency

## Read the mutation frequency file  (Mutation.Frequency.Gene)
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",Geneset,".Frequencies.RDATA"))
numMuts.Any     = data.frame(count=Mutation.Frequency.Patient$Freq.Any,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Any")

## Pick genes based on cutoff (Freq.Any.pct) (present in at least "cutoff" samples for each cluster)
gene.list = as.character(unique(Mutation.Frequency.Gene$Hugo_Symbol))
gene.list.selected = NULL
gene.list.selected = unique(Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Freq.Any.pct>Filter),"Hugo_Symbol"]) #filter one clutser has to have x% mutation
variation.table = NULL

## for each gene, pick the 4 clusters, corresponding Freq.Any.pct
for (gene in gene.list.selected){
  #print(gene)
  gene.data = Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Hugo_Symbol==gene),]   # select a gene
  gene.data = gene.data[order(gene.data$Cluster),]                                         # sort the cluster
  gene.data.pct = gene.data[which(gene.data$Hugo_Symbol==gene),"Freq.Any.pct"]             # add the percentages
  variation = max(gene.data.pct) - min(gene.data.pct)                                      # calculate max variation
  trend = sign(diff(gene.data.pct))                                                        # add direction of change
  trend.test = (trend[which(trend!=0)])                                                    # exclude 0's from trend
  flag = all(trend.test==trend.test[1])                                                    # test for trend 
  if(all(trend.test==0)){flag=FALSE}                                                       # 0,0,0 = no trend
  
  ## Add Chi-square
  gene.patients = gene.data$Freq.Any
  all.patients = gene.data$N
  trend.results = prop.trend.test(gene.patients, all.patients)
  trend.pval = trend.results[[3]]
  
  gene.data.row = data.frame(gene, paste0(gene.data.pct, collapse=" : "), variation, paste0((trend), collapse=" : "), flag, trend.pval)  
  variation.table = rbind(variation.table, gene.data.row)
}
colnames(variation.table) = c("Gene", "Cluster_Percentages", "Max_Variation",  "Direction", "Trend", "Trend_pVal_ChiSquared")

# significance filter
SL1 = 0.0005 #chisquare p for LOW OR trend = TRUE
SL2 = 4   #maxvar Filter multiplier for LOW
SH1 = SL1/10 #chisquare p for HIGH OR trend = TRUE
SH2 = SL2*2 #maxvar Filter multiplier for HIGH
low.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<SL1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*SL2), ]  
high.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<SH1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*SH2), ] 
# settings Table (SL1,SL2,SH1,SH2)
# BLCA  (0.01,2,0.005,4)
# COAD  (0.0005,4,0.00005,8)

#automatic significance filter
ASF1     = 0.001
ASF2     = 1
ASF.stop = 20
L.sig = nrow(variation.table)
while (L.sig > ASF.stop){
  auto.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<ASF1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*ASF2), ]
  ASF2 = ASF2+0.1
  L.sig = nrow(auto.significant.variation.table)
}

save(low.significant.variation.table, high.significant.variation.table,auto.significant.variation.table,variation.table, file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".VariationTables.RData"))
write.csv (variation.table,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".VariationTable.csv"))
