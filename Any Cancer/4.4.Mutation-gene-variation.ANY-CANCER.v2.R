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
Cancerset <- "BRCA.BSF"  # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset = "DBGS3.FLTR"   # SET GENESET HERE !!!!!!!!!!!!!!
K = 4                    # SET K here
Filter = 3               # at least one clutser has to have x% mutation frequency

## Read the mutation frequency file  (Mutation.Frequency.Gene)
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",Geneset,".Frequencies.RDATA"))
numMuts.Any     = data.frame(count=Mutation.Frequency.Patient$Freq.Any,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Any")

## Pick genes based on cutoff (Freq.Any.pct) (present in at least "cutoff" samples for each cluster)
gene.list = as.character(unique(Mutation.Frequency.Gene$Hugo_Symbol))
gene.list.selected = NULL
gene.list.selected = unique(Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Freq.Any.pct>Filter),"Hugo_Symbol"]) #filter one clutser has to have 10% mutation
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
low.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<.1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*0.75), ]
high.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<.05 | variation.table$Trend) & variation.table$Max_Variation>=Filter*1.5), ]

save(low.significant.variation.table, high.significant.variation.table,variation.table, file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".VariationTables.RData"))
write.csv (variation.table,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".VariationTable.csv"))
