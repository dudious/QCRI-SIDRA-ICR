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
Cancerset <- "BRCA.BSF2"        # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset = "DBGS3.FLTR"         # SET GENESET HERE !!!!!!!!!!!!!!
K = 4                          # SET K here
Filter = 1                     # at least one clutser has to have x% mutation frequency
mutation.type = "NonSilent"    # Alterantives "Any" , "Missense" , "NonSilent"
clusters = paste0(rep("ICR",4), 1:4)


## Read the mutation frequency file  (Mutation.Frequency.Gene)
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/Mutation.Data.TCGA.",Cancerset,".",Geneset,".Frequencies.IMS.RDATA"))
#numMuts.Any      = data.frame(count=Mutation.Frequency.Patient$Freq.Any,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Any")
#numMuts.Missense = data.frame(count=Mutation.Frequency.Patient$Freq.Missense,cluster=Mutation.Frequency.Patient$Cluster,mut.type = "Missense")

subtypes = levels(Mutation.Frequency.Patient$Subtype)

## Pick genes based on cutoff (Freq.Any.pct OR Freq.Missense.Any.pct) (present in at least "cutoff" samples for each cluster)
gene.list = as.character(unique(Mutation.Frequency.Gene$Hugo_Symbol))
gene.list.selected = NULL
if (mutation.type=="Any") {
  gene.list.selected = unique(Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Freq.Any.pct>Filter),"Hugo_Symbol"]) #filter one clutser has to have x% mutation
} 
if (mutation.type=="Missense") {
  gene.list.selected = unique(Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Freq.Missense.Any.pct>Filter),"Hugo_Symbol"]) #filter one clutser has to have x% mutation
} 
if (mutation.type=="NonSilent") {
  gene.list.selected = unique(Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Freq.NonSilent.Any.pct>Filter),"Hugo_Symbol"]) #filter one clutser has to have x% mutation
} 
variation.table = NULL

## for each gene, pick the 4 subtypes, corresponding Freq.Any.pct OR Freq.Missense.Any.pct
for (gene in gene.list.selected){
  gene.data = Mutation.Frequency.Gene[which(Mutation.Frequency.Gene$Hugo_Symbol==gene),]   # select a gene
  gene.data[is.na(gene.data)] = 0
  
  # Add missing subtypes data
  gene.subtypes = gene.data$Subtype
  if(length(gene.subtypes)<5){
   missing.subtypes = subtypes[which(!(subtypes %in% gene.subtypes))]
   gene.missingdata = data.frame(matrix(ncol=ncol(gene.data), nrow=length(missing.subtypes)))
   colnames(gene.missingdata) = colnames(gene.data)
   gene.missingdata$Hugo_Symbol = gene
   gene.missingdata$Subtype = missing.subtypes
   gene.missingdata[is.na(gene.missingdata)] = 0
   gene.data = rbind(gene.data, gene.missingdata)
  }
  
  if (mutation.type=="Any") {
    gene.data = gene.data[order(gene.data$Freq.Any.pct),]
    gene.data.pct = gene.data[which(gene.data$Hugo_Symbol==gene),"Freq.Any.pct"]           # add the percentages
    gene.patients = gene.data$Freq.Any
  }
  if (mutation.type=="Missense") {
    gene.data = gene.data[order(gene.data$Freq.Missense.Any.pct),]
    gene.data.pct = gene.data[which(gene.data$Hugo_Symbol==gene),"Freq.Missense.Any.pct"]  # add the percentages
    gene.patients = gene.data$Freq.Missense.Any
  }
  if (mutation.type=="NonSilent") {
    gene.data = gene.data[order(gene.data$Freq.NonSilent.Any.pct),]
    gene.data.pct = gene.data[which(gene.data$Hugo_Symbol==gene),"Freq.NonSilent.Any.pct"] # add the percentages
    gene.patients = gene.data$Freq.NonSilent.Any
  }
  gene.subtype.order = gene.data$Subtype
  ## FIsher exact test ICR1 vs ICR4
  gene.patients[is.na(gene.patients)] = 0
  gene.patients.wt = gene.data$N-gene.patients
  test.matrix = cbind(gene.patients[c(1,5)],gene.patients.wt[c(1,5)])
  res = fisher.test(test.matrix)
  #print(gene)
  #print( res$p.value)
  
  variation = max(gene.data.pct) - min(gene.data.pct)                                      # calculate max variation
  trend = sign(diff(gene.data.pct))                                                        # add direction of change
  trend.test = (trend[which(trend!=0)])                                                    # exclude 0's from trend
  flag = all(trend.test==trend.test[1])                                                    # test for trend 
  if(all(trend.test==0)){flag=FALSE}                                                       # 0,0,0 = no trend
  gene.data[is.na(gene.data)] <- 0
  ## Add Chi-square
  all.patients = gene.data$N
  trend.results = prop.trend.test(gene.patients, all.patients)
  trend.pval = trend.results[[3]]
  db.test=FALSE
  if (gene.data.pct[1]<=1) {
    if ((gene.data.pct[5] - gene.data.pct[1])>=4){
     db.test=TRUE
    }
  }
  if (gene.data.pct[4]<=1) {
    if ((gene.data.pct[1] - gene.data.pct[5])>=4){
      db.test=TRUE
    }
  }  
  gene.data.row = data.frame(gene, paste0(gene.data.pct, collapse=" : "),
                             paste0(gene.subtype.order, collapse=" : "),
                             variation,
                             #paste0((trend), collapse=" : "),
                             #flag,
                             #trend.pval,
                             res$p.value,
                             db.test)  
  variation.table = rbind(variation.table, gene.data.row)
}
colnames(variation.table) = c("Gene",
                              "Subtype_Percentages",
                              "Subtype_Order",
                              "Max_Variation",
                              #"Direction",
                              #"Trend",
                              #"Trend_pVal_ChiSquared",
                              "Fisher_HighestVsLowest",
                              "db.test_HighestVsLowest") 

# significance filter
#SL1 = 0.0005 #chisquare p for LOW OR trend = TRUE
#SL2 = 4   #maxvar Filter multiplier for LOW
#SH1 = SL1/10 #chisquare p for HIGH OR trend = TRUE
#SH2 = SL2*2 #maxvar Filter multiplier for HIGH
#low.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<SL1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*SL2), ]  
#high.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<SH1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*SH2), ] 
# settings Table (SL1,SL2,SH1,SH2)
# BLCA  (0.01,2,0.005,4)
# COAD  (0.0005,4,0.00005,8)

#automatic significance filter
#ASF1     = 0.001
#ASF2     = 1
#ASF.stop = 30
#L.sig = nrow(variation.table)
#while (L.sig > ASF.stop){
#  auto.significant.variation.table = variation.table[which((variation.table$Trend_pVal_ChiSquared<ASF1 | variation.table$Trend) & variation.table$Max_Variation>=Filter*ASF2), ]
#  ASF2 = ASF2+0.1
#  L.sig = nrow(auto.significant.variation.table)
#}

save(#low.significant.variation.table,
     #high.significant.variation.table,
     #auto.significant.variation.table,
     variation.table,
     file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".",mutation.type,".VariationTables.IMS.RData"))
write.csv (variation.table,file=paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",Geneset,".",mutation.type,".VariationTable.IMS.csv"))
