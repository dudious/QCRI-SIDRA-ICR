
# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")

# Dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("gage")
library (gage)
biocLite("gageData")
library(gageData)


## Parameters
Cancerset      = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF
Geneset        = "DBGS3.FLTR"  # SET GENESET HERE 
matrix.type    = "NonSilent"   # Alterantives "Any" , "Missense" , "NonSilent"
plot.type      = "db.test.strict"     # Alterantives "low" , "high" , "373genes"  ,"auto"," selected", "db.test", "db.test.strict"
IMS.filter     = "All"     # Alterantives "All" , "Luminal" , "Basal", "Her2" ,"LumA" ,"LumB"
cluster.select = "1vs4"

# Load Data
load (paste0("./3 ANALISYS/Mutations/",Cancerset,"/",Cancerset,".",IMS.filter,".",Geneset,".Mutation.Matrixes.",matrix.type,".Rdata"))


data(gse16873)
data(go.sets.hs)
data(go.subs.hs)
hn=(1:6)*2-1
dcis=(1:6)*2
gse16873.bp.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$BP],
                      ref = hn, samp = dcis)
