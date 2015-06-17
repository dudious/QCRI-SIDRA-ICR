# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("devtools")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
required.packages.BioC <- c("RTCGAToolbox")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
source("http://bioconductor.org/biocLite.R")
biocLite() 
if(length(missing.packages)) biocLite(missing.packages)
library("devtools") # You should install devtools package before loading it.
install_github("mksamur/RTCGAToolbox") # This function will download, build and install RTCGAToolbox.
library("RTCGAToolbox") # Package should be successfully installed!

getFirehoseDatasets()
getFirehoseRunningDates(last=3)
getFirehoseAnalyzeDates(last=3)

brcaData = getFirehoseData (dataset="BRCA",
                            runDate="20150402",
                            gistic2_Date="20141017",
                            Clinic=TRUE,
                            RNAseq_Gene=TRUE,
                            mRNA_Array=TRUE,
                            Mutation=TRUE)
