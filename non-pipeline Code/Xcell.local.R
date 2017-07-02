# Setup environment
rm(list=ls())
#setwd("~/Dropbox/BREAST_QATAR")
#setwd("/mnt3/wouter/BREAST-QATAR/")
setwd("~/Dropbox (TBI-Lab)/External Collaborations/BREAST_QATAR")
## dependencies
required.packages.BioC <- c("GSVA")
missing.packages <- required.packages.BioC[!(required.packages.BioC %in% installed.packages()[,"Package"])]
if(length(missing.packages)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(missing.packages)
}
required.packages <- c("devtools")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)


devtools::install_github('dviraran/xCell')
