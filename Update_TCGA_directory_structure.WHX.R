# Setup environment
rm(list=ls())
library (ggplot2)
setwd("~/Dropbox/BREAST_QATAR/")
source("./CODE/R tools/TCGA-Assembler/Module_A.r")
# Scan TCGA database
TraverseAllDirectories(entryPoint = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/",
                       fileLabel = "DirectoryTraverseResult");