# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
source("./1 CODE/R tools/TCGA-Assembler/Module_A.r")
# Scan TCGA database
setwd("~/Dropbox/BREAST_QATAR/2 DATA/")
TraverseAllDirectories(entryPoint = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/",
                       fileLabel = "DirectoryTraverseResult");
