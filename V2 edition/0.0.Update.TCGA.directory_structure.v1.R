# Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
## download TCGA assembler scripts http://www.compgenome.org/TCGA-Assembler/
source("~/Dropbox/R-projects/QCRI-SIDRA-ICR/R tools/TCGA-Assembler/Module_A.r")
# Scan TCGA database
setwd("./2 DATA/")
TraverseAllDirectories(entryPoint = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/",
                       fileLabel = "DirectoryTraverseResult");
