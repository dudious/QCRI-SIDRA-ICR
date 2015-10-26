#########################
## Script to perform mutation kind stats
## Input: Mutation .maf file, and the cluster assignment file (sample name, cluster assignment) AND ONCOTATOR
## Modify: Cancer Type (cancer)
##         Paths to mutation file, cluster assignment file, and output filename
## 
######


## Setup environment
rm(list=ls())
setwd("~/Dropbox/BREAST_QATAR/")
#Dependencies
required.packages <- c("ggplot2", "plyr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
library("ggplot2")
library("plyr")

## Parameters
Cancerset  = "BRCA.BSF2"   # FOR BRCA use BRCA.PCF or BRCA.BSF ,Dont use -GA or -hiseq
Geneset    = "DBGS3.FLTR"  # SET GENESET HERE !!!!!!!!!!!!!!
IMS.filter = "All"         # Alterantives "All" , "Luminal" , "Basal", "Her2" "Luminal A" "Luminal B"
stats      = "stats"       # Alterantives : "stats" ""
GOF        = "FALSE"

## Load Data
# Load MAF file
load ("./2 DATA/TCGA Mutations/BRCA/Somatic_Mutations/BRCA.TCGA.combined.Mutation.Data.maf.Rdata")
load ("./2 DATA/TCGA Mutations/BRCA/Somatic_Mutations/OncotatorResults.rdata")

dim (maf.merged.table)
dim (oncotatorR)
head (oncotatorR)
colnames(oncotatorR)

mutated.kind.count.byGene <- count(oncotatorR,vars=c("Hugo_Symbol","HGVS_protein_change"))
