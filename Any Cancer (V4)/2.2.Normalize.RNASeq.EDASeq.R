#################################################################
###
### This script normalizes ANY-CANCER cancer RNASeq Data 
### from the TCGA database
### It will process the data into a R data file. 
### Data is saved in:
### "./3_Dataprocessing/",download.method, "/", Cancer, "/RNASeqData/"
###
#################################################################

## Normalize RNASeq data using TCGA EDASeq

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                  # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                        # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_A.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_B.R"))

required.bioconductor.packages = c("SummarizedExperiment","EDASeq", "preprocessCore")
required.packages = c("base64enc", "HGNChelper","RCurl","httr","stringr","digest","bitops",
                      "rjson","Matrix","latticeExtra","matrixStats")
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                     # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"
assay.platform = "gene_RNAseq"
Log_file = paste0("./1_Log_Files/2.2_RNASeq_Normalization/RNASeq_Normalization_Log_File_",                              # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
TCGASampleTypeFile = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/TCGASampleType.txt")

# Load data
load(paste0(code_path, "R tools/geneInfo.July2017.RData"))
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

#update geneinfo file
#colnames(geneInfo)[which(colnames(geneInfo) == "entrezgene")] <- "EntrezID"
#gene.query <- queryMany(geneInfo$hgnc_symbol, scopes="symbol", fields=c("entrezgene", "go"), species="human")
#gene.table <- as.data.frame(gene.query)
#geneInfo$new_EntrezID = gene.table$entrezgene[match(geneInfo$hgnc_symbol,gene.table$query)]
#save(geneInfo,gene.table,file = paste0(code_path,"R tools/geneInfo.July2017.RData"))

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/2.2_RNASeq_Normalization/"), showWarnings = FALSE)
cat("This is a log file for RNASeq Normalization",                                                                      # Set-up logfile
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),                                                          
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    "",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

start.time.all <- Sys.time()

i=1
# Normalization
for (i in 1:N.sets) {
  start.time.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  cat(paste0(Cancer, ": "), file = Log_file, append = TRUE, sep = "\n")
  if (Cancer %in% Cancer_skip) {next}
  print(paste0 ("Processing ",Cancer,"."))
  
  load(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/RNASeqData/",Cancer, "_", 
              assay.platform, "_", "Processed.Rdata"))
  
  Des = as.data.frame(Des, StringsAsfactors = FALSE)
  Des$hgnc_symbol = geneInfo$hgnc_symbol[match(Des$EntrezID, geneInfo$new_EntrezID)]
  #lost.genes <- as.character(Des[which(is.na(Des$hgnc_symbol)),"GeneSymbol"])
  #lost.genes <- lost.genes[-which(lost.genes %in% Des$hgnc_symbol)]
  #which(duplicated(lost.genes))
  Des$hgnc_symbol[Des$hgnc_symbol=="?"] <- NA
  
  Raw.data = Data[, -(grep(".*\\.1$",colnames(Data)))]                                                                      # Drop columns with colname ending with .1 in data-file (these are scaled values) and save in Raw.data
  
  cat(paste0("Number of samples in Data is ", ncol(Raw.data), ".\n", "number of genes in Data is ", 
             nrow(Raw.data), ".", "\n", "Number of patients is ", 
             length(unique(substring(colnames(Raw.data),1,12))), "."), 
      file=Log_file, sep = "\n", append = TRUE)
  tissue.type = c("TP", "TR", "TAP", "TM", "TAM", "NT") 
  
  sink(Log_file, append=TRUE, split=TRUE)
  
  RNAseqData.to.normalize = ExtractTissueSpecificSamples(inputData = Raw.data,
                                                                     tissueType = tissue.type,
                                                                     singleSampleFlag = FALSE,
                                                                     sampleTypeFile = TCGASampleTypeFile)
  sink()
  
  # Drop genes without info in the RNAseq data. Drop genes without RNAseq data in the info dataframe.
  if(is.null(nrow(RNAseqData.to.normalize))) {next}
  RNAseqData.to.normalize = cbind(Des, RNAseqData.to.normalize, stringsAsFactors = FALSE)
  cat(paste0(length(which(is.na(RNAseqData.to.normalize$hgnc_symbol))), " genes that do not have a hgnc_symbol",
                    " were deleted."), file = Log_file, append = TRUE, sep = "\n")
  RNAseqData.to.normalize = RNAseqData.to.normalize[-which(is.na(RNAseqData.to.normalize$hgnc_symbol)),]                   # Delete genes without hgnc gene symbol 
  
  rownames(RNAseqData.to.normalize) = RNAseqData.to.normalize$hgnc_symbol                                                  # set row.names to hgnc gene symbol
  RNAseqData.to.normalize = RNAseqData.to.normalize[ , -which(names(RNAseqData.to.normalize) %in% c("GeneSymbol",          # delete first 3 columns of RNASeqData
                                                                                                    "EntrezID",
                                                                                                    "hgnc_symbol"))]
  RNAseqData.to.normalize = as.matrix(RNAseqData.to.normalize)

  # drop rows(genes) without approved gene-symbol
  RNAseq.genes = rownames(RNAseqData.to.normalize)
  info.genes = rownames(geneInfo)
  available.genes = unique(RNAseq.genes[which(RNAseq.genes %in% info.genes)])
  geneInfo = geneInfo[which(rownames(geneInfo) %in% available.genes),]                                                      # drop the genes without RNAseq.DATA
  RNAseqData.to.normalize = RNAseqData.to.normalize[which(rownames(RNAseqData.to.normalize) %in% available.genes),]         # drop the genes without info
  mode(RNAseqData.to.normalize) <- "numeric"
  
  geneInfo <- geneInfo[ order(row.names(geneInfo)), ]
  RNAseqData.to.normalize <- floor(RNAseqData.to.normalize[order(row.names(RNAseqData.to.normalize)),])
  
  cat(paste0(length(RNAseqData.to.normalize[,1]), " genes were available for normalization."), 
      file = Log_file, append = TRUE, sep = "\n")                                                                           # get information on number of genes included for normalization to logfile
  
  
  RNASeq.expr.set = newSeqExpressionSet(RNAseqData.to.normalize, featureData = geneInfo)                                    # Create a new SeqExpressionSet object.
  fData(RNASeq.expr.set)[, "gcContent"] = as.numeric(geneInfo[, "gcContent"])                                               # Make sure gcContent is numeric
  RNASeq.expr.set = withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE)                   # Removes lane gene specific effects, for example effects related to gene length or GC content
  RNASeq.expr.set = betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)                               # Removes effect related to in between lane distributional differences, as sequencing depth
  RNASeq.NORM = log(RNAseqData.to.normalize + .1) + offst(RNASeq.expr.set)                                                  # Apply the Edaseq Ofset
  RNASeq.NORM = floor(exp(RNASeq.NORM) - .1)                                                                                # Return non decimal values
  
  #Quantile normalisation RNA
  RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                                                                 # Quantile normalize
  RNASeq.NORM.quantiles <- floor(RNASeq.NORM.quantiles)                                                                     # Return non decimal values
  rownames(RNASeq.NORM.quantiles) <- rownames(RNASeq.NORM)
  colnames(RNASeq.NORM.quantiles) <- colnames(RNASeq.NORM)
  
  Rdata.file = paste0("./3_Dataprocessing/",download.method, "/", Cancer, "/RNASeqData/",Cancer, "_", 
         assay.platform, "_normalized.Rdata")
  save(RNASeq.NORM.quantiles,geneInfo, Des,
       file= Rdata.file)
  
  end.time.cancer = Sys.time()
  time = substring(as.character(capture.output(round(end.time.cancer - start.time.cancer, 2))),20,100)
  msg = paste0("Normalization time for ", Cancer, ": ",time, ".", "\n", "Outputfile is ", Rdata.file, "\n",
               " which contains RNASeq.NORM.quantiles, geneInfo and",
               " Des (which is the translation table from TCGA supplied GeneSymbols to updated Ensembl hgnc symbols).", "\n",
               "-----------------------------------------------------------------------------------------------------------")
  cat(msg)
  cat(msg, file= Log_file,sep = "\n",append=TRUE)
}
end.time.all = Sys.time ()
time = substring(as.character(capture.output(round(end.time.all - start.time.all, 2))),20,100)
msg = paste0("Normalization time for all cancertypes: ",time, "\n", 
             "---------------------------------------------------------------")
cat(msg)   
cat(msg, file= Log_file,sep = "\n",append=TRUE)

       