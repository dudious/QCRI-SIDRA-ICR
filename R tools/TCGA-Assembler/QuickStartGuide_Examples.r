#############################################################################
  
  # TCGA-Assembler : An open-source R program for downloading, processing and analyzing public TCGA data.
  # Copyright (C) <2014>  <Yitan Zhu>
  # This file is part of TCGA-Assembler.

  #  TCGA-Assembler is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.

  #  TCGA-Assembler is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.

  #  You should have received a copy of the GNU General Public License
  #  along with TCGA-Assembler.  If not, see <http://www.gnu.org/licenses/>.
  ############################################################################  


#####################################  Part I: Acquire Data  #########################################

# Clear workspace
rm(list = ls());

# Load module A functions.
source("Module_A.r");

# Download level-3 miRNA-seq data of six rectum adenocarcinoma (READ) samples
miRNASeqRawData = DownloadmiRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", 
  saveFolderName = "./QuickStartGuide_Results/RawData/", cancerType = "READ", 
  assayPlatform = "miRNASeq", inputPatientIDs = c("TCGA-EI-6884-01", 
  "TCGA-DC-5869-01", "TCGA-G5-6572-01", "TCGA-F5-6812-01", "TCGA-AF-2689-11", "TCGA-AF-2691-11")); 

# Download level-3 DNA copy number data of six READ samples
CNARawData = DownloadCNAData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", 
  saveFolderName = "./QuickStartGuide_Results/RawData/", cancerType = "READ", 
  assayPlatform = "genome_wide_snp_6", inputPatientIDs = c("TCGA-EI-6884-01",
  "TCGA-DC-5869-01", "TCGA-G5-6572-01", "TCGA-F5-6812-01", "TCGA-AF-2692-10", "TCGA-AG-4021-10")); 

# Download level-3 RNASeqV2 gene expression and exon expression data of six READ samples
RNASeqRawData = DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
  "./QuickStartGuide_Results/RawData/", cancerType = "READ", assayPlatform = "RNASeqV2", 
  dataType = c("rsem.genes.normalized_results", "exon_quantification"), inputPatientIDs = 
  c("TCGA-EI-6884-01", "TCGA-DC-5869-01", "TCGA-G5-6572-01", "TCGA-F5-6812-01", "TCGA-AG-3732-11", 
  "TCGA-AG-3742-11")); 

# Download level-3 HumanMethylation27 data of six READ samples
Methylation27RawData = DownloadMethylationData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
  "./QuickStartGuide_Results/RawData/", cancerType = "READ", assayPlatform = "humanmethylation27", 
  inputPatientIDs = c("TCGA-AG-3583-01", "TCGA-AG-A032-01", "TCGA-AF-2692-11", "TCGA-AG-4001-01", 
  "TCGA-AG-3608-01", "TCGA-AG-3574-01")); 

# Download level-3 HumanMethylation450 data of six READ samples
Methylation450RawData = DownloadMethylationData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
  "./QuickStartGuide_Results/RawData", cancerType = "READ", assayPlatform = "humanmethylation450", 
  inputPatientIDs = c("TCGA-EI-6884-01", "TCGA-DC-5869-01", "TCGA-G5-6572-01", "TCGA-F5-6812-01", 
  "TCGA-AG-A01W-11", "TCGA-AG-3731-11")); 

# Download level-3 RPPA protein expression data of six READ samples
RPPARawData = DownloadRPPAData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
  "./QuickStartGuide_Results/RawData", cancerType = "READ", assayPlatform = "mda_rppa_core", 
  inputPatientIDs = c("TCGA-EI-6884-01", "TCGA-DC-5869-01", "TCGA-G5-6572-01", "TCGA-F5-6812-01", 
  "TCGA-AG-3582-01", "TCGA-AG-4001-01"));  

# Download de-identified clinical information of READ patients
DownloadClinicalData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = 
  "./QuickStartGuide_Results/RawData", cancerType = "READ", clinicalDataType = c("patient", "drug", "follow_up"));



#####################################  Part II: Basic Data Processing  #########################################

# Clear workspace
rm(list = ls());

# Load module B functions.
source("Module_B.r");

# Process level-3 miRNA-seq data acquired by module A. 
miRNASeqData = ProcessmiRNASeqData(inputFilePath = 
  "./QuickStartGuide_Results/RawData/READ__bcgsc.ca__illuminahiseq_mirnaseq__GRCh37__Jul-08-2014.txt", 
  outputFileName = "READ__illuminahiseq_mirnaseq", outputFileFolder = 
  "./QuickStartGuide_Results/BasicProcessingResult");

# Process copy number data acquired by module A.
READ.GeneLevel.CNA = ProcessCNAData(inputFilePath = 
  "./QuickStartGuide_Results/RawData/READ__broad.mit.edu__genome_wide_snp_6__hg19__Jul-08-2014.txt", 
  outputFileName = "READ__genome_wide_snp_6__GeneLevelCNA", outputFileFolder = 
  "./QuickStartGuide_Results/BasicProcessingResult", refGenomeFile = "./SupportingFiles/Hg19GenePosition.txt");

# Process level-3 RNA-seq normalized gene expression data acquired by module A.
GeneExpData = ProcessRNASeqData(inputFilePath = 
  "./QuickStartGuide_Results/RawData/READ__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.normalized_results__Jul-08-2014.txt", 
  outputFileName = "READ__illuminahiseq_rnaseqv2__GeneExp", outputFileFolder = 
  "./QuickStartGuide_Results/BasicProcessingResult", dataType = "GeneExp", verType = "RNASeqV2");

# Process level-3 RNA-seq exon expression data acquired by module A.
ExonExpData = ProcessRNASeqData(inputFilePath = 
  "./QuickStartGuide_Results/RawData/READ__unc.edu__illuminahiseq_rnaseqv2__exon_quantification__Jul-08-2014.txt", 
  outputFileName = "READ__illuminahiseq_rnaseqv2__ExonExp", outputFileFolder = 
  "./QuickStartGuide_Results/BasicProcessingResult", dataType = "ExonExp", verType = "RNASeqV2");

# Process level-3 HumanMethylation27 data acquired by module A.
Methylation27Data = ProcessMethylation27Data(inputFilePath = 
  "./QuickStartGuide_Results/RawData/READ__jhu-usc.edu__humanmethylation27__Jul-08-2014.txt", 
  outputFileName = "READ__humanmethylation27", outputFileFolder = "./QuickStartGuide_Results/BasicProcessingResult");

# Process level-3 HumanMethylation450 data acquired by module A.
Methylation450Data = ProcessMethylation450Data(inputFilePath = 
  "./QuickStartGuide_Results/RawData/READ__jhu-usc.edu__humanmethylation450__Jul-08-2014.txt", 
  outputFileName = "READ__humanmethylation450", outputFileFolder = "./QuickStartGuide_Results/BasicProcessingResult");

# Process level-3 RPPA protein expression data acquired by module A.
RPPAData = ProcessRPPADataWithGeneAnnotation(
  inputFilePath = "./QuickStartGuide_Results/RawData/READ__mdanderson.org__mda_rppa_core__Jul-08-2014.txt", 
  outputFileName = "READ__mda_rppa_core", outputFileFolder = "./QuickStartGuide_Results/BasicProcessingResult");



#####################################  Part III: Advanced Data Processing  #########################################

# Clear workspace
rm(list = ls());

# Load module B functions.
source("Module_B.r");

# Load READ methylation27 and methylation450 data.
load("./QuickStartGuide_Results/BasicProcessingResult/READ__humanmethylation27.rda");
Methylation27Data = list(Des = Des, Data = Data);
load("./QuickStartGuide_Results/BasicProcessingResult/READ__humanmethylation450.rda");
Methylation450Data = list(Des = Des, Data = Data);

# Merge READ methylation27 and methylation450 data.
Methylation27_450_Merged = MergeMethylationData(input1 = Methylation27Data, input2 = 
  Methylation450Data, outputFileName = "READ__humanmethylation27_450_merged", outputFileFolder = 
  "./QuickStartGuide_Results/AdvancedProcessingResult");

# Calculate an average methylation value of CpG sites in each gene's region 
Methylation450_OverallAverage = CalculateSingleValueMethylationData(input = Methylation450Data, 
  regionOption = "All", DHSOption = "Both", outputFileName = "READ__humanmethylation450__SingleValue", 
  outputFileFolder = "./QuickStartGuide_Results/AdvancedProcessingResult");

# Calculate an average methylation value of CpG sites within 1500 base pairs of 
# transcription start site (TSS) and hypersensitive to DNAse. 
Methylation450_TSS1500_DHS = CalculateSingleValueMethylationData(input = Methylation450Data, 
  regionOption = "TSS1500", DHSOption = "DHS", outputFileName = "READ__humanmethylation450__SingleValue", 
  outputFileFolder = "./QuickStartGuide_Results/AdvancedProcessingResult");

# Extract methylation450 data of primary solid tumors.
ExtractedData_TP = ExtractTissueSpecificSamples(inputData = Methylation450Data$Data, tissueType = "TP", 
  singleSampleFlag = FALSE);

# Extract methylation450 data of primary solid tumors and solid normal tissues.
ExtractedData_TP_NT = ExtractTissueSpecificSamples(inputData = Methylation450Data$Data, 
  tissueType = c("TP", "NT"), singleSampleFlag = TRUE);

# Load multi-modal data of READ sampels
load("./QuickStartGuide_Results/AdvancedProcessingResult/READ__humanmethylation450__SingleValue__All__Both.rda");
READ__humanmethylation450__SingleValue__All__Both = list(Des = Des, Data = Data, dataType = "Methylation")
load("./QuickStartGuide_Results/BasicProcessingResult/READ__mda_rppa_core.rda");
READ__mda_rppa_core = list(Des = Des, Data = Data, dataType = "ProteinExp");
load("./QuickStartGuide_Results/BasicProcessingResult/READ__illuminahiseq_rnaseqv2__GeneExp.rda");
READ__illuminahiseq_rnaseqv2__GeneExp = list(Des = Des, Data = Data, dataType = "GeneExp");
load("./QuickStartGuide_Results/BasicProcessingResult/READ__genome_wide_snp_6__GeneLevelCNA.rda");
READ__genome_wide_snp_6__GeneLevelCNA = list(Des = Des, Data = Data, dataType = "CNA");
load("./QuickStartGuide_Results/BasicProcessingResult/READ__illuminahiseq_mirnaseq__RPM.rda");
READ__illuminahiseq_mirnaseq__RPM = list(Des = Des, Data = Data, dataType = "miRNAExp");

# Put multi-modal data in a vector of data list objects to be inputted into the data combination function.
inputDataList = vector("list", 5);
inputDataList[[1]] = READ__illuminahiseq_rnaseqv2__GeneExp;
inputDataList[[2]] = READ__humanmethylation450__SingleValue__All__Both;
inputDataList[[3]] = READ__genome_wide_snp_6__GeneLevelCNA;
inputDataList[[4]] = READ__mda_rppa_core;
inputDataList[[5]] = READ__illuminahiseq_mirnaseq__RPM;

# Merge multi-platform data
MergedData = CombineMultiPlatformData(inputDataList = inputDataList);

# Write the combined multi-platform data to into a tab-delimited txt file
write.table(cbind(MergedData$Des, MergedData$Data), 
  file = "./QuickStartGuide_Results/AdvancedProcessingResult/CombinedMultiPlatformData.txt", 
  quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE);