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

##############################################################################

# TCGA-Assembler Version 1.0.3 Module B

##############################################################################

library(RCurl);
library(httr);
library(stringr);
library(HGNChelper);



######################### Main Functions of Module B #############################################################

CombineMultiPlatformData<-function(inputDataList, combineStyle = "Intersect")
{
  options(warn=-1);
  
  for (i in 1:length(inputDataList))
  {
    
    # Keep only one sample for a tissue type of a patient. Usually, there is only one sample of a tissue type of a patient existing in data.
    inputDataList[[i]]$Data = ToPatientData(inputDataList[[i]]$Data);
    if (i == 1)
    {
      sample = colnames(inputDataList[[i]]$Data);
    }
    # For each genomic feature, keep only one row of data.
    Result = CombineRedundantFeature(Data = inputDataList[[i]]$Data, Des = inputDataList[[i]]$Des);
    
    inputDataList[[i]]$Des = Result$Des;
    inputDataList[[i]]$Data = Result$Data;
    if (combineStyle == "Intersect")
    {
      sample = sort(intersect(sample, colnames(inputDataList[[i]]$Data)));
    }
    if (combineStyle == "Union")
    {
      sample = sort(union(sample, colnames(inputDataList[[i]]$Data)));
    }    
    if (dim(inputDataList[[i]]$Des)[2] == 1)
    {
      inputDataList[[i]]$Des = cbind(inputDataList[[i]]$Des, Description = cbind(rep("", dim(inputDataList[[i]]$Des)[1])));
    }
    else
    {
      if ((dim(inputDataList[[i]]$Des)[2] == 3) && (inputDataList[[i]]$dataType == "CNA"))
      {
        inputDataList[[i]]$Des = cbind(inputDataList[[i]]$Des[, 1, drop = FALSE], 
        Description = paste(inputDataList[[i]]$Des[, 2], inputDataList[[i]]$Des[, 3], sep = ""));
      }
    }

    if (inputDataList[[i]]$dataType == "miRNAExp")
    {
      for (kk in 1:dim(inputDataList[[i]]$Des)[1])
      {
        inputDataList[[i]]$Des[kk, 1] = paste(substr(inputDataList[[i]]$Des[kk, 1], 5, 7), substr(inputDataList[[i]]$Des[kk, 1], 9, 100), sep = "");
      }
      inputDataList[[i]]$Des[, 1] = toupper(inputDataList[[i]]$Des[, 1]);
    }

    # for combining methylation data at CpG site level.
    if (inputDataList[[i]]$dataType == "Methylation")
    {    
      if (sum(toupper(c("REF", "GeneSymbol", "ChromosomeID", "CoordinateID")) %in% toupper(colnames(inputDataList[[i]]$Des))) == 4)
      {
        inputDataList[[i]]$Des = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
        Description = paste(inputDataList[[i]]$Des[, "REF"], inputDataList[[i]]$Des[, "ChromosomeID"], 
        inputDataList[[i]]$Des[, "CoordinateID"], sep = "|"));
        inputDataList[[i]]$Des[, "Description"] = gsub(" ", "", inputDataList[[i]]$Des[, "Description"]);
      }
    }
    
    inputDataList[[i]]$Des = switch(inputDataList[[i]]$dataType,
    GeneExp = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
    Platform = rep("GE", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]),
    ProteinExp = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
    Platform = rep("PE", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]),                                    
    Methylation = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
    Platform = rep("ME", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]),                                    
    CNA = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
    Platform = rep("CN", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]), 
    miRNAExp = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE],
    Platform = rep("miRExp", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]));

    # check NA gene. which(inputDataList[[i]]$Des[, "GeneSymbol"] != "NA")
    ID = intersect(which(!is.na(inputDataList[[i]]$Des[, "GeneSymbol"])), 
                   which(inputDataList[[i]]$Des[, "GeneSymbol"] != "?"));
    inputDataList[[i]]$Des = inputDataList[[i]]$Des[ID, , drop = FALSE];
    inputDataList[[i]]$Data = inputDataList[[i]]$Data[ID, , drop = FALSE];    

  }
  
  for (i in 1:length(inputDataList))
  {
    # Get the data of samples that should be kept
    if (combineStyle == "Intersect")
    {
      inputDataList[[i]]$Data = inputDataList[[i]]$Data[, sample, drop = FALSE];  
    }
    if (combineStyle == "Union")
    {
      tempData = matrix(NA, dim(inputDataList[[i]]$Data)[1], length(sample));
      colnames(tempData) = sample;
      tempData[, colnames(inputDataList[[i]]$Data)] = inputDataList[[i]]$Data;
      inputDataList[[i]]$Data = tempData;  
    }
    rownames(inputDataList[[i]]$Data) = NULL;
    rownames(inputDataList[[i]]$Des) = NULL;  
  }
  
  # Combine the datasets into matrix format
  Data = inputDataList[[1]]$Data;
  Des = inputDataList[[1]]$Des;  
  for (i in 2:length(inputDataList))
  {
    Data = rbind(Data, inputDataList[[i]]$Data);
    Des = rbind(Des, inputDataList[[i]]$Des);    
  } 

  OrderID = order(as.character(Des[, "GeneSymbol"]), as.character(Des[, "Platform"]), 
                  as.character(Des[, "Description"]), na.last = TRUE, decreasing = FALSE);
  Data = Data[OrderID, , drop = FALSE];
  Des = Des[OrderID, , drop = FALSE];
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  Result = list(Des = Des, Data = Data);

  options(warn=0);
  
  Result;
}



ExtractTissueSpecificSamples<-function(inputData, tissueType, singleSampleFlag, sampleTypeFile = "./SupportingFiles/TCGASampleType.txt")
{
  options(warn=-1);
  
  TCGABarcode = colnames(inputData);
  TCGABarcode = sapply(strsplit(TCGABarcode, split = "-"), function(x){
    Str = paste(x[1], x[2], x[3], substr(x[4], 1, 2), sep = "-");
    Str;
  });
  if (singleSampleFlag == TRUE)
  {
    DuplicatedLabel = duplicated(TCGABarcode);
    ID = which(DuplicatedLabel == FALSE);
    inputData = inputData[, ID, drop = FALSE];
    TCGABarcode = TCGABarcode[ID];
  }
  TCGADataLabel = sapply(strsplit(TCGABarcode, split = "-"), function(x)x[4]);
  TCGADataLabelValue = as.numeric(substr(TCGADataLabel, 1, 2));
  
  sampleType = read.table(file = sampleTypeFile, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE);  
  Code = sampleType[which(sampleType[, "Options"] %in% tissueType), "Code"];
  if (length(Code) == 0)
  {
    writeLines("Error: no valid tissue or cell type specified");
    inputData = NULL;
  }
  else
  {
    ID = which(TCGADataLabelValue %in% Code);
    if (length(ID) == 0)
    {
      writeLines(paste("There is no sample of ", paste(tissueType, collapse = ", "), sep = ""));
      inputData = NULL;
    }
    else
    {
      inputData = inputData[, ID, drop = FALSE];
      writeLines(paste(length(ID), " samples of ", paste(tissueType, collapse = ", "), " are extracted.", sep = ""));
    }
  }
  
  options(warn=0);
  
  inputData;
}



CalculateSingleValueMethylationData<-function(input, regionOption, DHSOption, outputFileName, outputFileFolder, chipAnnotationFile = "./SupportingFiles/MethylationChipAnnotation.rda")
{

  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  if (sum(! regionOption %in% c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "All")) > 0)
  {
    stop("Error: invalid regionOption.");
  }
  
  if (sum(! DHSOption %in% c("DHS", "notDHS", "Both")) > 0)
  {
    stop("Error: invalid DHSOption");
  }
    
  regionOptionO = regionOption;
  if (regionOptionO == "TSS1500")
  {
    regionOption = c(regionOption, "TSS200");
  }

  if ((regionOptionO != "All") || (DHSOption != "Both"))
  {
    load(chipAnnotationFile);
    ID = 1:dim(MethylAnno)[1];
    if (regionOptionO != "All")
    {
      ID = sort(unique(which(MethylAnno[, "UCSC_RefGene_Group"] %in% regionOption)));
    }
    if (DHSOption == "DHS")
    {
      ID = sort(intersect(ID, which(MethylAnno[, "DHS"] == TRUE)));
    }
    else
    {
      if (DHSOption == "notDHS")
      {
        ID = sort(intersect(ID, setdiff(1:dim(MethylAnno)[1], which(MethylAnno[, "DHS"] == TRUE))));
      }
    }
    MethylAnno = MethylAnno[ID, , drop = FALSE];
    ID = which(input$Des[, "REF"] %in% MethylAnno[, "IlmnID"]);
    input$Des = input$Des[ID, , drop = FALSE];    
    input$Data = input$Data[ID, , drop = FALSE];     
  }

  # Check whether conflict with NA gene
  ID = which(!is.na(input$Des[, "GeneSymbol"]));
  input$Des = input$Des[ID, , drop = FALSE];    
  input$Data = input$Data[ID, , drop = FALSE];  
  UniGene = unique(input$Des[, "GeneSymbol"]);
  Data = matrix(NA, length(UniGene), dim(input$Data)[2]);
  for (i in 1:length(UniGene))
  {
    if ((i %% 2000 == 0) || (i == length(UniGene)))
    {
      writeLines(paste("Calculation ", round(i/length(UniGene)*100, digits = 2), "% done.", sep = ""));      
    }
    Data[i, ] = colMeans(input$Data[which(input$Des[, "GeneSymbol"] == UniGene[i]), , drop = FALSE], na.rm = TRUE);
  }
  colnames(Data) = colnames(input$Data);
  
  Des = cbind(GeneSymbol = UniGene, SingleValueType = rep(paste(regionOptionO, DHSOption, sep = "|"), length(UniGene)));
  
  # draw and save a box plot of data
  png(filename = paste(outputFileFolder, "/", outputFileName, "__", regionOptionO, "__", DHSOption, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();
  
  # save data files.
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, "__", regionOptionO, "__", DHSOption, ".rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, "__", regionOptionO, "__", DHSOption, ".txt", sep = ""), 
              quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE);  
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des));
}



MergeMethylationData<-function(input1, input2, outputFileName, outputFileFolder)
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  input1_O = input1;
  input2_O = input2;  
  rownames(input1$Des) = paste(input1$Des[, "REF"], input1$Des[, "GeneSymbol"], sep = "||");
  rownames(input1$Data) = paste(input1$Des[, "REF"], input1$Des[, "GeneSymbol"], sep = "||");  
  rownames(input2$Des) = paste(input2$Des[, "REF"], input2$Des[, "GeneSymbol"], sep = "||");  
  rownames(input2$Data) = paste(input2$Des[, "REF"], input2$Des[, "GeneSymbol"], sep = "||");   
  CommonID = intersect(rownames(input1$Des), rownames(input2$Des));
  input1$Des = input1$Des[CommonID, , drop = FALSE];
  input1$Data = input1$Data[CommonID, , drop = FALSE];
  input2$Des = input2$Des[CommonID, , drop = FALSE];
  input2$Data = input2$Data[CommonID, , drop = FALSE];  
  if (dim(input1_O$Des)[1] >= dim(input2_O$Des)[1])
  {
    Des = input1$Des;
  }
  else
  {
    Des = input2$Des;
  }
  Data = cbind(input1$Data, input2$Data);
  
  # draw and save a box plot of combined data before normalization
  png(filename = paste(outputFileFolder, "/", outputFileName, "__BeforeNormalizationBoxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();

  Data = normalize.quantiles(Data);
  
  # draw and save a box plot of combined data after normalization
  png(filename = paste(outputFileFolder, "/", outputFileName, "__AfterNormalizationBoxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();
    
  # save data files.
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, ".rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE);  

  options(warn=0);
  
  return(list(Data = Data, Des = Des));
}



ProcessCNAData<-function(inputFilePath, outputFileName, outputFileFolder, refGenomeFile)
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  # Load reference genomes, which give the genomic locations of genes
  RefGenome = read.table(file = refGenomeFile, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE);
  NumGene = dim(RefGenome)[1];
  mode(RefGenome[, "txStarts"]) = "numeric";
  mode(RefGenome[, "txEnds"]) = "numeric";

  # Load copy number data
  InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);
  InData[, "Sample"] = toupper(InData[, "Sample"]);
  TCGAID = InData[, "Sample"];
  UniqueTCGAID = unique(TCGAID);  
  InData[, "Chromosome"] = as.character(InData[, "Chromosome"]);
  InData[, "Chromosome"] = paste("CHR", InData[, "Chromosome"], sep = "");
  IDID = which(InData[, "Chromosome"] == "CHR23");
  InData[IDID, "Chromosome"] = "CHRX";
  IDID = which(InData[, "Chromosome"] == "CHR24");
  InData[IDID, "Chromosome"] = "CHRY";  
  InData[, "Start"] = as.numeric(InData[, "Start"]);
  InData[, "End"] = as.numeric(InData[, "End"]);
  InData[, "Segment_Mean"] = as.numeric(InData[, "Segment_Mean"]);
  
  # Calculate gene-level copy number
  Data = matrix(NA, NumGene, length(UniqueTCGAID));
  for (i in 1:length(UniqueTCGAID))
  {
    IDi = which(InData[, "Sample"] == UniqueTCGAID[i]);
    InDatai = InData[IDi, c("Chromosome", "Start", "End", "Segment_Mean"), drop = FALSE];
    for (j in 1:NumGene)
    {
      IDj = which(InDatai[, "Chromosome"] == RefGenome[j, "Chromosome"]);
      if (length(IDj) > 0)
      {
        StartP = InDatai[IDj, "Start"];
        EndP = InDatai[IDj, "End"];
        Ratio = InDatai[IDj, "Segment_Mean"]; 
        Overlap = ifelse((EndP >= RefGenome[j, "txStarts"]) & (RefGenome[j, "txEnds"] >= StartP), pmin(EndP, RefGenome[j, "txEnds"]) - pmax(StartP, RefGenome[j, "txStarts"]) + 1, 0);
        IDNonNA = which(!is.na(Ratio));
        Ratio = Ratio[IDNonNA];
        Overlap = Overlap[IDNonNA];
        if (sum(Overlap, na.rm = TRUE) > 0)
        {
          Data[j, i] = sum(Overlap * Ratio, na.rm = TRUE)/sum(Overlap, na.rm = TRUE);  
        }
      }
    }
    if ((i %% 20 == 0) & (i < length(UniqueTCGAID)))
    {
      writeLines(paste("Calculating gene copy number , ", round(i/length(UniqueTCGAID)*100, digits = 2), "% done.", sep = ""));
    }
  }
  writeLines("Calculating gene copy number, 100% done.");  
  colnames(Data) = UniqueTCGAID;
  Des = as.matrix(RefGenome[, 1:3, drop = FALSE]);
  
  #Check and Correct gene symbol
  Des = CheckGeneSymbol(Des);
  
  # Draw and save a box plot
  png(filename = paste(outputFileFolder, "/", outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();
  
  # save data files.
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, ".rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE);  
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des));
}



ProcessRNASeqData<-function(inputFilePath, outputFileName, outputFileFolder, dataType, verType)
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  # Read in data.
  InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", stringsAsFactors = FALSE, quote = "", check.names = FALSE);
  InData = InData[2:dim(InData)[1], , drop = FALSE];  
  if ((dataType == "GeneExp") & (verType == "RNASeqV1"))
  {
    REF = InData[, 1];
    RPKM = as.matrix(InData[, seq(4, dim(InData)[2], 3), drop = FALSE]);
    mode(RPKM) = "numeric";
    GeneSymbol = sapply(strsplit(REF, split = "\\|"), function(x)x[1]);
    EntrezID = sapply(strsplit(REF, split = "\\|"), function(x)x[2]);
    Des = cbind(GeneSymbol = GeneSymbol, EntrezID = EntrezID);
    Data = RPKM;
    #Check and Correct gene symbol
    Des = CheckGeneSymbol(Des);    
  }
  if (dataType == "ExonExp")
  {
    REF = InData[, 1];
    RPKM = as.matrix(InData[, seq(4, dim(InData)[2], 3), drop = FALSE]);
    mode(RPKM) = "numeric";
    Des = cbind(ExonID = REF);
    Data = RPKM;    
  }
  if ((dataType == "GeneExp") & (verType == "RNASeqV2"))
  {
    REF = InData[, 1];
    NormalizedCount = as.matrix(InData[, 2:dim(InData)[2], drop = FALSE]);
    mode(NormalizedCount) = "numeric";
    GeneSymbol = sapply(strsplit(REF, split = "\\|"), function(x)x[1]);
    EntrezID = sapply(strsplit(REF, split = "\\|"), function(x)x[2]);
    Des = cbind(GeneSymbol = GeneSymbol, EntrezID = EntrezID);
    Data = NormalizedCount;    
    #Check and Correct gene symbol
    Des = CheckGeneSymbol(Des);     
  }
  if (verType == "Microarray")
  {
    dataType = "GeneExp";
    REF = InData[, 1];
    Data = as.matrix(InData[, 2:dim(InData)[2], drop = FALSE]);
    mode(Data) = "numeric";
    Des = cbind(GeneSymbol = REF, EntrezID = rep("", length(REF)));
    #Check and Correct gene symbol
    Des = CheckGeneSymbol(Des);     
  }
  
  if (dataType == "GeneExp")
  {
    # Draw and save a box plot
    png(filename = paste(outputFileFolder, "/", outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
    par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
    nonNaID = which(!is.na(Data));
    if ((sum(Data[nonNaID]<0)==0) && (max(Data[nonNaID])>50))
    {
      boxplot(log2(Data), main = "Boxplot drawn based on log2 tranformed data.");
    } else {
      boxplot(Data);
    }
    dev.off(); 
  }
  
  # save data files
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, ".rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des)); 
}



ProcessmiRNASeqData<-function(inputFilePath, outputFileName, outputFileFolder, fileSource = "TCGA-Assembler")
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  # Read in data.
  InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);
  InData = InData[2:dim(InData)[1], , drop = FALSE];
  
  # divide read cout data and RPM data
  REF = InData[, 1];  
  if (fileSource == "Firehose")
  {
    Count = as.matrix(InData[, seq(2, dim(InData)[2], 3), drop = FALSE]);
    RPM = as.matrix(InData[, seq(3, dim(InData)[2], 3), drop = FALSE]);
  }
  if (fileSource == "TCGA-Assembler")
  {
    Count = as.matrix(InData[, seq(2, dim(InData)[2], 2), drop = FALSE]);
    RPM = as.matrix(InData[, seq(3, dim(InData)[2], 2), drop = FALSE]);    
  }
#  rownames(Count) = REF;  
  mode(Count) = "numeric";
#  rownames(RPM) = REF;
  mode(RPM) = "numeric";
  
#   #   Draw and save box plot
#   png(filename = paste(outputFileFolder, "/", outputFileName, "__ReadCount.boxplot.png", sep = ""), width = 30*dim(Count)[2]+300, height = 1500, units = "px");
#   par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
#   boxplot(log2(Count));
#   dev.off();
#   png(filename = paste(outputFileFolder, "/", outputFileName, "__RPM.boxplot.png", sep = ""), width = 30*dim(RPM)[2]+300, height = 1500, units = "px");
#   par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
#   boxplot(log2(RPM));
#   dev.off();  
  
  #save data files
  Des = cbind(GeneSymbol = REF);
  Data = Count;
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, "__ReadCount.rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, "__ReadCount.txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 
  Data = RPM; 
  rownames(Data) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, "__RPM.rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, "__RPM.txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 

  options(warn=0);
  
  return(list(Data = Data, Des = Des)); 
}



ProcessMethylation450Data<-function(inputFilePath, outputFileName, outputFileFolder, fileSource = "TCGA-Assembler")
{

  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);  
  
  writeLines("Loading data.")
  
  if (fileSource == "Firehose")
  {
    # read in data
    GetMethylation450Data(inputFilePath, outputFileFolder, outputFileName);
    load(paste(outputFileFolder, "/", outputFileName, ".TempFiles/", outputFileName, ".TempData.rda", sep = ""));
    unlink(x = paste(outputFileFolder, "/", outputFileName, ".TempFiles", sep = ""), recursive = TRUE, force = TRUE);
    Methy450Des = as.matrix(Methy450Des);
  }
  
  if (fileSource == "TCGA-Assembler")
  {
    # read in data
    InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);  
    Methy450Des = as.matrix(InData[, 1:4, drop = FALSE]);
    colnames(Methy450Des) = c("REF", "GeneSymbol", "ChromosomeID", "CoordinateID");
    ID = which(Methy450Des[, "GeneSymbol"] == "DISCONTUNUED");
    Methy450Des[ID, "GeneSymbol"] = "";    
    Methy450Data = as.matrix(InData[, 5:dim(InData)[2], drop = FALSE]);
    mode(Methy450Data) = "numeric";
    rm("InData");
  }  
  
  # some probes have more than one gene symbols, collpase the data
  # AddDes = matrix("", dim(Methy450Des)[1], 4);
  AddDes = matrix("", 100000, 4);
  colnames(AddDes) = colnames(Methy450Des);
  # AddData = matrix(NA, dim(Methy450Data)[1], dim(Methy450Data)[2]);
  AddData = matrix(NA, 100000, dim(Methy450Data)[2]);
  colnames(AddData) = colnames(Methy450Data)
  AddIndex = 0;
  for (ProbeID in 1:dim(Methy450Des)[1])
  {
    GrepID = grep(";", Methy450Des[ProbeID, "GeneSymbol"]);
    if ((ProbeID %% 50000 == 0) || (ProbeID == dim(Methy450Des)[1]))
    {
      writeLines(paste("Process CpG sites corresponding to multiple genes, ", ceiling(ProbeID/dim(Methy450Des)[1]*100), "% done.", sep = ""));
    }
    if (length(GrepID) > 0)
    {
      TempStr = strsplit(Methy450Des[ProbeID, "GeneSymbol"], split = ";");
      Methy450Des[ProbeID, "GeneSymbol"] = TempStr[[1]][1];
      for (TempStri in 2:length(TempStr[[1]]))
      {
        TempAddDes = Methy450Des[ProbeID, , drop = FALSE];
        TempAddDes[1, "GeneSymbol"] = TempStr[[1]][TempStri];
        AddIndex = AddIndex + 1;
        AddDes[AddIndex, ] = TempAddDes;
        AddData[AddIndex, ] = Methy450Data[ProbeID, ];
      }
    }
  }

  # writeLines(paste("Number of added rows is ", AddIndex, sep = ""));

  Methy450Des = rbind(Methy450Des, AddDes[1:AddIndex, , drop = FALSE]);
  Methy450Data = rbind(Methy450Data, AddData[1:AddIndex, , drop = FALSE]);
  Methy450Des[, "ChromosomeID"] = gsub(" ", "", Methy450Des[, "ChromosomeID"]);
  Methy450Des[, "CoordinateID"] = gsub(" ", "", Methy450Des[, "CoordinateID"]);  
  OrderID = order(as.numeric(Methy450Des[, "ChromosomeID"]), as.numeric(Methy450Des[, "CoordinateID"]), Methy450Des[, "GeneSymbol"], na.last = TRUE, decreasing = FALSE);
  Des = Methy450Des[OrderID, , drop = FALSE];
  Data = Methy450Data[OrderID, , drop = FALSE];
  
  #Check and Correct gene symbol
  Des = CheckGeneSymbol(Des);
  
  # Draw and save box plot
  png(filename = paste(outputFileFolder, "/", outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();
  
  # save output data files
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Data, Des, file = paste(outputFileFolder, "/", outputFileName, ".rda",  sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), 
              quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE);    
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des));  
}



CheckGeneSymbol<-function(Des)
{
  data(hgnc.table);
  hgnc.table = rbind(hgnc.table, c("13-SEP", "SEPT7P2"));
  rID = intersect(which(toupper(hgnc.table[, "Symbol"]) == toupper("NCRNA00185")), 
                  which(toupper(hgnc.table[, "Approved.Symbol"]) == toupper("TTTY14")));
  hgnc.table = hgnc.table[sort(setdiff(1:dim(hgnc.table)[1], rID)), , drop = FALSE];
  hgnc.table = hgnc.table[which(!is.na(hgnc.table[, "Approved.Symbol"])), , drop = FALSE];
  regex = "[0-9]\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\.[0-9][0-9]E\\+[[0-9][0-9]";
  MonthID = grep(pattern = regex, hgnc.table[, 1], ignore.case = TRUE);
  MonthMappingTable = hgnc.table[MonthID, , drop = FALSE];
  Des[, "GeneSymbol"] = sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", Des[, "GeneSymbol"]);
  ID = intersect(which(!Des[, "GeneSymbol"] %in% hgnc.table[, "Approved.Symbol"]), which(Des[, "GeneSymbol"] %in% hgnc.table[, "Symbol"]));
  DesOrg = Des;
  if (length(ID) > 0)
  {
    Des[ID, "GeneSymbol"] = sapply(Des[ID, "GeneSymbol"], function(x)paste(hgnc.table[hgnc.table[, "Symbol"] == x, "Approved.Symbol"], collapse = "___"));
  }
#  writeLines("Changed genes are");
#  print(cbind(DesOrg[ID, "GeneSymbol"], Des[ID, "GeneSymbol"]));  
  ID = intersect(which(!Des[, "GeneSymbol"] %in% hgnc.table[, "Approved.Symbol"]), which(toupper(Des[, "GeneSymbol"]) %in% toupper(MonthMappingTable[, "Symbol"])));
  if (length(ID) > 0)
  {
    Des[ID, "GeneSymbol"] = sapply(Des[ID, "GeneSymbol"], function(x)paste(MonthMappingTable[toupper(MonthMappingTable[, "Symbol"]) == toupper(x), "Approved.Symbol"], collapse = "___"));
  }
#  writeLines("Changed genes are");
#  print(cbind(DesOrg[ID, "GeneSymbol"], Des[ID, "GeneSymbol"]));  
  Des;
}



ProcessRPPADataWithGeneAnnotation<-function(inputFilePath, outputFileName, outputFileFolder)
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  # Read in data
  InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);
  OutData = matrix(NA, 0, dim(InData)[2]-1);
  
  # Split the gene symbol and protein antibody into two separate columns
  # Duplicate rows for proteins coded by multiple genes.
  Gene = c();
  Antibody = c();
  for (i in 1:dim(InData)[1])
  {
    GeneAntibodyStr = InData[i, 1];
    GeneStr = strsplit(GeneAntibodyStr, split = "\\|")[[1]][1];
    AntibodyStr = strsplit(GeneAntibodyStr, split = "\\|")[[1]][2];
    GeneStr = unlist(strsplit(GeneStr, split = " "));
    Gene = c(Gene, GeneStr);
    Antibody = c(Antibody, rep(AntibodyStr, length = length(GeneStr)));
    OutData = rbind(OutData, InData[rep(i, length = length(GeneStr)), 2:dim(InData)[2], drop = FALSE]);
  }
  
  colnames(OutData) = colnames(InData)[2:dim(InData)[2]];
  Des = cbind(GeneSymbol = Gene, ProteinAntibody = Antibody);
  # check gene symbol
  Des = CheckGeneSymbol(Des);
  Data = as.matrix(OutData);
  mode(Data) = "numeric";
  
  # save output data files
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Data, Des, file = paste(outputFileFolder, "/", outputFileName, ".rda",  sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), 
              quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE);
  
  # draw and save boxplot
  png(filename = paste(outputFileFolder, "/", outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();
   
  options(warn=0);
  
  return(list(Data = Data, Des = Des));
}



ProcessMethylation27Data<-function(inputFilePath, outputFileName, outputFileFolder, fileSource  = "TCGA-Assembler")
{
  
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  if (fileSource == "TCGA-Assembler")
  {  
    # read in data    
    InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);  
    
    Methy27Des = as.matrix(InData[, 1:4, drop = FALSE]);
    colnames(Methy27Des) = c("REF", "GeneSymbol", "ChromosomeID", "CoordinateID");
    ID = which(Methy27Des[, "GeneSymbol"] == "DISCONTUNUED");
    Methy27Des[ID, "GeneSymbol"] = "";    
    Methy27Data = as.matrix(InData[, 5:dim(InData)[2], drop = FALSE]);
    mode(Methy27Data) = "numeric";
  }
  
  if (fileSource == "Firehose")
  {  
    # read in data    
    InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);  
    InData = InData[2:dim(InData)[1], , drop = FALSE];
    
    # remove redudant columns of CpG site descriptions    
    GeneSymbol = as.matrix(InData[, seq(3, dim(InData)[2], 4), drop = FALSE]);
    ID = which(GeneSymbol == "DISCONTUNUED");
    GeneSymbol[ID] = "";
    GeneSymbol = t(unique(t(GeneSymbol)));  
    ChromID = as.matrix(InData[, seq(4, dim(InData)[2], 4), drop = FALSE]);
    ChromID = t(unique(t(ChromID)));
    CoordID = as.matrix(InData[, seq(5, dim(InData)[2], 4), drop = FALSE]);
    CoordID = t(unique(t(CoordID)));
    Methy27Data = as.matrix(InData[, seq(2, dim(InData)[2], 4), drop = FALSE]);
    mode(Methy27Data) = "numeric";
    REF = InData[, 1];
    # Check whether the probe information is the same between samples
    if (dim(GeneSymbol)[2] != 1)
    {
      stop("Some probe's gene symbols are not consistent between samples");      
    }
    if (dim(ChromID)[2] != 1)
    {
      stop("Some probe's chromosome IDs are not consistent between samples");      
    }
    if (dim(CoordID)[2] != 1)
    {
      stop("Some probe's genmoe coordinates are not consistent between samples");      
    }
    ChromID = toupper(ChromID);
    Methy27Des = cbind(REF, GeneSymbol, ChromID, CoordID);
    colnames(Methy27Des) = c("REF", "GeneSymbol", "ChromosomeID", "CoordinateID");
  }
  
  # some probes have more than one gene symbols, collpase the data
  AddGeneDes = matrix("", dim(Methy27Des)[1]*3, 4);
  AddData = matrix(NA, dim(Methy27Data)[1]*3, dim(Methy27Data)[2]);
  AddIndex = 0;  
  for (ProbeID in 1:dim(Methy27Des)[1])
  {
    GrepID = grep(";", Methy27Des[ProbeID, "GeneSymbol"]);
    if (length(GrepID) > 0)
    {
      TempStr = strsplit(Methy27Des[ProbeID, "GeneSymbol"], split = ";");
      Methy27Des[ProbeID, "GeneSymbol"] = TempStr[[1]][1];
      for (TempStri in 2:length(TempStr[[1]]))
      {
        TempAddGeneDes = c(Methy27Des[ProbeID, "REF"], TempStr[[1]][TempStri], Methy27Des[ProbeID, "ChromosomeID"], Methy27Des[ProbeID, "CoordinateID"]);
        AddIndex = AddIndex + 1;
        AddGeneDes[AddIndex, ] = TempAddGeneDes;
        AddData[AddIndex, ] = Methy27Data[ProbeID, ];        
      }
    }
  }
  
  colnames(AddGeneDes) = colnames(Methy27Des);
  colnames(AddData) = colnames(Methy27Data);
  Methy27Des = rbind(Methy27Des, AddGeneDes[1:AddIndex, , drop = FALSE]);
  Methy27Data = rbind(Methy27Data, AddData[1:AddIndex, , drop = FALSE]);
  Methy27Des[, "ChromosomeID"] = gsub(" ", "", Methy27Des[, "ChromosomeID"]);
  Methy27Des[, "CoordinateID"] = gsub(" ", "", Methy27Des[, "CoordinateID"]);
  OrderID = order(as.numeric(Methy27Des[, "ChromosomeID"]), as.numeric(Methy27Des[, "CoordinateID"]), Methy27Des[, "GeneSymbol"], na.last = TRUE, decreasing = FALSE);
  Methy27Des = Methy27Des[OrderID, , drop = FALSE];
  Methy27Data = Methy27Data[OrderID, , drop = FALSE];
  Data = Methy27Data;
  Des = Methy27Des;
  # check gene symbol
  Des = CheckGeneSymbol(Des);
    
  # Draw box plot
  png(filename = paste(outputFileFolder, "/", outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
  par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  boxplot(Data);
  dev.off();
    
  # save output data files
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Data, Des, file = paste(outputFileFolder, "/", outputFileName, ".rda",  sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), 
              quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE);    
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des));  
}



######################### Auxiliary Functions of Module B #############################################################

ToPatientData<-function(TCGAData)
{
  TCGABarcode = colnames(TCGAData);
  TCGABarcode = sapply(strsplit(TCGABarcode, split = "-"), function(x){
    Str = paste(x[1], x[2], x[3], substr(x[4], 1, 2), sep = "-");
    Str;
  });
  DuplicatedLabel = duplicated(TCGABarcode);
  ID = which(DuplicatedLabel == FALSE);
  TCGAData = TCGAData[, ID, drop = FALSE];
  TCGABarcode = TCGABarcode[ID];
  colnames(TCGAData) = TCGABarcode;
  TCGAData;
}



CombineRedundantFeature <- function(Data, Des)
{
  RowName = Des[, 1];
  if (dim(Des)[2] > 1)
  {
    for (i in 2:dim(Des)[2])
    {
      RowName = paste(RowName, Des[, i], sep = "||");
    }
  }
  
  UniqueRowName = unique(RowName);
  if (length(RowName) == length(UniqueRowName))
  {
    Result = list(Data = Data, Des = Des);
  }
  else
  {
    IDKeep = c();
    for (i in 1:length(UniqueRowName))
    {
      IDi = which(RowName == UniqueRowName[i]);
      if (length(IDi) > 1)
      {
        Data[IDi[1], ] = apply(Data[IDi, , drop = FALSE], 2, mean, na.rm = TRUE);
      }
      IDKeep = c(IDKeep, IDi[1]);
    }
    IDKeep = sort(IDKeep);
    Result = list(Data = Data[IDKeep, , drop = FALSE], Des = Des[IDKeep, , drop = FALSE]);
  }
  Result;
}



normalize.quantiles <- function(Data)
{
  VectorData = as.vector(Data);
  VectorData = VectorData[!is.na(VectorData)];
  VectorData = sort(VectorData, decreasing = FALSE);
  NumVectorData = length(VectorData);
  NormData = matrix(NA, dim(Data)[1], dim(Data)[2]);
  
  for (i in 1:dim(Data)[2])
  {
    if ((i %% 10 == 0) || (i == dim(Data)[2]))
    {
      writeLines(paste("Normalizing data, ", round(i/dim(Data)[2]*100, digits = 2), "% done.", sep = ""));
    }
    IDNonNA = which(!is.na(Data[, i]));
    RankNonNA = rank(Data[IDNonNA, i], ties.method = "average"); 
    NumNonNA = length(IDNonNA);
    Sep = round(seq(1, NumVectorData, (NumVectorData-1)/NumNonNA));
    NormData[IDNonNA, i] = sapply(RankNonNA, function(x)
    {
      NumTie = sum(RankNonNA == x);
      mean(VectorData[(Sep[ceiling(x-NumTie/2)]+1):Sep[floor(x+NumTie/2)+1]]);
    })
  }
  rownames(NormData) = rownames(Data);
  colnames(NormData) = colnames(Data);
  NormData;
}



GetMethylation450Data<-function(InputFileName, OutputFileFolder, outputFileName)
{
  OriginalOutputFileFolder = OutputFileFolder;
  OutputFileFolder = paste(OutputFileFolder, "/", outputFileName, ".TempFiles/", sep = "");
  dir.create(OutputFileFolder, showWarnings = TRUE, recursive = FALSE, mode = "0777");
  
  FileHandle = file(InputFileName, "rt");
  on.exit(close(FileHandle));
  ReadIn = readLines(FileHandle, 1);
  ColumnName1 = strsplit(ReadIn, split = "\t")[[1]];
  ReadIn = readLines(FileHandle, 1);  
  ColumnName2 = strsplit(ReadIn, split = "\t")[[1]];  
  
  LineNum = 10000;
  Iter = 0;
  Flag = TRUE;
  while (Flag)
  {
    ReadIn = readLines(FileHandle, LineNum);
    if (length(ReadIn) > 0)
    {
      Iter = Iter + 1;
      if (Iter %% 5 == 0)
      {
        writeLines(paste("Load ", floor(Iter/49*100), "% of the humanmethylation450 data set.", sep = ""));
      }
      if (length(ReadIn) == 1)
      {
        TempMatrix = rbind(strsplit(ReadIn, split = "\t")[[1]]);
        colnames(TempMatrix) = ColumnName1;
        TempMethyData = rbind(TempMatrix[, seq(2, dim(TempMatrix)[2], 4), drop = FALSE]);      
        GeneSymbol = unique(TempMatrix[, seq(3, dim(TempMatrix)[2], 4)]);
        ChromID = unique(TempMatrix[, seq(4, dim(TempMatrix)[2], 4)]);
        CoordID = unique(TempMatrix[, seq(5, dim(TempMatrix)[2], 4)]);      
        REF = unique(TempMatrix[, 1]);
        ChromID = toupper(ChromID);
        # Check whether the probe information is the same between samples
        if (length(GeneSymbol) != 1)
        {
          stop("Some probe's gene symbols are not consistent between samples");      
        }
        if (length(ChromID) != 1)
        {
          stop("Some probe's chromosome IDs are not consistent between samples");      
        }
        if (length(CoordID) != 1)
        {
          stop("Some probe's genmoe coordinates are not consistent between samples");      
        }      
        TempDes = cbind(REF, GeneSymbol, ChromID, CoordID);
      }
      else
      {
        ReadInSep = strsplit(ReadIn, split = "\t");
        TempMatrix = matrix(NA, length(ReadIn), length(ColumnName1));
        for (i in 1:length(ReadInSep))
        {
          TempMatrix[i, ] = ReadInSep[[i]];
        }
        colnames(TempMatrix) = ColumnName1;
        GeneSymbol = TempMatrix[, seq(3, dim(TempMatrix)[2], 4), drop = FALSE];
        ChromID = TempMatrix[, seq(4, dim(TempMatrix)[2], 4), drop = FALSE];
        CoordID = TempMatrix[, seq(5, dim(TempMatrix)[2], 4), drop = FALSE];
        TempMethyData = as.matrix(TempMatrix[, seq(2, dim(TempMatrix)[2], 4), drop = FALSE]);
        GeneSymbol = t(unique(t(GeneSymbol)));
        ChromID = t(unique(t(ChromID)));
        CoordID = t(unique(t(CoordID)));
        REF = TempMatrix[, 1];
        ChromID = toupper(ChromID);
        # Check whether the probe information is the same between samples
        if (dim(GeneSymbol)[2] != 1)
        {
          stop("Some probe's gene symbols are not consistent between samples");      
        }
        if (dim(ChromID)[2] != 1)
        {
          stop("Some probe's chromosome IDs are not consistent between samples");      
        }
        if (dim(CoordID)[2] != 1)
        {
          stop("Some probe's genmoe coordinates are not consistent between samples");      
        }
        TempDes = cbind(REF, GeneSymbol, ChromID, CoordID);
      }
      colnames(TempDes) = c("REF", "GeneSymbol", "ChromosomeID", "CoordinateID");  
      mode(TempMethyData) = "numeric";
      TempMethyData = round(x = TempMethyData, digits = 4);
      save(TempMethyData, TempDes, file = paste(OutputFileFolder, "Methy450Data", Iter, ".rda", sep = ""));  
    }
    if (length(ReadIn) < LineNum)
    {
      Flag = FALSE;
    }
  }
  writeLines(paste("Load 100% of the humanmethylation450 data set.", sep = ""));
  
  rm(ReadIn, ReadInSep, TempMatrix, GeneSymbol, ChromID, CoordID, TempMethyData, REF, TempDes)
  load(paste(OutputFileFolder, "Methy450Data", 1, ".rda", sep = ""));
  Methy450Des = TempDes;
  Methy450Data = TempMethyData;
  for (i in 2:Iter)
  {
    load(paste(OutputFileFolder, "Methy450Data", i, ".rda", sep = ""));
    Methy450Des = rbind(Methy450Des, TempDes);
    Methy450Data = rbind(Methy450Data, TempMethyData);
  }
  
  save(Methy450Data, Methy450Des, file = paste(OutputFileFolder, outputFileName, ".TempData.rda", sep = ""));
}



######################### Check whether this is the most updated version of TCGA-Assembler #####################################################

VCwebContent = try(content(GET("http://health.bsd.uchicago.edu/yji/soft.html"), as = "text"), silent = TRUE);
if (class(VCwebContent) == "try-error")
{
  rm(VCwebContent);
} else {
  VCstartID = str_locate_all(string = VCwebContent, pattern = "LinkToCheckTCGA-AssebmlerVersionNumber");
  if (dim(VCstartID[[1]])[1] == 0)
  {
    rm(VCwebContent, VCstartID);    
  } else {
    VCstartID = VCstartID[[1]][1, "end"]+1;
    VCwebContent = substr(VCwebContent, VCstartID, nchar(VCwebContent));
    VCstartID = str_locate_all(string = VCwebContent, pattern = "href=\"")[[1]][1, "end"]+1;
    VCendID = str_locate_all(string = VCwebContent, pattern = "\">")[[1]][1, "start"]-1;
    VCurl = substr(VCwebContent, VCstartID, VCendID);
    VCwebContent = try(content(GET(VCurl), as = "text"), silent = TRUE);
    if (class(VCwebContent) == "try-error")
    {
      rm(VCstartID, VCwebContent, VCendID, VCurl);
    } else {
      VCstartID = str_locate_all(string = VCwebContent, pattern = "CheckVersionNumber1");
      if (dim(VCstartID[[1]])[1] == 0)
      {
        rm(VCstartID, VCwebContent, VCendID, VCurl);
      } else {
        VCstartID = VCstartID[[1]][1, "end"]+1;
        VCwebContent = substr(VCwebContent, VCstartID, nchar(VCwebContent));
        VCstartID = str_locate_all(string = VCwebContent, pattern = "\">")[[1]][1, "end"]+1;
        VCendID = str_locate_all(string = VCwebContent, pattern = "</span>")[[1]][1, "start"]-1;
        VCnewestVersionNum = substr(VCwebContent, VCstartID, VCendID);
        if (VCnewestVersionNum != "1.0.3")
        {
          writeLines("\n");
          writeLines("***************************************************************");
          writeLines("A new version of TCGA-Assembler is available!")
          writeLines(paste("Please download version ", VCnewestVersionNum, " at ", VCurl, sep = ""));
          writeLines("***************************************************************");      
          writeLines("\n");      
        }
        rm(VCstartID, VCwebContent, VCendID, VCnewestVersionNum, VCurl);
      }
    }
  }
}
