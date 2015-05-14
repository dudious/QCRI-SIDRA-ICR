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

# TCGA-Assembler Version 1.0.3 Module A

##############################################################################

library(RCurl);
library(httr);
library(stringr);



######################### Main Functions of Module A #####################################################

TraverseAllDirectories <- function(entryPoint, fileLabel)
{
  options(warn=-1);
  time1 = proc.time();
  
  #   get all the files/directories under the entryPoint
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines("Traverse data directories of TCGA DCC to gather URLs of TCGA data files");
  
  URLcontent = getTCGA_URL(entryPoint);
  file_url = URLcontent$file_url;
  dir_url = URLcontent$dir_url;
  
  counter = 1;
  while (counter<=length(dir_url))
  {
    URLcontent = getTCGA_URL(dir_url[counter]);
    file_url = c(file_url, URLcontent$file_url);
    dir_url = c(dir_url, URLcontent$dir_url);
    if ((counter %% 500) == 0)
    {
      writeLines(paste("Identified ", counter, " directories and ", length(file_url), " files.", sep = ""));
    }
    counter = counter + 1;
  }
  writeLines(paste("IN TOTAL, identified ", length(dir_url), " directories and ", length(file_url), " files.", sep = ""));

  #   generate the filename to store results
  date_string = date();
  date_string = strsplit(date_string, split = " ")[[1]];
  date_string = paste(date_string[2], date_string[3], date_string[5], sep = "-");
  filename = paste(fileLabel, '_', date_string, '.rda', sep = "");
  orderID = order(file_url);
  file_url = file_url[orderID];
  upper_file_url = toupper(file_url);
  save(file_url, dir_url, upper_file_url, file = filename);
  
  time = proc.time() - time1;
  writeLines(paste("Total elapsed time is ", round(time[3]/3600, digits = 3), " hours.", sep = ""));
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
}



DownloadClinicalData <- function(traverseResultFile, saveFolderName, cancerType, clinicalDataType, outputFileName = "")
{
  options(warn=-1);

  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download clinical information of ", cancerType, " patients.", sep = ""));
  writeLines("Load information of TCGA data files.");  
  load(traverseResultFile);
  
  Platform = "bio/clin";
  Institution = "nationwidechildrens.org";
  dir.create(path = saveFolderName, recursive = TRUE);
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }
  
  for (IDIn in 1:length(Institution))
  {
    for (IDPl in 1:length(Platform))
    {
      DirURL = paste("/", cancerType, "/bcr/", Institution[IDIn], "/", Platform[IDPl], "/", sep = "");
      ID_DirURL = grep(pattern = toupper(DirURL), x = upper_file_url, ignore.case = FALSE);
      ID_DirURL = ID_DirURL[grep(pattern = toupper("Level_2"), x = upper_file_url[ID_DirURL], ignore.case = FALSE)];
      for (IDFile in 1:length(clinicalDataType))
      {
        FileName = paste("clinical_", clinicalDataType[IDFile], sep = "");
        ind = ID_DirURL[grep(pattern = toupper(FileName), x = upper_file_url[ID_DirURL], ignore.case = FALSE)];
        FileName = unique(sapply(strsplit(file_url[ind], split = "/"), function(x) x[length(x)]));
        if (length(FileName) > 0)
        {
          for (FileNameIndex in 1:length(FileName))
          {        
            ind = ID_DirURL[grepEnd(pattern = toupper(FileName[FileNameIndex]), x = upper_file_url[ID_DirURL], ignore.case = FALSE)];
            if (length(ind) == 0)
            {
              next;
            }else{
              if (length(ind) > 1)
              {
                URL = GetNewestURL(AllURL = file_url[ind]);
              }else{
                URL = file_url[ind];
              }
            }          
            SaveFileName = paste(saveFolderName, "/", outputFileName, FileName[FileNameIndex], sep = "");
            downloadFile(url = URL, saveFileName = SaveFileName);
          }
        }
      }
    }
  }
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
}



DownloadBiospecimenData <- function(traverseResultFile, saveFolderName, cancerType, biospecimenDataType, outputFileName = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download biospecimen information of ", cancerType, " patients.", sep = ""));
  writeLines("Load information of TCGA data files.");  
  load(traverseResultFile);
  
  Platform = "bio/clin";
  Institution = "nationwidechildrens.org";
  dir.create(path = saveFolderName, recursive = TRUE);
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }
  
  for (IDIn in 1:length(Institution))
  {
    for (IDPl in 1:length(Platform))
    {
      DirURL = paste("/", cancerType, "/bcr/", Institution[IDIn], "/", Platform[IDPl], "/", sep = "");
      ID_DirURL = grep(pattern = toupper(DirURL), x = upper_file_url, ignore.case = FALSE);
      ID_DirURL = ID_DirURL[grep(pattern = toupper("Level_2"), x = upper_file_url[ID_DirURL], ignore.case = FALSE)];
      for (IDFile in 1:length(biospecimenDataType))
      {
        FileName = paste("biospecimen_", biospecimenDataType[IDFile], sep = "");
        ind = ID_DirURL[grep(pattern = toupper(FileName), x = upper_file_url[ID_DirURL], ignore.case = FALSE)];
        FileName = unique(sapply(strsplit(file_url[ind], split = "/"), function(x) x[length(x)]));
        if (length(FileName) > 0)
        {
          for (FileNameIndex in 1:length(FileName))
          {        
            ind = ID_DirURL[grepEnd(pattern = toupper(FileName[FileNameIndex]), x = upper_file_url[ID_DirURL], ignore.case = FALSE)];
            if (length(ind) == 0)
            {
              next;
            }else{
              if (length(ind) > 1)
              {
                URL = GetNewestURL(AllURL = file_url[ind]);
              }else{
                URL = file_url[ind];
              }
            }          
            SaveFileName = paste(saveFolderName, "/", outputFileName, FileName[FileNameIndex], sep = "");
            downloadFile(url = URL, saveFileName = SaveFileName);
          }
        }
      }
    }
  }
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
}



DownloadCNAData <- function(traverseResultFile, saveFolderName, cancerType, assayPlatform = "genome_wide_snp_6", tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download copy number data of ", cancerType, " patients.", sep = ""));
  
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  load(traverseResultFile);
  dir.create(path = saveFolderName, recursive = TRUE);
  SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/broad\\.mit\\.edu/", assayPlatform, "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
  IDLevel_3InFile_Url = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  SpecificName = paste(cancerType, "__broad.mit.edu__", assayPlatform, sep = "");
  
  # download and process Sample and Data Relationship Format (SDRF) file
  ind = SpecificID[grepEnd(pattern = toupper("\\.sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Program existed due to missing SDRF file."); 
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    return();
  }

  sdrf = toupper(downloadResult$data);
  level_3_filename_column = which(sdrf[1, ] == toupper("Derived Array Data File"));
  DataLevelColID = which(sdrf[1, ] == toupper("Comment [TCGA Data Level]"));
  DataLevelColID = DataLevelColID[(length(DataLevelColID)-length(level_3_filename_column)+1):length(DataLevelColID)];
  BarcodeColID = which(sdrf[1, ] == toupper("Comment [TCGA Barcode]"));
  colnames(sdrf) = sdrf[1, ];
  sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
  Level3_ID = 1:dim(sdrf)[1];
  for (DataLevelColIDIndex in 1:length(DataLevelColID))
  {
    Level3_ID = sort(intersect(Level3_ID, union(which(sdrf[, DataLevelColID[DataLevelColIDIndex]] == "LEVEL_3"),
                    which(sdrf[, DataLevelColID[DataLevelColIDIndex]] == "LEVEL 3"))), decreasing = FALSE);
  }
  if (length(Level3_ID) == 0)
  {
    writeLines("Error: there are no Level 3 data");
    return();
  }
  sdrf = unique(sdrf[Level3_ID, sort(c(BarcodeColID, level_3_filename_column, DataLevelColID), decreasing = FALSE), drop = FALSE]);  
  
  # If specific patient TCGA barcodes are inputted, only download data of the specified samples.
  if (!is.null(inputPatientIDs))
  {
    indInputPatientID = c();
    for (i in 1:length(inputPatientIDs))
    {
      indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
    }
    if (length(indInputPatientID) == 0)
    {
      writeLines("No Level 3 data for the inputted TCGA barcodes.");
      return();      
    }else{
      sdrf = sdrf[indInputPatientID, , drop = FALSE];
    }
  }
  
  # Download data of specified tissue
  if (!is.null(tissueType))
  {
    SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                       Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
    sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
  }
  
  if (dim(sdrf)[1] == 0)
  {
    writeLines("No available data.");
    return();
  }
  
  CopyNumberData = vector("list", 4); 
  names(CopyNumberData) = c("hg18", "hg19", "nocnv_hg18", "nocnv_hg19");
  for (i in 1:length(CopyNumberData))
  {
    CopyNumberData[[i]] = matrix("", 0, 6);
  }

  for (i in 1:dim(sdrf)[1])
  {
    time1 = proc.time();
    sample_TCGA_id = sdrf[i, 1];
    for (DataFileID in 1:floor(dim(sdrf)[2]/2))
    {
      DataFileName_ID = sdrf[i, DataFileID*2];
      DataFileType_ID = strsplit(DataFileName_ID, split = "\\.")[[1]];
      DataFileType_ID = DataFileType_ID[length(DataFileType_ID)-2];
      CellID = which(toupper(names(CopyNumberData)) == DataFileType_ID);
      
      ind = IDLevel_3InFile_Url[grepEnd(pattern = DataFileName_ID, x = upper_file_url[IDLevel_3InFile_Url], ignore.case = FALSE)];
      if (length(ind) == 0)
      {
        next;
      }else{
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        s = downloadResult$data;
        colnames(s) = s[1, ];
        s = s[2:dim(s)[1], , drop = FALSE];
        ChromosomeID = grep(pattern = "Chromosome", x = colnames(s), ignore.case = TRUE);
        StartID = grep(pattern = "Start", x = colnames(s), ignore.case = TRUE);
        EndID = grep(pattern = "End", x = colnames(s), ignore.case = TRUE);
        NumProbesID = grep(pattern = "Num_Probes", x = colnames(s), ignore.case = TRUE);
        SegmentMeanID = grep(pattern = "Segment_Mean", x = colnames(s), ignore.case = TRUE);
        IDX = which(toupper(s[, ChromosomeID]) == "X");
        s[IDX, ChromosomeID] = "23";
        IDY = which(toupper(s[, ChromosomeID]) == "Y");
        s[IDY, ChromosomeID] = "24";       
        CopyNumberData[[CellID]] = rbind(CopyNumberData[[CellID]], cbind(Sample = rep(sample_TCGA_id, dim(s)[1]),
        s[, c(ChromosomeID, StartID, EndID, NumProbesID, SegmentMeanID), drop = FALSE]));
      }
    }
    time = proc.time() - time1;
    writeLines(paste("Downloaded - ", SpecificName, " - Sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
  }

  writeLines("Save data to local disk.");
  ID = str_locate_all(traverseResultFile, "_")[[1]];
  ID = ID[dim(ID)[1], 2];
  for (FileID in 1:length(CopyNumberData))
  {
    filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__", names(CopyNumberData)[FileID], "__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
    write.table(CopyNumberData[[FileID]], file = filename, quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 
    
    # For returning downloaded data
    CopyNumberData[[FileID]] = rbind(colnames(CopyNumberData[[FileID]]), CopyNumberData[[FileID]]);
    colnames(CopyNumberData[[FileID]]) = NULL;
    rownames(CopyNumberData[[FileID]]) = NULL;    
  }  
  
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  # Return downloaded data
  names(CopyNumberData) = paste(SpecificName, names(CopyNumberData), sep = "__");
  return(CopyNumberData);     
}



DownloadMethylationData <- function(traverseResultFile, saveFolderName, cancerType, assayPlatform, tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download DNA methylation data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  load(traverseResultFile);  
  dir.create(path = saveFolderName, recursive = TRUE);
  SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/jhu-usc\\.edu/", assayPlatform, "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
  MeLevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];

  # search for the Sample and Data Relationship Format (SDRF) file of the specified platform and cancer type
  ind = SpecificID[grepEnd(pattern = toupper("\\.sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Program existed due to missing SDRF file."); 
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    return();
  }
  sdrf = toupper(downloadResult$data);
  
  # Process SDRF file, identify the columns of level 3 data file name and TCGA sample barcodes.
  level_3_filename_column = max(grep(pattern = "Data Matrix File", x = sdrf[1, ], ignore.case = TRUE));
  DataLevelColID = max(grep(pattern = "TCGA Data Level", x = sdrf[1, ], ignore.case = TRUE));
  TCGABarcodeID = min(grep(pattern = "TCGA Barcode", x = sdrf[1, ], ignore.case = TRUE));
  colnames(sdrf) = sdrf[1, ];
  sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
  sdrf = sdrf[!duplicated(sdrf[, level_3_filename_column]), c(TCGABarcodeID, DataLevelColID, level_3_filename_column), drop = FALSE];
  SDRFID = sort(union(which(sdrf[, 2] == "LEVEL_3"), which(sdrf[, 2] == "LEVEL 3")), decreasing = FALSE);
  if (length(SDRFID) == 0)
  {
    writeLines("Error: there are no Level 3 data");
    return();
  }
  sdrf = sdrf[SDRFID, , drop = FALSE];  
  
  # If specific patient TCGA barcodes are inputted, only download the specified samples.
  if (!is.null(inputPatientIDs))
  {
    indInputPatientID = c();
    for (i in 1:length(inputPatientIDs))
    {
      indInputPatientID = c(indInputPatientID, grepBeginning(pattern = inputPatientIDs[i], x = sdrf[, 1], ignore.case = FALSE));
    }
    if (length(indInputPatientID) == 0)
    {
      writeLines("No Level 3 data for the inputted TCGA barcodes.");
      return();      
    }else{
      sdrf = sdrf[indInputPatientID, , drop = FALSE];
    }
  }
  
  # Download data of specified tissue
  if (!is.null(tissueType))
  {
    SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                       Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
    sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
  }
  
  if (dim(sdrf)[1] == 0)
  {
    writeLines("No available data.");
    return();
  }
  
  # Download data files of all samples.
  left_columns = NULL;
  AllPosition = NULL;
  exp_names = NULL;
  data = NULL;
  for (i in 1:dim(sdrf)[1])
  {
    time1 = proc.time();
    sample_TCGA_id = sdrf[i, 1];
    ind = MeLevel3ID[grepEnd(pattern = toupper(sdrf[i, 3]), x = upper_file_url[MeLevel3ID], ignore.case = FALSE)];
    if (length(ind) == 0)
    {
      next;
    } 
    if (length(ind) > 1)
    {
      URL = GetNewestURL(AllURL = file_url[ind]);
    }else{
      URL = file_url[ind];
    }
    downloadResult = urlReadTable(url = URL);
    if (downloadResult$errorFlag != 0)
    {
      next;
    }
    s = downloadResult$data;
    s = s[2:dim(s)[1], , drop = FALSE];
    
    chr = rep(0, dim(s)[1]);
    for (j in 1:22)
    {
      IDj = which(s[, 4] == as.character(j));
      chr[IDj] = j;
    }
    IDj = which(toupper(s[, 4]) == "X");
    chr[IDj] = 23;      
    IDj = which(toupper(s[, 4]) == "Y");
    chr[IDj] = 24; 
    position = rep(0, dim(s)[1]);
    position[2:length(position)] = as.numeric(s[2:dim(s)[1], 5]);
    Yj = chr*(10e+10) + position;
    orderIDj = order(Yj, decreasing = FALSE);
    Yj = Yj[orderIDj];
    s = s[orderIDj, , drop = FALSE];

    # Need to check whether every data file has the same methylation probes
    if (is.null(left_columns))
    {
      left_columns = s[, c(1, 3, 4, 5), drop = FALSE];
      AllPosition = Yj;
      exp_names = sample_TCGA_id;
      data = s[2:dim(s)[1], 2, drop = FALSE];
    }else{
      if (sum(AllPosition != Yj) > 0)
      {
        next;
      }
      exp_names = c(exp_names, sample_TCGA_id);
      data = cbind(data, s[2:dim(s)[1], 2, drop = FALSE]);
    }
    
    time = proc.time() - time1;
    writeLines(paste("Downloaded - ", cancerType, "__jhu-usc.edu__", assayPlatform, " - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
  }

  writeLines("Save data to local disk.");
  ID = str_locate_all(traverseResultFile, "_")[[1]];
  ID = ID[dim(ID)[1], 2];
  filename = paste(saveFolderName, "/", outputFileName, cancerType, "__jhu-usc.edu__", assayPlatform, "__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");  
  header = c(left_columns[1, ], exp_names);
  left_columns = left_columns[2:dim(left_columns)[1], , drop = FALSE];
  ID = grepBeginning(pattern = "NA", x = left_columns[, 3], ignore.case = TRUE)  
  ID = sort(setdiff(1:dim(left_columns)[1], ID), decreasing = FALSE);
  left_columns = left_columns[ID, , drop = FALSE];
  data = data[ID, , drop = FALSE];
  write.table(rbind(header, cbind(left_columns, data)), file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  # Return downloaded data
  downloadedData = rbind(header, cbind(left_columns, data));
  rownames(downloadedData) = NULL;
  colnames(downloadedData) = NULL;
  downloadedData;
}



DownloadRNASeqData <-function (traverseResultFile, saveFolderName, cancerType, assayPlatform, dataType = "", tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
{
  
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download Gene Expression data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }  
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  load(traverseResultFile);  
  dir.create(path = saveFolderName, recursive = TRUE);
  if (assayPlatform == "RNASeqV1")
  {
    platform = c("illuminaga_rnaseq", "illuminahiseq_rnaseq"); 
    Institution = c("unc.edu", "bcgsc.ca");    
  }
  if (assayPlatform == "RNASeqV2")
  {
    platform = c("illuminaga_rnaseqv2", "illuminahiseq_rnaseqv2"); 
    Institution = c("unc.edu");    
  }
  if (assayPlatform == "Microarray")
  {
    platform = c("agilentg4502a_07_3", "ht_hg-u133a", "agilentg4502a_07_1", "agilentg4502a_07_2", "hg-u133_plus_2");
    Institution = c("unc.edu", "broad.mit.edu", "genome.wustl.edu");
    dataType = "";
  }
  
  # For returning downloaded data
  downloadedData = vector("list", 0);     
  dataIndex = 0;     
  
  # download RNASeqV2 data
  if (assayPlatform == "RNASeqV2")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # Only keep the file information for the data types that should be downloaded.
        keepID = c();
        for (keep_i in 1:length(dataType))
        {
          keepID = c(keepID, grep(pattern = dataType[keep_i], x = sdrf[, 2], ignore.case = TRUE))
        }
        sdrf = sdrf[sort(unique(keepID), decreasing = FALSE), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
        
        # Download data of specified tissue
        if (!is.null(tissueType))
        {
          SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                             Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
          sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
        }
        
        if (dim(sdrf)[1] == 0)
        {
          next;
        }
        
        # Start to download data files
        exp_gene = NULL;
        column_gene = NULL;
        gene_left_column = NULL;
        gene_data = NULL;
        exp_gene_normalized = NULL;
        column_gene_normalized = NULL;
        gene_normalized_left_column = NULL;
        gene_normalized_data = NULL;    
        
        exp_isoform = NULL;
        column_isoform = NULL;
        isoform_left_column = NULL;
        isoform_data = NULL;
        exp_isoform_normalized = NULL;
        column_isoform_normalized = NULL;
        isoform_normalized_left_column = NULL;
        isoform_normalized_data = NULL;    
        
        exp_exon = NULL;
        column_exon = NULL;
        exon_left_column = NULL;
        exon_data=NULL;
        exp_junction = NULL;
        column_junction = NULL;
        junction_left_column = NULL;
        junction_data=NULL;
        
        for (i in 1:dim(sdrf)[1])
        {
          time1 = proc.time();
          
          sample_TCGA_id = sdrf[i, 1];
          ind = RNALevel3ID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[RNALevel3ID], ignore.case = FALSE)];
          if (length(ind) == 0)
          {
            next;
          }
          if (length(ind) > 1)
          {
            URL = GetNewestURL(AllURL = file_url[ind]);
          }else{
            URL = file_url[ind];
          }  
          downloadResult = urlReadTable(url = URL);
          if (downloadResult$errorFlag != 0)
          {
            next;
          }
          s = downloadResult$data;
          
          # read gene expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.genes.results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(gene_left_column))
              {
                gene_left_column = s[, c(1, 4), drop = FALSE];
                gene_data = s[, c(2, 3), drop = FALSE];
                exp_gene = c(sample_TCGA_id, sample_TCGA_id);
                column_gene = c("raw_count", "scaled_estimate");
              }else{
                if (sum(gene_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                gene_data = cbind(gene_data, s[, c(2, 3), drop = FALSE]);
                exp_gene = c(exp_gene, sample_TCGA_id, sample_TCGA_id);
                column_gene = c(column_gene, "raw_count", "scaled_estimate");
              }
            }
          }
          
          # read gene normalized data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.genes.normalized_results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(gene_normalized_left_column))
              {
                gene_normalized_left_column = s[, 1, drop = FALSE];
                gene_normalized_data = s[, 2, drop = FALSE];
                exp_gene_normalized = sample_TCGA_id;
                column_gene_normalized = "normalized_count";
              }else{
                if (sum(gene_normalized_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                gene_normalized_data = cbind(gene_normalized_data, s[, 2, drop = FALSE]);
                exp_gene_normalized = c(exp_gene_normalized, sample_TCGA_id);
                column_gene_normalized = c(column_gene_normalized, "normalized_count");
              }
            }
          }          
          
          # read isoform data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.isoforms.results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(isoform_left_column))
              {
                isoform_left_column = s[, 1, drop = FALSE];
                isoform_data = s[, c(2, 3), drop = FALSE];
                exp_isoform = c(sample_TCGA_id, sample_TCGA_id);
                column_isoform = c("raw_count", "scaled_estimate");
              }else{
                if (sum(isoform_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                isoform_data = cbind(isoform_data, s[, c(2, 3), drop = FALSE]);
                exp_isoform = c(exp_isoform, sample_TCGA_id, sample_TCGA_id);
                column_isoform = c(column_isoform, "raw_count", "scaled_estimate");
              }
            }
          }
          
          # read isoform normalized data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.isoforms.normalized_results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(isoform_normalized_left_column))
              {
                isoform_normalized_left_column = s[, 1, drop = FALSE];
                isoform_normalized_data = s[, 2, drop = FALSE];
                exp_isoform_normalized = sample_TCGA_id;
                column_isoform_normalized = "normalized_count";
              }else{
                if (sum(isoform_normalized_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                isoform_normalized_data = cbind(isoform_normalized_data, s[, 2, drop = FALSE]);
                exp_isoform_normalized = c(exp_isoform_normalized, sample_TCGA_id);
                column_isoform_normalized = c(column_isoform_normalized, "normalized_count");
              }
            }
          }          
          
          # read exon data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("exon_quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(exon_left_column))
              {
                exon_left_column = s[, 1, drop = FALSE];
                exon_data = s[, 2:4, drop = FALSE];
                exp_exon = c(sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c("raw_counts", "median_length_normalized", "RPKM");
              }else{
                if (sum(exon_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                exon_data = cbind(exon_data, s[, 2:4, drop = FALSE]);
                exp_exon = c(exp_exon, sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c(column_exon, "raw_counts", "median_length_normalized", "RPKM");
              }
            }
          }
          
          # read junction data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("junction_quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(junction_left_column))
              {
                junction_left_column = s[, 1, drop = FALSE];
                junction_data = s[, 2, drop = FALSE];
                exp_junction = sample_TCGA_id;
                column_junction = "raw_counts";
              }else{
                if (sum(junction_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                junction_data = cbind(junction_data, s[, 2, drop = FALSE]);
                exp_junction = c(exp_junction, sample_TCGA_id);
                column_junction = c(column_junction, "raw_counts");
              }
            }
          }          
          
          time = proc.time() - time1;
          writeLines(paste("Downloaded - ", SpecificName, " - file ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        if (length(grep(pattern = "rsem.genes.results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.genes.results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", "Hybridization REF", exp_gene), c("gene_id", "transcript_id", column_gene), cbind(gene_left_column, gene_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", "Hybridization REF", exp_gene), c("gene_id", "transcript_id", column_gene), cbind(gene_left_column, gene_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.genes.results", sep = "__");
        }
        if (length(grep(pattern = "rsem.genes.normalized_results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.genes.normalized_results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_gene_normalized), c("gene_id", column_gene_normalized), cbind(gene_normalized_left_column, gene_normalized_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_gene_normalized), c("gene_id", column_gene_normalized), cbind(gene_normalized_left_column, gene_normalized_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.genes.normalized_results", sep = "__");
        }
        if (length(grep(pattern = "rsem.isoforms.results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.isoforms.results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_isoform), c("isoform_id", column_isoform), cbind(isoform_left_column, isoform_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_isoform), c("isoform_id", column_isoform), cbind(isoform_left_column, isoform_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.isoforms.results", sep = "__");
        }
        if (length(grep(pattern = "rsem.isoforms.normalized_results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.isoforms.normalized_results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_isoform_normalized), c("isoform_id", column_isoform_normalized), cbind(isoform_normalized_left_column, isoform_normalized_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_isoform_normalized), c("isoform_id", column_isoform_normalized), cbind(isoform_normalized_left_column, isoform_normalized_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.isoforms.normalized_results", sep = "__");
        }
        if (length(grep(pattern = "exon_quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__exon_quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_exon), c("exon", column_exon), cbind(exon_left_column, exon_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_exon), c("exon", column_exon), cbind(exon_left_column, exon_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "exon_quantification", sep = "__");
        }
        if (length(grep(pattern = "junction_quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__junction_quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_junction), c("junction", column_junction), cbind(junction_left_column, junction_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_junction), c("junction", column_junction), cbind(junction_left_column, junction_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "junction_quantification", sep = "__");
        }
      }
    }
  }
  
  # download RNASeqV1 data
  if (assayPlatform == "RNASeqV1")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # Only keep the file information for the data types that should be downloaded.
        keepID = c();
        for (keep_i in 1:length(dataType))
        {
          keepID = c(keepID, grep(pattern = dataType[keep_i], x = sdrf[, 2], ignore.case = TRUE))
        }
        sdrf = sdrf[sort(unique(keepID), decreasing = FALSE), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
        
        # Download data of specified tissue
        if (!is.null(tissueType))
        {
          SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                             Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
          sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
        }
        
        if (dim(sdrf)[1] == 0)
        {
          next;
        }
        
        exp_names_gene = NULL;
        column_gene = NULL;
        gene_left_column = NULL;
        gene_RPKM = NULL;
        
        exp_names_exon = NULL;
        column_exon = NULL;
        exon_left_column = NULL;
        exon_RPKM = NULL;
        
        junction_count = NULL;
        exp_names_junction = NULL;    
        column_junction = NULL;     
        junction_left_column = NULL;        
        
        for (i in 1:dim(sdrf)[1])
        {
          time1 = proc.time();
          
          sample_TCGA_id = sdrf[i, 1];
          ind = RNALevel3ID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[RNALevel3ID], ignore.case = FALSE)];
          if (length(ind) == 0)
          {
            next;
          }
          if (length(ind) > 1)
          {
            URL = GetNewestURL(AllURL = file_url[ind]);
          }else{
            URL = file_url[ind];
          }  
          downloadResult = urlReadTable(url = URL);
          if (downloadResult$errorFlag != 0)
          {
            next;
          }
          s = downloadResult$data;
          
          # read gene expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("gene.quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(gene_left_column))
              {
                gene_left_column = s[, 1, drop = FALSE];
                gene_RPKM = s[, 2:4, drop = FALSE];
                exp_names_gene = c(sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_gene = c("raw_counts", "median_length_normalized", "RPKM");
              }else{
                if (sum(gene_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                gene_RPKM = cbind(gene_RPKM, s[, 2:4, drop = FALSE]);
                exp_names_gene = c(exp_names_gene, sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_gene = c(column_gene, "raw_counts", "median_length_normalized", "RPKM");
              }
            }
          }
          
          # read exon expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("exon.quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(exon_left_column))
              {
                exon_left_column = s[, 1, drop = FALSE];
                exon_RPKM = s[, 2:4, drop = FALSE];
                exp_names_exon = c(sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c("raw_counts", "median_length_normalized", "RPKM");
              }else{
                if (sum(exon_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                exon_RPKM = cbind(exon_RPKM, s[, 2:4, drop = FALSE]);
                exp_names_exon = c(exp_names_exon, sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c(column_exon, "raw_counts", "median_length_normalized", "RPKM");
              }
            }
          }
          
          # read junction expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("spljxn.quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(junction_left_column))
              {
                junction_left_column = s[, 1, drop = FALSE];
                junction_count = s[, 2, drop = FALSE];
                exp_names_junction = sample_TCGA_id;
                column_junction = "raw_counts";
              }else{
                if (sum(junction_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                junction_count = cbind(junction_count, s[, 2, drop = FALSE]);
                exp_names_junction = c(exp_names_junction, sample_TCGA_id);
                column_junction = c(column_junction, "raw_counts");
              }
            }
          }
          
          time = proc.time() - time1;
          writeLines(paste("Downloaded - ", SpecificName, " - file ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        if (length(grep(pattern = "gene.quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__gene.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names_gene), c("gene", column_gene), cbind(gene_left_column, gene_RPKM)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1;          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_gene), c("gene", column_gene), cbind(gene_left_column, gene_RPKM));
          names(downloadedData)[dataIndex] = paste(SpecificName, "gene.quantification", sep = "__");
        }
        if (length(grep(pattern = "exon.quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__exon.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names_exon), c("exon", column_exon), cbind(exon_left_column, exon_RPKM)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1;
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_exon), c("exon", column_exon), cbind(exon_left_column, exon_RPKM));
          names(downloadedData)[dataIndex] = paste(SpecificName, "exon.quantification", sep = "__");          
        }
        if (length(grep(pattern = "spljxn.quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__spljxn.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names_junction), c("junction", column_junction), cbind(junction_left_column, junction_count)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1;
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_junction), c("junction", column_junction), cbind(junction_left_column, junction_count));
          names(downloadedData)[dataIndex] = paste(SpecificName, "spljxn.quantification", sep = "__");             
        }
      }
    }
  }
  
  # download Microarray data
  if (assayPlatform == "Microarray")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        # Derived Array Data Matrix File
        level_3_filename_column = max(grep(pattern = "Derived Array Data Matrix File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
        
        # Download data of specified tissue
        if (!is.null(tissueType))
        {
          SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                             Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
          sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
        }
        
        if (dim(sdrf)[1] == 0)
        {
          next;
        }
        
        exp_names_gene = NULL;
        column_gene = NULL;
        gene_left_column = NULL;
        gene_RPKM = NULL;
        for (i in 1:dim(sdrf)[1])
        {
          time1 = proc.time();
          
          sample_TCGA_id = sdrf[i, 1];
          ind = RNALevel3ID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[RNALevel3ID], ignore.case = FALSE)];
          if (length(ind) == 0)
          {
            next;
          }
          if (length(ind) > 1)
          {
            URL = GetNewestURL(AllURL = file_url[ind]);
          }else{
            URL = file_url[ind];
          }  
          downloadResult = urlReadTable(url = URL);
          if (downloadResult$errorFlag != 0)
          {
            next;
          }
          s = downloadResult$data;
          
          column_gene = c(column_gene, s[2, 2]);
          s = s[3:dim(s)[1], , drop = FALSE];
          I_order_probes = order(s[, 1], decreasing = FALSE);
          s = s[I_order_probes, , drop = FALSE];
          if (is.null(gene_left_column))
          {
            gene_left_column = s[, 1, drop = FALSE];
            gene_RPKM = s[, 2, drop = FALSE];
            exp_names_gene = c(sample_TCGA_id);
          }else{
            if (sum(gene_left_column[, 1] != s[, 1]) > 0)
            {
              next;
            }
            gene_RPKM = cbind(gene_RPKM, s[, 2, drop = FALSE]);
            exp_names_gene = c(exp_names_gene, sample_TCGA_id);
          }
          
          time = proc.time() - time1;
          writeLines(paste("Downloaded - ", SpecificName, " - file ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__gene.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
        write.table(rbind(c("Hybridization REF", exp_names_gene), c("Composite Element REF", column_gene), cbind(gene_left_column, gene_RPKM)), 
                    file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
        
        # For returning downloaded data
        dataIndex = dataIndex + 1;          
        downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_gene), c("Composite Element REF", column_gene), cbind(gene_left_column, gene_RPKM));
        names(downloadedData)[dataIndex] = paste(SpecificName, "gene.quantification", sep = "__");
        
      }
    }
  }
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  # Return downloaded data
  if (length(downloadedData) > 0)
  {
    for (i in 1:length(downloadedData))
    {
      rownames(downloadedData[[i]]) = NULL;
      colnames(downloadedData[[i]]) = NULL;
    }
  }
  downloadedData;
}



DownloadRPPAData <- function(traverseResultFile, saveFolderName, cancerType, assayPlatform = "mda_rppa_core", tissueType = NULL, inputPatientIDs = NULL, outputFileName = "") 
{
  options(warn=-1);

  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download RPPA protein expression data of ", cancerType, " patients.", sep = ""));

  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }  
    
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");  
  load(traverseResultFile);
  dir.create(path = saveFolderName, recursive = TRUE);
  DirURL_ID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/mdanderson\\.org/", assayPlatform, "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
  IDLevel_3InFile_Url = DirURL_ID[grep(pattern = toupper("Level_3"), x = upper_file_url[DirURL_ID], ignore.case = FALSE)];

  # download protein antibody annotation file
  FileEndURL = paste(cancerType, "\\.MDA_RPPA_Core\\.antibody_annotation\\.txt", sep = "");
  ind = DirURL_ID[grepEnd(pattern = toupper(FileEndURL), x = upper_file_url[DirURL_ID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Error: no antibody annotation file.");
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }  
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading antibody annotation file.");
    return();
  }
  Annotation = downloadResult$data;
  colnames(Annotation) = Annotation[1, ];
  Annotation = Annotation[2:dim(Annotation)[1], , drop = FALSE];
  GeneNameColID = grepBeginning(pattern = toupper("Gene Name"), x = colnames(Annotation), ignore.case = TRUE);
  REFColID = grepBeginning(pattern = toupper("Composite Element REF"), x = colnames(Annotation), ignore.case = TRUE);
  Annotation = unique(Annotation[, c(GeneNameColID, REFColID), drop = FALSE]);
  replaceTable = cbind(wrongSymbol = c("ACACAACACB", "AKT1AKT2 AKT3", "GSK3AGSK3B", "MAPK1MAPK3", "RAB11ARAB11B", "BIRC2 "), 
                      correctSymbol = c("ACACA ACACB", "AKT1 AKT2 AKT3", "GSK3A GSK3B", "MAPK1 MAPK3", "RAB11A RAB11B", "BIRC2"));
  for (iReplace in 1:dim(replaceTable)[1])
  {
    jReplace = which(Annotation[, 1] == replaceTable[iReplace, 1]);
    Annotation[jReplace, 1] = replaceTable[iReplace, 2];
  }
  idTmp = grepBeginning(pattern = "CDK1-", x = Annotation[, 2], ignore.case = TRUE);
  Annotation[idTmp, 1] = "CDK1";
  
  # download Sample and Data Relationship Format (SDRF) file
  writeLines("Download and process Sample and Data Relationship Format (SDRF) file.");
  FileEndURL = paste(cancerType, "\\.MDA_RPPA_Core\\.sdrf\\.txt", sep = "");
  ind = DirURL_ID[grepEnd(pattern = toupper(FileEndURL), x = upper_file_url[DirURL_ID], ignore.case = FALSE)];
  if (length(ind) == 0)
  {
    writeLines("Error: no SDRF file.");
    return();
  }
  if (length(ind) > 1)
  {
    URL = GetNewestURL(AllURL = file_url[ind]);
  }else{
    URL = file_url[ind];
  }  
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    return();
  }
  sdrf = toupper(downloadResult$data);  
  
  level_3_filename_column = max(grep(pattern = "Derived Array Data Matrix File", x = sdrf[1, ], ignore.case = TRUE));
  LevelInfoColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
  SampleNameColID = min(grep(pattern = "Sample Name", x = sdrf[1, ], ignore.case = TRUE));
  ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));
  colnames(sdrf) = sdrf[1, ];
  sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
  sdrf = sdrf[!duplicated(sdrf[, level_3_filename_column]), , drop = FALSE];
  SDRFID = sort(union(which(sdrf[, LevelInfoColID] == "LEVEL_3"), which(sdrf[, LevelInfoColID] == "LEVEL 3")), decreasing = FALSE);
  if (length(SDRFID) == 0)
  {
    writeLines("Error: there are no Level 3 data");
    return();
  }
  sdrf = sdrf[SDRFID, c(ExtractNameColID, SampleNameColID, level_3_filename_column), drop = FALSE];  
  
  # If specific patient TCGA barcodes are inputted, only download data of the specified samples.
  if (!is.null(inputPatientIDs))
  {
    indInputPatientID = c();
    for (i in 1:length(inputPatientIDs))
    {
      indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
    }
    if (length(indInputPatientID) == 0)
    {
      writeLines("No Level 3 data for the inputted TCGA barcodes.");
      return();      
    }else{
      sdrf = sdrf[indInputPatientID, , drop = FALSE];
    }
  }
  
  # Download data of specified tissue
  if (!is.null(tissueType))
  {
    SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                       Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
    sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
  }
  
  if (dim(sdrf)[1] == 0)
  {
    writeLines("No available data.");
    return();
  }
  
  left_columns = NULL;
  exp_names = NULL;
  data = NULL;
  for (i in 1:dim(sdrf)[1])
  {
    time1 = proc.time();
    sample_TCGA_id = strsplit(sdrf[i, 1], split = "\\.")[[1]][1];
    ind = IDLevel_3InFile_Url[grepEnd(pattern = toupper(sdrf[i, 3]), x = upper_file_url[IDLevel_3InFile_Url], ignore.case = FALSE)];
    if (length(ind) == 0)
    {
      next;
    }
    if (length(ind) > 1)
    {
      URL = GetNewestURL(AllURL = file_url[ind]);
    }else{
      URL = file_url[ind];
    }  
    downloadResult = urlReadTable(url = URL);
    if (downloadResult$errorFlag != 0)
    {
      next;
    }
    s = downloadResult$data[3:dim(downloadResult$data)[1], , drop = FALSE];
    I_order_probes = order(s[, 1], decreasing = FALSE);
    s = s[I_order_probes, , drop = FALSE];
    if (is.null(left_columns))
    {
      left_columns = s[, 1];
      exp_names = sample_TCGA_id;
      data = s[, 2, drop = FALSE];
    }else{
      if (sum(left_columns != s[ ,1]) > 0)
      {
        next;
      }
      exp_names = c(exp_names, sample_TCGA_id);
      data = cbind(data, s[, 2, drop = FALSE]);
    }
    time = proc.time() - time1;
    writeLines(paste("Downloaded - ", cancerType, "__mdanderson.org__", assayPlatform, " - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
  }
  colnames(data) = exp_names;

  IX = which(toupper(Annotation[, 2]) %in% toupper(left_columns));
  if(length(IX) < length(left_columns))
  {
    writeLines("Downloaded data have antibodies not included in the antibody annotation file.");  
  }
  left_columns_temp = c();
  for (i in 1:length(left_columns))
  {
    IDi = which(toupper(Annotation[, 2]) == toupper(left_columns)[i]);
    if (length(IDi) == 0)
    {
      left_columns_temp = c(left_columns_temp, paste("GeneSymbolNotFound|", left_columns[i], sep = ""));      
    } else {
      IDi = IDi[1];
      left_columns_temp = c(left_columns_temp, paste(Annotation[IDi, 1], "|", Annotation[IDi, 2], sep = ""));
    }
  }
  left_columns = left_columns_temp;   

  writeLines("Save data to local disk.");
  ID = str_locate_all(traverseResultFile, "_")[[1]];
  ID = ID[dim(ID)[1], 2];
  filename = paste(saveFolderName, "/", outputFileName, cancerType, "__mdanderson.org__", assayPlatform, "__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
  write.table(cbind(Composite.Element.REF = left_columns, data), file = filename, quote = FALSE, sep = "\t", na = "", col.names = TRUE, row.names = FALSE);  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  # Return downloaded data
  downloadedData = cbind(Composite.Element.REF = left_columns, data);
  downloadedData = rbind(colnames(downloadedData), downloadedData);
  rownames(downloadedData) = NULL;
  colnames(downloadedData) = NULL;
  downloadedData;
}



DownloadmiRNASeqData <- function(traverseResultFile, saveFolderName, cancerType, assayPlatform = "miRNASeq", tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download miRNA-seq data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }  
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }
  
  # For returning downloaded data
  downloadedData = vector("list", 0);     
  dataIndex = 0;     
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  load(traverseResultFile);
  if (assayPlatform == "miRNASeq")
  {
    platform = c("illuminaga_mirnaseq", "illuminahiseq_mirnaseq");
  }
  dir.create(path = saveFolderName, recursive = TRUE);
  miRNALevel3ID = grep(pattern = toupper("miRNASeq\\.Level_3"), x = upper_file_url, ignore.case = FALSE);
  
  for (IDpl in 1:length(platform))
  {
    SpecificName = paste(cancerType, "__", "bcgsc.ca", "__", platform[IDpl], sep = "");
    SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/bcgsc\\.ca/", platform[IDpl], sep = "")), x = upper_file_url, ignore.case = FALSE);
    ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
    if (length(ind) == 0)
    {
      next;
    }
    if (length(ind) > 1)
    {
      URL = GetNewestURL(AllURL = file_url[ind]);
    }else{
      URL = file_url[ind];
    }  
    downloadResult = urlReadTable(url = URL);
    if (downloadResult$errorFlag != 0)
    {
      next;
    }
    sdrf = toupper(downloadResult$data);  
  
    # Process SDRF file, identify the columns of level 3 data file name and TCGA sample barcodes.
    level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
    DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
    ExtractNameColID = min(grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE));
    RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
    if (length(ExtractNameColID) == 0)
    {
      ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));
    }
    colnames(sdrf) = sdrf[1, ];
    sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
    sdrf = sdrf[!duplicated(sdrf[, level_3_filename_column]), , drop = FALSE];
    
    Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
    Level3_ID = sort(intersect(Level3_ID, grep(pattern = toupper("mirna\\.quantification"), 
                x = sdrf[, level_3_filename_column], ignore.case = FALSE)), decreasing = FALSE);
    if (length(Level3_ID) == 0)
    {
      next;
    }
    sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];    

    # If specific patient TCGA barcodes are inputted, only download the specified samples.
    if (!is.null(inputPatientIDs))
    {
      indInputPatientID = c();
      for (i in 1:length(inputPatientIDs))
      {
        indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
      }
      if (length(indInputPatientID) == 0)
      {
        next;      
      }else{
        sdrf = sdrf[indInputPatientID, , drop = FALSE];
      }
    }
    
    # Download data of specified tissue
    if (!is.null(tissueType))
    {
      SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                         Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
      sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
    }
    
    if (dim(sdrf)[1] == 0)
    {
      next;
    }
    
    sdrfO = sdrf;
    for (RefGId in 1:2)
    {
      if (RefGId == 1)
      {
        sdrf = sdrfO[grep(pattern = toupper("NCBI36"), x = sdrfO[, 3], ignore.case = FALSE), , drop = FALSE];  
      }
      if (RefGId == 2)
      {
        sdrf = sdrfO[grep(pattern = toupper("GRCh37"), x = sdrfO[, 3], ignore.case = FALSE), , drop = FALSE];  
      }
      if (dim(sdrf)[1] == 0)
      {
        next;
      }
      
      exp_names = NULL;
      gene_left_column = NULL;
      gene_RPM = NULL;
      column_gene = NULL;
      for (i in 1:dim(sdrf)[1])
      {
        time1 = proc.time();
        
        sample_TCGA_id = sdrf[i, 1];
        ind = intersect(miRNALevel3ID, SpecificID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[SpecificID], ignore.case = FALSE)]);
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }  
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        s = downloadResult$data[2:dim(downloadResult$data)[1], , drop = FALSE];
        I_order_probes = order(s[, 1], decreasing = FALSE);
        s = s[I_order_probes, , drop = FALSE];
        if (is.null(gene_left_column))
        {
          gene_left_column = s[, 1, drop = FALSE];
          exp_names = c(sample_TCGA_id, sample_TCGA_id);
          gene_RPM = s[, 2:3, drop = FALSE];
          column_gene = c("read_count", "reads_per_million_miRNA_mapped");
        }else{
          if (sum(gene_left_column[, 1] != s[ ,1]) > 0)
          {
            next;
          }
          exp_names = c(exp_names, sample_TCGA_id, sample_TCGA_id);
          gene_RPM = cbind(gene_RPM, s[, 2:3, drop = FALSE]);
          column_gene = c(column_gene, "read_count", "reads_per_million_miRNA_mapped");
        }
        
        time = proc.time() - time1;
        if (RefGId == 1)
        {
          writeLines(paste("Downloaded - ", SpecificName, " - NCBI36 - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        if (RefGId == 2)
        {
          writeLines(paste("Downloaded - ", SpecificName, " - GRCh37 - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
      }

      if (!is.null(gene_RPM))
      {
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        if (RefGId == 1)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__NCBI36__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM)), file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex+1;     
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM));
          rownames(downloadedData[[dataIndex]]) = NULL;  
          colnames(downloadedData[[dataIndex]]) = NULL;          
          names(downloadedData)[dataIndex] = paste(SpecificName, "__NCBI36", sep = "");
        }
        if (RefGId == 2)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__GRCh37__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM)), file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  

          # For returning downloaded data
          dataIndex = dataIndex+1;     
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM));
          rownames(downloadedData[[dataIndex]]) = NULL;  
          colnames(downloadedData[[dataIndex]]) = NULL;
          names(downloadedData)[dataIndex] = paste(SpecificName, "__GRCh37", sep = "");
        }
      }
    
    }
  }
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  downloadedData;
}



######################### Auxiliary Functions of Module A #####################################################

# urlReadTable is the function to read a data table from a website and pass it to a variable.

# Input arguments:
# url: URl of the website from which data table will be obtained.

# Output argument:
# data: a character matrix holding the data table obtained from website. 

urlReadTable <- function(url)
{
  
  data = try(content(GET(url), as = "text"), silent = TRUE);
  if (class(data) == "try-error")
  {
    return(list(data = data, errorFlag = 1));    
  }
  if (length(grep(pattern = "HTTP 404: Page Not Found\n\nThe page you requested was not found.", x = data, ignore.case = TRUE)) > 0)
  {
    return(list(data = data, errorFlag = 2)); 
  }
  data = gsub(pattern = "\r", replacement = "", x = data);
  data = as.matrix(read.table(text = data, sep="\t", fill = TRUE, quote = NULL, check.names = FALSE));
  
#   data = strsplit(data, split = "\n")[[1]];
#   numCol = length(strsplit(data[1], split = "\t")[[1]]);
#   data = unlist(strsplit(data, split = "\t"));
#   if ((length(data) %% numCol) != 0)
#   {
#     return(list(data = data, errorFlag = 3)); 
#   }
#   data = t(matrix(data, numCol, length(data)/numCol));
  
  return(list(data = data, errorFlag = 0));
}



# downloadFile is a function to download content from a website and save it as a local file.

# Input arguments:
# url: URl of the website whose content to be obtained.
# saveFileName: path and name of the file to store the web content.

downloadFile <- function(url, saveFileName)
{
  data = try(content(GET(url), as = "text"), silent = TRUE);
  if (class(data) == "try-error")
  {
    return(errorFlag = 1);    
  }
  if (length(grep(pattern = "HTTP 404: Page Not Found\n\nThe page you requested was not found.", x = data, ignore.case = TRUE)) > 0)
  {
    return(errorFlag = 2); 
  }
  write(x = data, file = saveFileName);
  return(errorFlag = 0);
}



# grepEnd is a function similar to grep but identifies the strings with pattern at the end of the strings.
grepEnd <- function(pattern, x, ignore.case = FALSE)
{
  ind = grep(pattern = pattern, x = x, ignore.case = ignore.case);
  if (ignore.case)
  {
    ind = ind[nchar(x[ind]) == sapply(str_locate_all(toupper(x[ind]), pattern = toupper(pattern)), function(y){y[dim(y)[1], 2]})];    
  }else{
    ind = ind[nchar(x[ind]) == sapply(str_locate_all(x[ind], pattern = pattern), function(y){y[dim(y)[1], 2]})];    
  }
}



# grepBeginning is a function similar to grep but identifies the strings with pattern at the Beginning of the strings.
grepBeginning <- function(pattern, x, ignore.case = FALSE)
{
  ind = grep(pattern = pattern, x = x, ignore.case = ignore.case);
  if (ignore.case)
  {
    ind = ind[rep(1, length(ind)) == sapply(str_locate_all(toupper(x[ind]), pattern = toupper(pattern)), function(y)y[1, 1])];    
  }else{
    ind = ind[rep(1, length(ind)) == sapply(str_locate_all(x[ind], pattern = pattern), function(y)y[1, 1])];    
  }
}



# This function analyzes a TCGA webpage and identify all files and directories on it. 
# The URLs of files and directories are returned using two character vectors.

# Input arguments:
# TCGA_link: a string of URL for the TCGA webpage to be analyzed.

# Output arguments:
# file_url: a character vector, including the URLs of all files on the webpage.
# dir_url: a character vector, including the URLs of all directories on the webpage.

getTCGA_URL<-function(TCGA_link)
{
  file_url = c();
  file_name = c();
  is_dir = c();
  
  if (substr(TCGA_link, nchar(TCGA_link), nchar(TCGA_link)) != "/")
  {
    writeLines(paste("Add / to the end of TCGA_link for ", TCGA_link, sep = ""));
    TCGA_link = paste(TCGA_link, "/", sep = "");
  }
  
  s = try(content(GET(TCGA_link), as = "text"), silent = TRUE);
  if (class(s) == "try-error")
  {
    return(list(file_url = c(), dir_url = c()));    
  }
  if (length(grep(pattern = "HTTP 404: Page Not Found\n\nThe page you requested was not found.", x = s, ignore.case = TRUE)) > 0)
  {
    return(list(file_url = c(), dir_url = c())); 
  }
  s = substr(s, str_locate(string = s, pattern = "Parent Directory")[1, 2] + 5, nchar(s));
  
  start_points = str_locate_all(toupper(s), "<A HREF")[[1]][, 1];
  end_points = str_locate_all(toupper(s), "</A>")[[1]][, 1];
  if (length(start_points) != length(end_points))
  {
    stop(paste("Error in parsing URL", TCGA_link, sep = ", "));
  }  
  if (length(end_points) == 0)
  {
    return(list(file_url = file_url, is_dir = is_dir));
  }
  
  for (i in 1:length(start_points))
  {
    tmp = substr(s, start_points[i]+1, end_points[i]);
    ind = str_locate_all(string = tmp, pattern = "\"")[[1]][, 1];
    file_url[i] = substr(tmp, ind[1]+1, ind[2]-1);
    is_dir[i] = (substr(file_url[i], nchar(file_url[i]), nchar(file_url[i])) == "/");
  }
  
  for (i in 1:length(file_url))
  {
    HTTP_i = str_locate_all(string = toupper(file_url[i]), pattern = "HTTP")[[1]];
    if (dim(HTTP_i)[1] == 0)
    {
      file_url[i] = paste(TCGA_link, file_url[i], sep = "");
    }
  }
  
  return(list(file_url = file_url[!is_dir], dir_url = file_url[is_dir]));
}



# This function selects the newest file URL from all the input URLs by considering the
# last folder name tail numbers. The number format is *.*.*

GetNewestURL <- function(AllURL)
{
  SeriesNum = matrix(rep("", length(AllURL)*3), length(AllURL), 3);
  NumLength = rep(0, 3);
  SN = rep("", length(AllURL));
  for (i in 1:length(AllURL))
  {
    Str = AllURL[i];
    SepID = str_locate_all(string = Str, pattern = "/")[[1]];
    SepID = SepID[(dim(SepID)[1]-1):dim(SepID)[1], 1];
    Str = substr(Str, SepID[1]+1, SepID[2]-1);
    Str = strsplit(Str, split = "\\.")[[1]];
    SeriesNum[i, ] = Str[(length(Str)-2) : length(Str)];
    NumLength[1] = max(NumLength[1], nchar(SeriesNum[i, 1]));
    NumLength[2] = max(NumLength[2], nchar(SeriesNum[i, 2]));
    NumLength[3] = max(NumLength[3], nchar(SeriesNum[i, 3]));    
  }
  
  for (i in 1:length(AllURL))
  {
    for (j in 1:3)
    {
      if (nchar(SeriesNum[i, j]) < NumLength[j])
      {
        SeriesNum[i, j] = paste(paste(rep("0", NumLength[j] - nchar(SeriesNum[i, j])), collapse = ""), SeriesNum[i, j], sep = "");
      }
    }
    SN[i] = paste(SeriesNum[i, ], collapse = "");
  }
  AllURL[which.max(as.numeric(SN))];
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
