## Create function that installs and loads the required packages 
## as specified in character vector (required.packages)



ipak <- function(required.packages){
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) {install.packages(missing.packages, dependencies = TRUE)}
  invisible(sapply(required.packages, library, character.only = TRUE))
}

ibiopak <- function(required.bioconductor.packages){
  source("https://bioconductor.org/biocLite.R")
  missing.packages <- required.bioconductor.packages[!(required.bioconductor.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) {biocLite(missing.packages, dependencies = TRUE)}
  invisible(sapply(required.bioconductor.packages, library, character.only = TRUE))
}