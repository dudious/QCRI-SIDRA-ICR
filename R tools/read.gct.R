##Function to read the file.gct
read.gct <- function(file) {
  if (is.character(file)) 
    if (file == "") 
      file <- stdin()
  else {
    file <- file(file, "r")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
    stop("argument `file' must be a character string or connection")
  
  # line 1 version
  version <- readLines(file, n=1) 
  
  # line 2 dimensions
  dimensions <- scan(file, what=list("integer", "integer"), nmax=1, quiet=TRUE)   
  rows <- dimensions[[1]]
  columns <- dimensions[[2]]
  # line 3 Name\tDescription\tSample names...
  column.names <- read.table(file, header=FALSE, nrows=1, sep="\t", fill=FALSE) 
  column.names <-column.names[3:length(column.names)]
  
  
  if(length(column.names)!=columns) {
    stop("Number of sample names not equal to the number of columns.")  
  }
  
  colClasses <- c(rep(c("character"), 2), rep(c("double"), columns))
  
  x <- read.table(file, header=FALSE, quote="", row.names=NULL, comment.char="", sep="\t", colClasses=colClasses, fill=FALSE)
  row.descriptions <- as.character(x[,2]) 
  data <- x[seq(from=3, to=dim(x)[2], by=1)]
  
  column.names <- column.names[!is.na(column.names)]
  names(data) <- column.names
  
  row.names(data) <- x[,1]
  return(list(row.descriptions=row.descriptions, data=data))
}
