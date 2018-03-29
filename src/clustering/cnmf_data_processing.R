# Filename: CNMF_Data_Processing.R
#  Authors: Pablo Tamayo (tamayo@genome.wi.mit.edu), R. Zupko II (rzupko@broadinstitute.org)
#
#  Purpose: Used by GDAC_CNMF.R to load data and preprocess the data sets.
options(warn = -1)

CNMF.read.dataset <- function(file) {
  # Read CNMF datasets (either gct file or res file)
  result <- regexpr(paste(".gct","$",sep=""), tolower(file))
  if(result[[1]] != -1) {
      return(CNMF.read.gct(file))
  }
  result <- regexpr(paste(".res","$",sep=""), tolower(file))
  if(result[[1]] != -1) {
      return(CNMF.read.res(file))
  }
  stop("Input is not a res or gct file.")
}


CNMF.read.gct <- function(filename = "NULL") {
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T,check.name=F)
  ds <- ds[-1]
  return(ds)
}

CNMF.read.res <- function(filename = "NULL") {
  # Reads a gene expression dataset in RES format and converts it into an R data frame
  header.cont <- readLines(filename, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 2)]
  ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
  colst <- length(ds[1,])
  cols <- (colst - 1)/2
  rows <- length(ds[,1])
  A <- matrix(nrow=rows - 1, ncol=cols)
  A <- ds[1:rows, seq(2, colst, 2)]
  table1 <- data.frame(A)
  names(table1) <- header.labels
  return(table1)
}

CNMF.write.gct <- function (gct, filename) {
  # Write gct formated file to filename
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct)[1], "     ", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "       ", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  names <- names(gct)
  cat("       ", names[1], file = f, append = TRUE, sep = "")
  for (j in 2:length(names)) {
      cat("   ", names[j], file = f, append = TRUE, sep = "")
  }
  cat("\n", file = f, append = TRUE, sep = "")
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] +  2)
  m[, 1] <- row.names(gct)
  m[, 2] <- row.names(gct)
  index <- 3
  for (i in 1:dim(gct)[2]) {
      m[, index] <- gct[, i]
      index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "      ", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
  return(gct)
}

