
# $Id$


## NOTE: If using with gct-io.r, SOURCE gct-io.r FIRST
##      Functions names for write.gctv2 and read.gctv2 will be decided appropriately

write.gct.exists <- existsFunction ('write.gct')


write.gctv2 <- function (data, file, ...) {
  # write out gct file
  version <- '#1.2'
  size <- dim (data)
  size[2] <- size[2] - 2  # adjust for "Name" & "Description"
  write.table (version, file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  write.table (rbind(size), file, append=TRUE, quote=FALSE, sep='\t',
               row.names=FALSE, col.names=FALSE, ...)
  write.table (data, file, append=TRUE, quote=FALSE, sep='\t', row.names=FALSE, ...)
}
if (write.gct.exists) write.gct2 <- write.gctv2 else write.gct <- write.gctv2


write.cls <- function (cls, file, respect_factor_order = FALSE) {
  # write out cls file

  # create a new cls vector -- otherwise assignments may result in <NA> if cls is a factor
  orig.cls <- as.character (cls)   # de-factor to handle NA items
  orig.cls [ is.na(orig.cls) ] <- "NA"

  classes_unique = sort ( unique (cls) ) # get unique classes
  if(respect_factor_order) {
    fct_order = levels(cls) # get factor levels
    classes = c(intersect(fct_order, classes_unique), # retain order of original classes, but only take those that appear in data
                setdiff(classes_unique,fct_order)) # add any classes that were missing
  } else { classes <- classes_unique } # otherwise, take resorted unique class array as factor order

  line1 <- c ( length(cls), length(classes), 1 )
  line2 <- paste ( '#',
                   paste (classes, collapse=" ") )

  cls <- rep ( 0, length(orig.cls) )
  for ( i in 1:length(classes) )
    cls [ orig.cls==classes[i] ] <- i-1

  write.table (rbind(line1), file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table (line2, file, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table (cls, file, eol=' ', append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)  # single line
  write.table ('', file, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)            # end above line 
}


read.gctv2 <- function (file, preproc=FALSE, ...) {
  # reads a data file in gct format
  data <- read.table (file, header=TRUE, sep='\t', skip=2,
                      comment.char='', quote='"', ...)
  if (preproc) {
    # preprocess data by setting rownames and eliminating the information columns
    rownames (data) <- data [,1]
    data <- data [, c(-1,-2)]
  }
  invisible (data)
}
if (write.gct.exists) read.gct2 <- read.gctv2 else read.gct <- read.gctv2


read.res <- function (file, ...) {
  # read a data file in res format, and creates a data frame after
  # stripping the A/P call columns
  col.names <- scan (file=file, what=list("string"), nmax=1, sep='\t')[[1]]
  n.cols <- length(col.names)
  if ( n.cols %% 2 == 0 ) {
    # res file with an extra \t at the end of each line
    # (usually from res files created within GeneCluster)
    table <- read.table (file, header=FALSE, sep='\t', skip=3) [, 1:n.cols]
    colnames (table) <- col.names
  } else {
    # res file with no extra \t at the end of each line
    col.names <- c (col.names, "") # scan skips the last empty field for A/P
    table <- read.table (file, header=FALSE, sep='\t', col.names=col.names, skip=3)
  }
    data <- table [, c(1,2,seq(3,ncol(table),2))]

  invisible (data)
}


read.cls <- function (file, labels=TRUE, respect_factor_order = FALSE, ...) {
  # reads a cls file
  cls <- scan (file, skip=2, quiet=TRUE, ...)
  if (labels) {
    cls.labels <- scan (file, skip=1, "string", nlines=1) [-1]    # [-1] to remove '#' at start 
    if (length (cls.labels) > 0) {
      for (i in (length(cls.labels)-1):0) cls [cls==i] <- cls.labels [i+1]
      if (respect_factor_order) cls <- factor(cls, levels=cls.labels) # optionally convert to factor
    }
  }
  invisible (cls)
}


read.cls.mapped <- function (file, ...) {
  # reads a cls file, and maps the classes to names listed in the file
  stream <- file (file, open='r')
  line1 <- scan (stream, what=list(integer(0)), nmax=1, ...)
  names <- scan (stream, what=list('string'), nmax=1, ...)
  cls <- scan (stream, ...)
  close (stream)

  n.samples <- line1[[1]][1]
  n.classes <- line1[[1]][2]
  for (i in 0:(n.classes-1)) cls [ cls==i ] <- names[[1]][i+2]  # names[[1]] == '#'
  invisible (cls)
}



permute.cls <- function (file, ...) {
  # generates a random permutation of the cls vector in file
  cls <- read.cls (file, ...)
  write.cls (sample (cls), file)
}


read.dataset <- function (file, read.cls=TRUE, name.prefix="MZ.") {
  # reads both gct and cls file, transposes gct file and
  # returns dataset and class vector
  # reading class vector is skipped if read.cls=FALSE
  
  d <- read.gct2 ( paste(file,'.gct',sep="") )
  d.t <- t ( d[, 3:dim(d)[2]] )
  names <- unlist ( lapply (d[,'Name'],
                            function (x) { paste (name.prefix, toString(x), sep="") }) )
  colnames (d.t) <- names

  if (read.cls) {
    c <- read.cls ( paste(file,'.cls',sep="") )
    c.t <- as.factor (c)
  } else c.t <- NULL

  invisible ( list (data=data.frame(d.t), class.vector=c.t) )
}





# ODF (GenePattern) format file processing

read.odf <- function (file) {
  # reads GenePattern ODF files
  # most of the header information is discarded;
  # only the column names are slavaged

  
  con <- file (file, open="r")
  # header lines
  line1 <- scan (con, what=list('string'), nmax=1)
  if ( line1[[1]][1] != 'ODF' ) {
    close (con)
    stop ( paste ("Invalid ODF format:",file) )
  }
  line2 <- scan (con, what=list('string',integer(0)), sep='=', nmax=1)
  header.lines <- line2[[2]][1]
  # column names
  line.num <- 2
  col.names <- NULL
  repeat {
    line.num <- line.num + 1
    if (line.num > header.lines+2) break
    line <- scan (con, what=list('string'), nmax=1, sep='\t')[[1]]
    CN <- "COLUMN_NAMES:"
    if ( substr (line[1], 1, nchar(CN)) == "COLUMN_NAMES:" ) {
      line[1] <- substr (line[1], nchar(CN)+1, nchar(line[1]))
      col.names <- make.names (line)
    }
  }
  if ( is.null(col.names) ) {
    close (con)
    stop ( paste ("Can't find COLUMN_NAMES:", file) )
  }

  # read table
  table <- read.table (con, header=FALSE, sep='\t', col.names=col.names, as.is=TRUE)
  close (con)

  return (table)
}


