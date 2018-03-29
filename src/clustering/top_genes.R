library(optparse)
options(stringsAsFactors=FALSE)

# Determine the top genes according to MAD (Median Absolute Deviation) ranking across samples.
#
# @param expfile A tab-delimited expression file
# @param output_prefix A prefix to use on generated output
# @param selected_genes A value indicating gene selection criteria, ALL (all genes), Auto (findpeaks method), or
#                       a numeric value indicating the fraction of available genes to select.

option.list <- list(
  make_option(c("-m", "--expfile"), help="A tab delimited expression file"),
  make_option(c("-o", "--output_prefix"), help="A prefix to use on generated output"),
  make_option(c("-u", "--selected_genes"), help="A value indicating gene selection criteria"),
  make_option(c("-y", "--yaml"), help="Optional yaml file specifying some or all config options. Precedence given to CL.")
)

main <- function(opt) {
  exp.file       <- opt$expfile
  output.prefix  <- opt$output_prefix
  selected.genes <- opt$selected_genes

  exp <- read.data(opt$expfile)

  # Remove up to 2 leading non-numeric columns.
  if (all(is.na(as.numeric(exp[,1])))) {
    exp <- exp[, -1, drop=FALSE]
    if (all(is.na(as.numeric(exp[,1])))) {
      exp <- exp[, -1, drop=FALSE]
    }
  }

  if(dim(exp)[1] > length(unique(rownames(exp)))) {
    tnames <- rownames(exp)[duplicated(rownames(exp))]
    tt     <- which(rownames(exp) %in% tnames)
    exp    <- exp[-tt,]
  }

  exp[exp == -Inf] <- NA

  exp.rownames <- rownames(exp)

  exp <- apply(exp, 2, as.numeric)

  exp <- sweep(exp, 1, apply(exp, 1, median, na.rm=TRUE))

  # colnames(exp) <- gsub("\\.","-",colnames(exp))

  rownames(exp) <- exp.rownames

  # filter columns that are 80% NA
  na.number  <- apply(exp, 2, function(x) sum(is.na(x)))
  t          <- which(na.number > 0.8 * dim(exp)[1])

  if(length(t) > 0 ) {
    message(paste0("[Filter] removed ", length(t), " columns more than 80% NA."))
    exp <- exp[,-t]
  }

  # filter rows that are 70% NA
  na.number  <- apply(exp, 1, function(x) sum(is.na(x)))
  t          <- which(na.number >= 0.7 * dim(exp)[2])

  if(length(t) > 0) {
    message(paste0("[Filter] removed ", length(t), " rows at least or more than 70% NA."))
    exp <- exp[-t, ]
  }

  if (selected.genes == "ALL") {
    convert2gct(exp, output.prefix)
  } else {
    if (selected.genes < 1) {
      selected.genes <- as.numeric(selected.genes) * nrow(exp)
    }
    top.exp <- top_genes(as.data.frame(exp), selected.genes)
    convert2gct(top.exp, output.prefix)
  }
}

convert2gct <- function(exp, output_prefix) {
  # gct format matrix using gene name and description columns
  exp <- cbind("GeneName"    = as.character(rownames(exp)),
               "DESCRIPTION" = as.character(rownames(exp)),
               exp)

  # should not be required
  t <- which( is.na(exp[,1]) )
  if (length(t) > 0 ) {
    exp <- exp[-t,]
  }

  gct.rows <- nrow(exp)
  gct.cols <- ncol(exp)-2
  # write out gct file
  output.fn <- paste0(output_prefix, ".expclu.gct")
  write_tsv("#1.2", filename=output.fn, col.names=FALSE)
  write_tsv(paste(gct.rows, gct.cols, sep="\t"), filename=output.fn, col.names=FALSE, append=TRUE)
  write_tsv(exp, filename=output.fn, append=TRUE)
  print(paste("\n\nWrote", paste0(output_prefix,".expclu.gct")))
  print(setNames(data.frame(t(dim(exp))), c("Rows","Columns")))
}


top_genes <- function(unified.scale, selected.genes) {
  mads  <- sort(apply(unified.scale, 1, mad, na.rm=TRUE), decreasing=TRUE, index.return=TRUE)
  count <- 1
  if (selected.genes == "Auto") {
    selected.genes <- find_peaks(mads$x)
  } else {
    selected.genes <- as.numeric(selected.genes)
  }
  dataset <-matrix(nrow=0, ncol=ncol(unified.scale))
  while (count <= selected.genes){
    tmp     <- unified.scale[mads$ix[count], ]
    dataset <- rbind(dataset, tmp )
    count   <- count + 1
  }
  return(dataset)
}


find_peaks <- function(tmad) {
  den <- density(tmad)
  id  <- which(den$y == max(den$y))
  X   <- den$x
  Y   <- den$y
  slope <- sapply((id+1):(length(X)-1),function(i){
    tslope <- (Y[i+1]-Y[i])/(X[i+1]-X[i])
    return(tslope)
  })
  t <- which(slope >= 0)
  if(length(t) == 0) {
    t    <- which(tmad >= quantile(tmad,0.9))
    tlen <- length(t)
  } else {
    peak <- X[id+t[1]-1]
    tid  <- which(tmad >= X[id+t[1]-1])
    tlen <- length(tid)
  }
  return(tlen)
}

read.data <- function(fn, suffix.dup=FALSE) {
  # loads data from one of two formats
  check  <- readLines(fn, n=1)
  check0 <- strsplit(check, '\t')[[1]][1]
  if (check0 == '#1.2') {
    check <- readLines(fn, n=3)[3]
    mat   <- read.delim(fn, as.is=TRUE, skip=3, header=FALSE)
    gnms  <- toupper(mat[,1])
    mat   <- mat[, -(1:2),drop=FALSE]
    dnms  <- toupper(strsplit(check, '\t')[[1]][-(1:2)])
  } else {
    check2 <- readLines(fn, n=2)[2]
    check2 <- toupper(strsplit(check2, '\t')[[1]])
    nskip  <- ifelse(grepl('COMPOSITE *ELEMENT REF', check2[1]), 2, 1)
    mat    <- read.delim(fn, as.is=TRUE, skip=nskip, header=FALSE)
    gnms   <- toupper(mat[,1])
    mat    <- mat[, -1, drop=F]
    dnms   <- toupper(strsplit(check, '\t')[[1]][-1])
  }
  # fix excel gene errors
  ind <- grep('^\\d+-SEP$', gnms)
  if(length(ind) > 0) {
    tmp <- gnms[ind]
    tmp <- sub('-SEP','',tmp)
    tmp <- sprintf('SEPT%s', tmp)
    gnms[ind] <- tmp
  }
  ind <- grep('^\\d+-MAR$', gnms)
  if(length(ind) > 0) {
    tmp <- gnms[ind]
    tmp <- sub('-MAR','',tmp)
    tmp <- sprintf('MARCH%s', tmp)
    gnms[ind] <- tmp
  }
  ind <- grep('^\\d+-APR$', gnms)
  if(length(ind) > 0) {
    tmp <- gnms[ind]
    tmp <- sub('-APR','',tmp)
    tmp <- sprintf('APR-%s', tmp)
    gnms[ind] <- tmp
  }
  mat <- as.matrix(mat)
  # not sure about the required use-case
  if (suffix.dup) {
    dup  <- rep(0, nrow(mat))
    gnm2 <- gnms
    ii   <- which(duplicated(gnm2))
    dup[ii] <- dup[ii]+1
    while(length(ii)>0) {
      gnm2 <- paste(gnms, dup, sep='__')
      ii   <- which(duplicated(gnm2))
      dup[ii] <- dup[ii]+1
    }
    gnm2 <- sub('\\__0$', '', gnm2)
    rownames(mat) <- gnm2
  } else {
    rownames(mat)=gnms
  }

  colnames(mat)=dnms
  return(mat)
}

write_tsv <- function(data, filename, col.names=TRUE, append=FALSE, fmt=FALSE, ...) {
  if(fmt) {
    is.num        <- sapply(data, is.numeric)
    data[is.num]  <- lapply(data[is.num], round, 6)
  }
  suppressWarnings(write.table(data, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=col.names, append=append, ...))
}

opts  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

main(opts)

