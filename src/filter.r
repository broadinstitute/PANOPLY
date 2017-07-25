
source ('config.r')


filter.dataset <- function (file.prefix, numratio.file=NULL, out.prefix=NULL,
                            na.max=NULL, no.na=TRUE, separate.QC.types=TRUE, cls=NULL,
                            min.numratio=1, n.min.numratio=NULL, sd.threshold=NULL) {
  
  write.out <- function (d, f) {
    ds@mat <- d
    if (separate.QC.types && length (unique (cls)) > 1)   # if there is more than one type
      for (cl in unique (cls)) {
        d.t <- subset.gct (ds, cls==cl)
        if (ncol(d.t@mat) > 0) write.gct (d.t, paste (f, '-', cl, '.gct', sep=''),
                                          ver=3, precision=ndigits, appenddim=FALSE)
      }
    else write.gct (ds, paste (f, '.gct', sep=''),
                    ver=input.ver, precision=ndigits, appenddim=FALSE)
  }
  
  
  if (is.null (out.prefix)) out.prefix <- file.prefix
  
  
  ds <- parse.gctx ( paste (file.prefix, '.gct', sep='') )
  data <- ds@mat
  input.ver <- ifelse (ds@version=="#1.3", 3, 2)
  # cls for QC types -- input cls takes precedence
  # if input not specified, check if qc.col exisits in gct3 file
  # if none are specified, generate one with single class
  if (is.null (cls)) {
    if (qc.col %in% colnames (ds@cdesc)) cls <- ds@cdesc [,qc.col]
    else cls <- rep ('QC.pass', ncol (ds@mat))
  }
  
  
  ## filters
  ##
  # exclude rows with numratio less than threshold in too many samples
  # d.numratios matches data -- this must be the first filter
  if (! is.null (numratio.file)) {
    # only if numratio file is specified
    numratios <- parse.gctx (numratio.file)
    d.numratios <- numratios@mat
    if ( is.null (n.min.numratio) ) min.numratio.count <- ifelse (is.null (na.max), ncol (d.numratios), na.max )
    else min.numratio.count <- n.min.numratio
    data <- data [ unlist (apply (d.numratios, 1, 
                                  function (x) sum (x < min.numratio, na.rm=TRUE) < min.numratio.count)), ]
  }
  
  if (!is.null (sd.threshold))
    # exclude rows with SD less than sd.threshold
    data <- data [ apply (data, 1, sd, na.rm=TRUE) > sd.threshold, ]

  if (! is.null (na.max)) {
    # keep rows with no more than na.max missing values
    data <- data [ apply (data, 1, function (x) sum (!is.finite (x)) < na.max), ]
    write.out (data, paste (out.prefix, '-NAmax', na.max, sep=''))  
  } 
  
  if (no.na) {
    # keep only if observed in all samples: exclude any rows with missing values
    # this must occur after  if (! is.null (na.max))  ...
    data <- data [ apply (data, 1, function (x) all (is.finite (x))), ]
    write.out (data, paste (out.prefix, '-noNA', sep=''))
  }

}


## apply filtering
nr.file <- ifelse (apply.SM.filter, 
		   file.path (pre.dir, paste (type, '-num-ratio.gct', sep='')),
		   NULL)
# separate bimodal and unimodal samples (no effect here), with SD filter
filter.dataset (paste (type, '-ratio-norm', sep=''), 
                numratio.file=nr.file,
                na.max=na.max, min.numratio=min.numratio, sd.threshold=sd.filter.threshold)
# do not separate bi/unimodal samples (no effect here) -- no SD filter (to retain max proteins, for mRNA correlation)
filter.dataset (paste (type, '-ratio-norm', sep=''), 
                numratio.file=nr.file,
                out.prefix=paste (type, '-ratio-norm-nosdfilter', sep=''),
                na.max=na.max, no.na=FALSE, min.numratio=min.numratio, sd.threshold=NULL,
                separate.QC.types=FALSE)

