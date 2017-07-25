
source ('preamble.r')


filter.dataset <- function (file.prefix, numratio.file=NULL, out.prefix=NULL,
                            na.max=NULL, no.na=TRUE, separate.QC.types=TRUE, cls=NULL,
                            min.numratio=1, n.min.numratio=NULL, sd.threshold=NULL) {
  
  write.out <- function (d, f) {
    # ds@mat <- d
    # if (separate.QC.types && length (unique (cls)) > 1)   # if there is more than one type
    #   for (cl in unique (cls)) {
    #     d.t <- subset.gct (ds, cls==cl)
    #     if (ncol(d.t@mat) > 0) write.gct (d.t, paste (f, '-', cl, '.gct', sep=''), precision=ndigits)
    #   }
    # else write.gct (ds, paste (f, '.gct', sep=''), precision=ndigits)
    if (separate.QC.types && length (unique (cls)) > 1)   # if there is more than one type
      for (cl in unique (cls)) {
        d.t <- subset.gct (d, cls==cl)
        if (ncol(d.t@mat) > 0) write.gct (d.t, paste (f, '-', cl, '.gct', sep=''), precision=ndigits)
      }
    else write.gct (d, paste (f, '.gct', sep=''), precision=ndigits)
    
  }
  
  
  if (is.null (out.prefix)) out.prefix <- file.prefix
  
  
  ds <- parse.gctx ( paste (file.prefix, '.gct', sep='') )
  # replace Description with gene names
  map <- read.delim ( file.path (data.dir, protein.gene.map) )
  id.table <- merge ( cbind (ds@rdesc, refseq_protein=sub ("(.*?)\\..?", "\\1", ds@rdesc[,'id'])), map,
                  by='refseq_protein', all.x=TRUE, sort=FALSE )
  ds@rdesc [,'Description'] <- as.character (id.table [,'gene_name'])
  
  #data <- ds@mat
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
    numratios <- parse.gctx (numratio.file)
    d.numratios <- numratios@mat
    if ( is.null (n.min.numratio) ) min.numratio.count <- ifelse (is.null (na.max), ncol (d.numratios), na.max )
    else min.numratio.count <- n.min.numratio
    ds <- row.subset.gct (ds, unlist (apply (d.numratios, 1, 
                                             function (x) sum (x < min.numratio, na.rm=TRUE) < min.numratio.count)))
    # data <- data [ unlist (apply (d.numratios, 1, 
    #                               function (x) sum (x < min.numratio, na.rm=TRUE) < min.numratio.count)), ]
  }
  
  if (!is.null (sd.threshold))
    # exclude rows with SD less than sd.threshold
    ds <- row.subset.gct (ds, apply (ds@mat, 1, sd, na.rm=TRUE) > sd.threshold)
    # data <- data [ apply (data, 1, sd, na.rm=TRUE) > sd.threshold, ]

  if (! is.null (na.max)) {
    # keep rows with no more than na.max missing values
    ds <- row.subset.gct (ds, apply (ds@mat, 1, function (x) sum (!is.finite (x)) < na.max))
    # data <- data [ apply (data, 1, function (x) sum (!is.finite (x)) < na.max), ]
    write.out (ds, paste (out.prefix, '-NAmax', na.max, sep=''))  
  } 
  
  if (no.na) {
    # keep only if observed in all samples: exclude any rows with missing values
    # this must occur after  if (! is.null (na.max))  ...
    ds <- row.subset.gct (ds, apply (ds@mat, 1, function (x) all (is.finite (x))) )
    # data <- data [ apply (data, 1, function (x) all (is.finite (x))), ]
    write.out (ds, paste (out.prefix, '-noNA', sep=''))
  }

}


## apply filtering
# separate bimodal and unimodal samples (no effect here), with SD filter
filter.dataset (paste (type, '-ratio-norm', sep=''), 
                na.max=na.max, min.numratio=min.numratio, sd.threshold=sd.filter.threshold)
# do not separate bi/unimodal samples (no effect here) -- no SD filter (to retain max proteins, for mRNA correlation)
filter.dataset (paste (type, '-ratio-norm', sep=''), 
                out.prefix=paste (type, '-ratio-norm-nosdfilter', sep=''),
                na.max=na.max, no.na=FALSE, min.numratio=min.numratio, sd.threshold=NULL,
                separate.QC.types=FALSE)

