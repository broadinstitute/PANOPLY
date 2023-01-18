#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')
pacman::p_load (dplyr)  # transition to dplyr to make data wranging more consistent and bug-free


filter.dataset <- function (file.prefix, numratio.file=NULL, out.prefix=NULL,
                            na.max=NULL, no.na=TRUE, separate.QC.types=TRUE, cls=NULL,
                            min.numratio=1, n.min.numratio=NULL, sd.threshold=NULL,
                            combine.replicates=NULL) {
  ## filter input dataset
  ## separate QC pass/fail by default
  
  write.out <- function (f) {
    # write out current ds (gct) to file f
    if (separate.QC.types && length (unique (cls)) > 1)   # if there is more than one type
      for (cl in unique (cls)) {
        d.t <- col.subset.gct (ds, cls==cl)
        if (ncol(d.t@mat) > 0) {
          if (cl == qc.pass.label) file.name <- paste (f, '.gct', sep='')  # QC passes cases
          else file.name <- paste (f, '-', cl, '.gct', sep='')
          write.gct (d.t, file.name, ver=3, precision=ndigits, appenddim=FALSE)
        }
      }
    else write.gct (ds, paste (f, '.gct', sep=''),
                    ver=input.ver, precision=ndigits, appenddim=FALSE)
  }
  
  
  if (is.null (out.prefix)) out.prefix <- file.prefix
  
  
  ds <- parse.gctx ( file.path( norm.dir, paste (file.prefix, '.gct', sep='') ) )
  input.ver <- ifelse (ds@version=="#1.3", 3, 2)
  # replace Description with gene names
  if (input.ver == 3 && any (grepl ('gene[.-_]?(name|id|symbol)s?$', colnames (ds@rdesc), ignore.case=TRUE))) {
    # gene symbol is already present as an annotation column
    genesym.col <- grep ('gene[.-_]?(name|id|symbol)s?$', colnames (ds@rdesc), ignore.case=TRUE)
    if (length (genesym.col) > 1) {
      genesym.col <- genesym.col[1]
      warning ( paste ('Identified multiple gene symbol columns. Using', genesym.col) )
    }
    ds@rdesc [,gene.id.col] <- as.character (ds@rdesc [,genesym.col])
  } else if (input.ver == 3 && any (grepl (gene.id.col, colnames (ds@rdesc), ignore.case=TRUE))) {
    # gene symbol is already present as an annotation column with name given in yaml parameters file
    genesym.col <- grep (gene.id.col, colnames (ds@rdesc), ignore.case=TRUE)
    if (length (genesym.col) > 1) {
      genesym.col <- genesym.col[1]
      warning ( paste ('Identified multiple gene symbol columns. Using', genesym.col) )
    }
    ds@rdesc [,gene.id.col] <- as.character (ds@rdesc [,genesym.col])

  } else {
    # map protein id to gene symbols
    map <- read.delim ( file.path (data.dir, protein.gene.map) )
    id.table <- left_join ( cbind (ds@rdesc, refseq_protein=sub ("(.*?)\\..?", "\\1", ds@rdesc[,'id'])), map,
                            by='refseq_protein' )
    ds@rdesc [,gene.id.col] <- as.character (id.table [,'gene_name'])
  }  
  
  
  # convert na.max to integer if a fraction 
  if (na.max < 1) na.max <- ceiling (ncol(ds@mat) * na.max)
  
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
    if ( is.null (n.min.numratio) ) {
      # normally, ALL samples must have num ratio >= min.numratio
      keep <- unlist (apply (d.numratios, 1, function (x) all (x >= min.numratio, na.rm=TRUE)))
    } else {
      # if n.min.numratio is specified, then at least that many samples should have num ratios >= min.numratios
      if (n.min.numratio < 1) n.min.numratio <- ceiling (ncol(ds@mat) * n.min.numratio)    # convert to integer if a fraction
      keep <- unlist (apply (d.numratios, 1, 
                             function (x) sum (x >= min.numratio, na.rm=TRUE) >= n.min.numratio))
    }
    ds <- row.subset.gct (ds, index=keep)
  }
  
  if (!is.null (sd.threshold))
    # exclude rows with SD less than sd.threshold
    ds <- row.subset.gct (ds, index=apply (ds@mat, 1, sd, na.rm=TRUE) > sd.threshold)

  if (!is.null (combine.replicates) && input.ver == 3) {
    # combine replicate samples; possible options -- first, mean, median
    # replicates are identified by identical (Participant,Type) 
    # or (Participant) only if Type does not exist
    # executed only for GCT v1.3 inputs
    avail.cols <- intersect (colnames (ds@cdesc), c ('Participant', 'Type'))
    if (length (avail.cols) > 0 && 'Participant' %in% avail.cols) {
      dup.table <- ds@cdesc %>% mutate (sample.num=row_number()) %>% 
        mutate (group.num=group_indices (., !!!rlang::syms (avail.cols))) %>%
        group_by (group.num) %>% filter (n() > 1)
      # alternative for above suggested by Karl to avoid deprecated warning
      # dup.table <- ds@cdesc %>% mutate (sample.num=row_number()) %>%
      #   group_by (!!!rlang::syms (avail.cols)) %>%  #!!!rlang to get actual content of vector: Participant, Type (if present)
      #   mutate(group.num = cur_group_id()) %>%      #dplyr 1.0 (R4.1.0) cur_group_id() replaces group_indices()
      #   filter (QC.status == "QC.pass") %>%         #if needed 
      #   filter (n() > 1)
      for (g in unlist (unique (dup.table [,'group.num']))) {
        rows <- unlist (dup.table [ dup.table[,'group.num']==g, 'sample.num'])
        # replace first occurance of replicate with combined value
        ds@mat[,rows[1]] <- apply (ds@mat[,rows], 1, function (x) { 
          switch (combine.replicates,
                  "first"=x[1],
                  "mean"=mean (x, na.rm=T),
                  "median"=median (x, na.rm=T)) })
      }
      # keep first occurance only, and fix cls
      keep.samples <- !duplicated (ds@cdesc[, avail.cols])
      ds <- col.subset.gct (ds, keep.samples)   
      cls <- cls [keep.samples]
    }    
  }
  
  if (! is.null (na.max)) {
    # keep rows with no more than na.max missing values
    ds <- row.subset.gct (ds, index=apply (ds@mat, 1, function (x) sum (!is.finite (x)) < na.max))
    write.out (paste (out.prefix, '-NArm', sep=''))  
  } 
  
  if (no.na) {
    # keep only if observed in all samples: exclude any rows with missing values
    # this must occur after  if (! is.null (na.max))  ...
    ds <- row.subset.gct (ds, index=apply (ds@mat, 1, function (x) all (is.finite (x))))
    write.out (paste (out.prefix, '-noNA', sep=''))
  }

}


## apply filtering
# enable num-ratio filter only if apply.SM.filter is set and the file is present
nr.file.path <- file.path (parse.dir, paste (type, '-num-ratio.gct', sep=''))
if (apply.SM.filter && file.exists(nr.file.path)) {
  nr.file <- nr.file.path
} else {
  nr.file <- NULL
}

## Change default to No SD filter -- filter removes very few items
filter.dataset (paste (type, '-ratio-norm', sep=''),
                numratio.file=nr.file,
                na.max=na.max, min.numratio=min.numratio, sd.threshold=NULL,
                combine.replicates='mean', n.min.numratio=min.numratio.fraction)

# # With SD filter
# filter.dataset (paste (type, '-ratio-norm', sep=''), 
#                 numratio.file=nr.file,
#                 na.max=na.max, min.numratio=min.numratio, sd.threshold=sd.filter.threshold,
#                 combine.replicates='mean', n.min.numratio=min.numratio.fraction)
# # No SD filter (to retain max proteins, for mRNA correlation)
# filter.dataset (paste (type, '-ratio-norm', sep=''), 
#                 numratio.file=nr.file,
#                 out.prefix=paste (type, '-ratio-norm-nosdfilter', sep=''),
#                 na.max=na.max, no.na=FALSE, min.numratio=min.numratio, sd.threshold=NULL,
#                 combine.replicates='mean', n.min.numratio=min.numratio.fraction)

