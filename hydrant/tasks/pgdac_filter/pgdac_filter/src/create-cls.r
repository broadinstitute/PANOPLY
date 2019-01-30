

source ('config.r')



create.cls <- function (out.prefix, expt.design, type.cls=NULL) {
  # expt.design must be either a gct3 file, or a data frame with sample annotations
  
  # get proper column names and class information
  if ( (grepl(".gct$", expt.design) || grepl(".gctx$", expt.design)) && 
       scan (file=expt.design, what="string", n=1)=="#1.3" ) {
    # input file is gct3
    ds <- parse.gctx (expt.design)
    tumor.info <- ds@cdesc
    sample.order <- ds@cid
    # create expt.design file for use downstream in the pipeline
    expt.desn <- data.frame (ds@cdesc)
    write.csv (expt.desn, expt.design.file, row.names=FALSE, quote=FALSE)
  } else {
    tumor.info <- read.csv (expt.design)
    rownames (tumor.info) <- tumor.info[,'Sample.ID']
    
    sample.order <- scan (paste (master.prefix, 'gct', sep='.'),
                          what='char', nlines=1, skip=2, sep='\t')[c(-1,-2)]
  }
  for (g in setdiff (colnames (tumor.info), c('id', 'Sample.ID'))) {
    tumor.group <- make.names (tumor.info [ sample.order, g ])  # use make.names to convert text to proper class labels
    n.groups <- nlevels (factor (tumor.group))                  # and ensure there are a reasonable number of groups
    if (n.groups < 2 || n.groups > 5) next
    # global class labels
    write.cls (tumor.group, paste (out.prefix, '-', g, '.cls', sep=''))
    
    # class labels separated by uni/bi-modal groups
    if (! is.null (type.cls)) {
      cls <- read.cls (type.cls)
      for (type in unique (cls)) {
        t.cls <- cls == type
        if (any (t.cls)) {
          out <- paste (out.prefix, '-', type, sep='')
          # write all class labels
          write.cls (tumor.group[t.cls], paste (out, '-', g, '.cls', sep=''))
        }
      }  
    }
  }
  
  
}



# different "type"s result in the same cls files, 
# but different file prefix is convenient for batch scripts
# use the master data file (normalized+filtered) to account for filtered samples
create.cls (type,
            expt.design=paste (master.prefix, 'gct', sep='.'))

