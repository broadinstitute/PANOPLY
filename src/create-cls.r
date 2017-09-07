

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
    tumor.group <- tumor.info [ sample.order, g ]
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
# use exptdesign file if it exists; else use the master data file (normalized+filtered)
create.cls (type,
            expt.design=ifelse (file.exists (expt.design.file), expt.design.file, 
                                paste (master.prefix, 'gct', sep='.')))

