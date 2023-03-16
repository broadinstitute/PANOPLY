#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')
Source ('normalization.r')


if (normalize.proteomics) { # if normalization is on
  file.prefix <- paste (type, '-ratio', sep='')
  
  ## Check if qc.col exists
  ds=parse.gctx(file.path (pre.dir, paste (file.prefix, '.gct', sep='')))
  if ( is.null(ds@cdesc[[qc.col]]) ) { # if qc.col doesn't exists in ds@cdesc
    qc.col=NULL # set qc.col as NULL for the sake of normalization
  }
  
  ## Run Normalization and plot results
  normalize.dataset (input.gct=file.path (pre.dir, paste (file.prefix, '.gct', sep='')),
                     out.prefix=file.prefix, 
                     method=norm.method,
                     qc.col=qc.col, ndigits=ndigits)
  
  if (!is.null (alt.method)) 
    normalize.dataset (input.gct=file.path (pre.dir, paste (file.prefix, '.gct', sep='')),
                       out.prefix=paste (type, '-ratio-', alt.method, sep=''), 
                       method=alt.method,
                       qc.col=qc.col, ndigits=ndigits)
} else { # otherwise
  # copy old file to new file location
  print( "No normalization applied." )
  file.prefix <- paste (type, '-ratio', sep='')
  file.copy( file.path (pre.dir, paste (file.prefix, '.gct', sep='')), paste (file.prefix, '-norm.gct', sep='') )
}

