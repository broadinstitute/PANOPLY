
source ('config.r')
Source ('normalization.r')

## Run Normalization and plot results
file.prefix <- paste (type, '-ratio', sep='')
normalize.dataset (input.gct=file.path (pre.dir, paste (file.prefix, '.gct', sep='')),
                   out.prefix=file.prefix, 
                   method=norm.method,
                   qc.col=qc.col, ndigits=ndigits)

if (!is.null (alt.method)) 
  normalize.dataset (input.gct=file.path (pre.dir, paste (file.prefix, '.gct', sep='')),
                     out.prefix=paste (type, '-ratio-', alt.method, sep=''), 
                     method=alt.method,
                     qc.col=qc.col, ndigits=ndigits)
