

source ('config.r')


sd.threshold <- rna.sd.threshold
out.prefix <- rna.output.prefix



d <- parse.gctx (rna.data.file)
d@rdesc [,'Description'] <- d@rdesc [,'id']   # make id and Description identical
samples <- d@cid

if (file.exists(expt.design.file)) {
  # obtain sample order from experiment design file
  tumor.info <- read.csv (expt.design.file)
  samples.order <- tumor.info[,'Sample.ID']
} else {
  # obtain sample order from master data file (gct3)
  filtered.data <- parse.gctx (master.file)
  samples.order <- filtered.data@cdesc [, 'Sample.ID']
}

run.order <- sapply (samples.order,
                     function (x) {
                       pos <- which (samples %in% x)
                       # missing samples point to NA column add below
                       return ( ifelse (length(pos)==0, length(samples)+1, pos) )
                     })
d <- add.cols.gct (d, data.frame (na.col=rep(NA, nrow(d@mat))))   # add NA column at the end -- use for missing samples
d <- rearrange.gct (d, run.order, new.cid=samples.order)
write.gct (d, paste (out.prefix, '.gct', sep=''))


## filters
# exclude rows with SD less than sd.threshold
if ( !is.na (sd.threshold) ) {
  sds <- apply (d@mat, 1, sd, na.rm=TRUE)
  d.subset <- row.subset.gct (d, sds > sd.threshold)
  write.gct (d.subset, paste (out.prefix, '-sdfilter.gct', sep=''))
}

