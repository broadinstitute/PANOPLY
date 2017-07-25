

source ('preamble.r')


sd.threshold <- 1
out.prefix <- 'rna-seq'



d <- parse.gctx ( file.path (data.dir, rna.data.file) )
d@rdesc [,'Description'] <- d@rdesc [,'id']   # make id and Description identical
samples <- d@cid

tumor.info <- read.csv (file.path (data.dir, 'exptdesign.csv'))
samples.order <- tumor.info[,'Sample.ID']
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
sds <- apply (d@mat, 1, sd, na.rm=TRUE)
d.subset <- row.subset.gct (d, sds > sd.threshold)
write.gct (d.subset, paste (out.prefix, '-sdfilter.gct', sep=''))

