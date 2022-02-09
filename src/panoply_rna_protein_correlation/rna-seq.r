#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


source ('config.r')


sd.threshold <- rna.sd.threshold
out.prefix <- rna.output.prefix



d <- parse.gctx (rna.data.file)
d@rdesc[,gene.id.col] <- d@rdesc [,'Description'] <- d@rdesc [,'id']   # make id, Description and gene.id.col identical
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

if (rna.row.norm.method == "mean") { # normalize sample mrna to a 'common reference' expression level for each gene
  center <- apply (d@mat, 1, mean, na.rm=TRUE) # compute mean for each row (i.e. gene), excluding NA values
} else if (rna.row.norm.method == "median") {
  center <- apply (d@mat, 1, median, na.rm=TRUE) # compute median for each row (i.e. gene), excluding NA values
} else {
  warning("The method specified with 'rna.row.norm.method' is not a valid option. Defaulting to median.")
  center <- apply (d@mat, 1, median, na.rm=TRUE) #or else write NA
}
d@mat <- sweep(d@mat, 1, center) #normalize data according to chosen center (subtracted out)

write.gct (d, paste (out.prefix, '.gct', sep=''), ver=3, appenddim=FALSE)  # write out GCT v1.3

## filters
# exclude rows with SD less than sd.threshold
if ( !is.na (sd.threshold) ) {
  sds <- apply (d@mat, 1, sd, na.rm=TRUE)
  d.subset <- row.subset.gct (d, sds > sd.threshold)
  write.gct (d.subset, paste (out.prefix, '-sdfilter.gct', sep=''), appenddim=FALSE)
}

