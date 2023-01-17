#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')
library (pamr)


data <- parse.gctx (file.path (filt.dir, sprintf ("%s-ratio-norm-NArm.gct", type)))

# if data has missing values, remove items with > na.threshold NAs and impute
if (sum (is.na(data@mat)) > 0) {
  nas <- apply (data@mat, 1, function (x) sum (is.na(x)))
  keep <- nas < (ncol(data@mat) * clustering.na.threshold)
  data.NAthresh <- row.subset.gct (data, keep)
  data.imputed <- pamr.knnimpute (data=list (x=as.matrix(data.NAthresh@mat)), colmax=sample.na.max)
  data.NAthresh@mat <- data.imputed$x
  data <- data.NAthresh
}

# eliminate features with not enough variation
feature.sd <- apply (data@mat, 1, sd, na.rm=TRUE)
keep <- feature.sd > clustering.sd.threshold
data <- row.subset.gct (data, keep)

cluster.data <- data.frame (Name=data@rid, Description=data@rid,  data@mat)
write.gct2 (cluster.data, sprintf ("%s-data.gct", type))

