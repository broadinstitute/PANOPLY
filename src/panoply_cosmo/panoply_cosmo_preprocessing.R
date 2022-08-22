# preprocess PANOPLY gct files for input to COSMO

###############################################################################
# extract command line inputs

args = commandArgs(trailingOnly = T)

if (length(args) == 4) {
  d1_gct_path <- args[1]
  d2_gct_path <- args[2]
  sample_csv_path <- args[3]
  sample_label <- args[4]
}

print(paste("d1_gct_path:", d1_gct_path))
print(paste("d2_gct_path:", d2_gct_path))
print(paste("sample_csv_path:", sample_csv_path))
print(paste("sample_label:", sample_label))


###############################################################################
library(cmapR)
library(dplyr)

data_dir <- "cosmo_preprocessed_data/"
dir.create(data_dir)

sample_names <- list()
for (file in c(d1_gct_path, d2_gct_path)) {
  ## read data
  data_gct <- parse_gctx(file)
  data_ids <- data.frame(ID = meta(data_gct, dimension='row')$geneSymbol)
  data_out <- cbind(data_ids, data.frame(mat(data_gct)))
  
  ## collapse to the gene level
  warning("collapsing protein data by geneSymbol")
  data_out <- aggregate(data_out[,-(1)], list(data_out$ID),
                        function(x) mean(x, na.rm=T))
  
  # convert gene symbols to row names
  rownames(data_out) <- data_out[,1]
  data_out <- data_out[,-1]
  
  # write data
  file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file))
  file_path <- paste(data_dir, file_basename, ".tsv", sep='')
  write.table(data_out,
              file = file_path,
              row.names=TRUE, sep="\t")
  
  sample_names[[file]] <- names(data_out)
}

## read and process sample annotations
sample_anno <- read.csv(sample_csv_path)
sample_label_list <- unlist(strsplit(x=sample_label, split=','))
if (length(setdiff(sample_label_list, names(sample_anno))) > 0) {
  stop("At least one sample label is not a column in sample annotation file")
}
# check if any columns have only one observation in a level
for (label in sample_label_list) {
  if (any(base::table(sample_anno[, label]) == 1)) {
    stop(paste("Sample label '", label, "' has level(s) '",
                  names(which(base::table(sample_anno[,label]) == 1)),
                  "' with only one sample observation. ",
                  "This is not allowed in the cosmo function. ",
                  "Remove this label from the list or edit the data.",
                  sep=''))
  }
}
sample_anno_out <- sample_anno[, c('Sample.ID', sample_label_list)]
names(sample_anno_out)[1] <- 'sample'

# write sample annotations to tsv
file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(sample_csv_path))
file_path <- paste(data_dir, file_basename, ".tsv", sep='')
write.table(sample_anno_out,
            file = file_path,
            row.names=FALSE, sep="\t")