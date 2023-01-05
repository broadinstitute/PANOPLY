# preprocess PANOPLY files for input to COSMO

###############################################################################
# extract command line inputs

args = commandArgs(trailingOnly = T)

if (length(args) == 4) {
  d1_path <- args[1]
  d2_path <- args[2]
  sample_csv_path <- args[3]
  yaml_file <- args[4]
  
}

cat("d1_path:", d1_path, '\n')
cat("d2_path:", d2_path, '\n')
cat("sample_csv_path:", sample_csv_path, '\n')
cat("yaml_file:", yaml_file, '\n\n\n')


###############################################################################
library(dplyr)
library(yaml)
library(cmapR)

data_dir <- "cosmo_preprocessed_data/"
dir.create(data_dir)

## extract from yaml file
yaml_out <- read_yaml(yaml_file)
gene.id.col <- yaml_out$global_parameters$gene_mapping$gene_id_col
sample_label <- yaml_out$cosmo.params$sample_label

## read data from d1 file and d2 file
for (file in c(d1_path, d2_path)) {
  ## read data
  ex <- strsplit(basename(file), split="\\.")[[1]][-1]
  if (ex == 'csv') {
    data_out <- read.csv(file)
    # change geneSymbol column name to ID
    stopifnot("First column of data not geneSymbol" = names(data_out)[1] %in% c('geneSymbol', gene.id.col))
    names(data_out)[1] <- 'ID'
  } else if (ex == 'gct') {
    data_gct <- parse_gctx(file)
    data_ids <- data.frame(ID = meta(data_gct, dimension='row')[, gene.id.col])
    data_out <- cbind(data_ids, data.frame(mat(data_gct)))
  } else {
    stop("Invalid file input")
  }
  
  
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
}

## read and process sample annotations
sample_anno <- read.csv(sample_csv_path)
sample_label_list <- unlist(strsplit(x=sample_label, split=','))

good_sample_labels <- c()
for (label in sample_label_list) {
  if (!(label %in% names(sample_anno))) {
    warning(paste("Sample label", label, "is not in sample annotation file"))
  } else if (min(base::table(sample_anno[, label])) < min(10, dim(sample_anno)[1] / 3)) {
    warning(paste("Sample label", label, "is not well-balanced. It is being excluded."))
  } else if (any(is.na(sample_anno[, label]))) {
    warning(paste("Sample label", label, "has NAs. It is being excluded."))
  } else if (length(unique(sample_anno[, label])) < 2) {
    warning("Sample label", label, "only has one level. It is being excluded.")
  } else if (length(unique(sample_anno[, label])) > 2) {
    warning(paste("Sample label", label, "being excluded because COSMO does not handle classes with more than 2 levels."))
  } else {
    good_sample_labels <- c(good_sample_labels, label)
  }
}
if (length(good_sample_labels) == 0) {
  stop("No valid sample annotation columns were found")
} else {
  cat("Sample labels used:", good_sample_labels)
}

sample_anno_out <- sample_anno[, c('Sample.ID', good_sample_labels)]
names(sample_anno_out)[1] <- 'sample'

# write sample annotations to tsv
file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(sample_csv_path))
file_path <- paste(data_dir, file_basename, ".tsv", sep='')
write.table(sample_anno_out,
            file = file_path,
            row.names=FALSE, sep="\t")

# save good sample labels
writeLines(paste(good_sample_labels, collapse=','), 
           con = paste(data_dir, "sample_label.txt", sep=''))
