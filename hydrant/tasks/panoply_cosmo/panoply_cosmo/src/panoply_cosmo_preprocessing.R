# preprocess PANOPLY gct files for input to COSMO

###############################################################################
# extract command line inputs

args = commandArgs(trailingOnly = T)

if (length(args) == 5) {
  d1_gct_path <- args[1]
  d2_gct_path <- args[2]
  sample_csv_path <- args[3]
  sample_label <- args[4]
  yaml_file <- args[5]
  
  if (sample_label == "none") {
    sample_label <- NULL
  }
}

print(paste("d1_gct_path:", d1_gct_path))
print(paste("d2_gct_path:", d2_gct_path))
print(paste("sample_csv_path:", sample_csv_path))
print(paste("sample_label:", sample_label))
print(paste("yaml_file:", yaml_file))


###############################################################################
library(dplyr)
library(yaml)
library(cmapR)

data_dir <- "cosmo_preprocessed_data/"
dir.create(data_dir)

## extract from yaml file
yaml_out <- read_yaml(yaml_file)
gene.id.col <- yaml_out$global_parameters$gene_mapping$gene_id_col

## read data from d1 file and d2 file
for (file in c(d1_gct_path, d2_gct_path)) {
  ## read data
  data_gct <- parse_gctx(file)
  data_ids <- data.frame(ID = meta(data_gct, dimension='row')[, gene.id.col])
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
}

## read and process sample annotations
sample_anno <- read.csv(sample_csv_path)

if (is.null(sample_label)) {
  warning("using groups from yaml file as default sample attributes")
  sample_label_list <- yaml_out$groups.cols
} else {
  sample_label_list <- unlist(strsplit(x=sample_label, split=','))
}

good_sample_labels <- c()
for (label in sample_label_list) {
  if (!(label %in% names(sample_anno))) {
    warning(paste(label, "is not in sample annotation file"))
  } else if (any(base::table(sample_anno[, label]) == 1)) {
    warning(paste("Sample label '", label, "' has level(s) '",
               names(which(base::table(sample_anno[,label]) == 1)),
               "' with only one sample observation. ",
               "This is not allowed in the cosmo function. ",
               "It is currently being removed as an input.",
               sep=''))
  } else {
    good_sample_labels <- c(good_sample_labels, label)
  }
}
if (length(good_sample_labels) == 0) {
  stop("No valid sample annotation columns were found")
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
