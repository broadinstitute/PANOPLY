###############################################################################
# FUNCTION: preprocessing for OmicsEV on PANOPLY
# AUTHOR: Stephanie Vartany
###############################################################################



################################################################################
## handle command line arguments
cat("\n\nExtracting command line arguments...\n")

args = commandArgs(trailingOnly = T)

if (length(args) == 8) {
  
  # all inputs provided
  data_files <- args[1]
  sample_anno_file <- args[2]
  
  if (args[3] == 'no_rna') {
    rna_file <- NULL
  } else {
    rna_file <- args[3]
  }
  
  yaml_file <- args[4]
  class_colname <- args[5]
  batch_colname <- args[6]
  data_log_transformed <- as.logical(args[7])
  rna_log_transformed <- as.logical(args[8])
  
  stopifnot(!is.na(data_log_transformed))
  stopifnot(!is.na(rna_log_transformed))
  
} else {
  stop("Incorrect number of inputs")
}

cat("data_files:", data_files, '\n')
cat("sample_anno_file:", sample_anno_file, '\n')
cat("rna_file:", rna_file, '\n')
cat("yaml_file:", yaml_file, '\n')
cat("class_colname:", class_colname, '\n')
cat('batch_colname:', batch_colname, '\n')
cat('data_log_transformed:', data_log_transformed, '\n')
cat('rna_log_transformed:', rna_log_transformed, '\n\n')

################################################################################

library(cmapR)
library(dplyr)
library(yaml)


# define directory for dataset
data_dir <- "dataset/"
dir.create(data_dir)

## extract from yaml file
yaml_out <- read_yaml(yaml_file)
gene.id.col <- yaml_out$global_parameters$gene_mapping$gene_id_col

data_files <- strsplit(data_files, ',')[[1]]

sample_names <- list()
column_descriptions <- list()
data_outs <- list()

for (file in data_files) {

  ## read data
  data_gct <- parse_gctx(file)
  data_ids <- data.frame(ID = meta(data_gct, dimension='row')[,gene.id.col])
  
  sample_names[[file]] <- colnames(mat(data_gct))
  column_descriptions[[file]] <- meta(data_gct, dimension='column')
  
  if (any(mat(data_gct) < 0, na.rm=T) | data_log_transformed) { #log2-transform was performed
    data_out <- cbind(data_ids, data.frame(2^mat(data_gct))) # undo log2-transform
  } else {
    data_out <- cbind(data_ids, data.frame(mat(data_gct)))
  }
  
  ## collapse to the gene level
  warning("collapsing protein data by geneSymbol")
  data_out <- aggregate(data_out[,-(1)], list(data_out$ID),
                                  function(x) mean(x, na.rm=T))
  names(data_out)[1] <- 'ID'
  
  # save data in list
  data_outs[[file]] <- data_out
}


## read and process sample annotations
sample_anno <- read.csv(sample_anno_file)

## subset to only samples that are in the datasets files
all_samples <- unique(unlist(sample_names))
sample_anno <- filter(sample_anno, Sample.ID %in% all_samples)

# find the batch column, make sure it is valid
if (batch_colname %in% names(sample_anno)) {
  batch <- sample_anno[, batch_colname]
  names(batch) <- sample_anno$Sample.ID
  sample_anno$batch <- batch
} else if (any(sapply(column_descriptions, function(x) batch_colname %in% names(x)))) {
  # see which omes have batch cdesc
  ome_idx_with_batch <- 
    which(sapply(column_descriptions, function(x) batch_colname %in% names(x)))
  
  # get all the samples with batch info available
  all_samples_with_batch <- unique(unlist(sample_names[ome_idx_with_batch]))
  
  # find which -ome has the full batch info
  ome_with_all_samples <- which(sapply(column_descriptions,
                                       function(x) setequal(x$Sample.ID, all_samples_with_batch)))
  
  stopifnot("There is no data file with batch info for all samples" = 
              length(ome_with_all_samples) > 0)
  
  # extract batch
  batch <- column_descriptions[[ome_with_all_samples[1]]][, batch_colname]
  names(batch) <- column_descriptions[[ome_with_all_samples[1]]]$Sample.ID
  
  # make sure the batch for all input data files are the same
  for (i in 1:length(column_descriptions)) {
    df1 <- column_descriptions[[1]]
    df2 <- column_descriptions[[i]]
    
    # subset to samples that overlap
    samples_intersect <- intersect(df1$Sample.ID, df2$Sample.ID)
    df1 <- filter(df1, Sample.ID %in% samples_intersect)
    df2 <- filter(df2, Sample.ID %in% samples_intersect)
    
    stopifnot('batch numbers do not correspond across different omes. Try runing each ome independently' = 
                setequal(df1[, batch_colname], df2[, batch_colname]))
    
  }
  
  # add batch to sample_anno
  sample_anno <- sample_anno[sample_anno$Sample.ID %in% names(batch), ]
  sample_anno <- sample_anno[order(sample_anno$Sample.ID), ]
  batch <- batch[order(names(batch))]
  sample_anno$batch <- batch
  
} else {
  warning("Cannot find batch in sample annotation file or any of the dataset files. Defaulting to all samples in same batch")
  sample_anno$batch <- rep(1L, dim(sample_anno)[1])
  batch <- sample_anno$batch
}

# check that batch is all integers
if (!all(is.integer(sample_anno$batch))) {
  stop("Something is wrong with the batch column. Make sure all values are integers.")
}

if ("Experiment" %in% names(sample_anno) & "Channel" %in% names(sample_anno)) {
  sample_anno$order <- order(sample_anno$Experiment, sample_anno$Channel)
} else {
  warning("Not reordering data based on experiment/channel. This should only matter for how the data is displayed.")
  sample_anno$order <- 1:nrow(sample_anno)
}
sample_anno_out <- sample_anno[, c('Sample.ID', class_colname, 'batch', 'order')]
names(sample_anno_out) <- c('sample', 'class', 'batch', 'order')


## read and process rna file (x2)
if (!is.null(rna_file)) {
  rna_gct <- parse_gctx(rna_file)
  rna_ids <- data.frame(ID = meta(rna_gct, dimension='row')[,gene.id.col])
  if (any(mat(rna_gct) < 0, na.rm=T) | rna_log_transformed) { #log2-transform was performed
    rna_out <- cbind(rna_ids, data.frame(2^mat(rna_gct))) # undo log2-transform
  } else {
    rna_out <- cbind(rna_ids, data.frame(mat(rna_gct)))
  }
}


## subset to only samples present in all data files, write data
all_samples <- sample_anno_out$sample
for (file in data_files) {
  all_samples <- intersect(all_samples, sample_names[[file]])
}
if (!is.null(rna_file)) {
  all_samples <- intersect(all_samples, setdiff(names(rna_out), 'ID'))
}
if (length(all_samples) == 0) {
  stop("No samples found that are present in all inputted data files.")
}

# write data files
for (file in data_files) {
  file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file))
  file_path <- paste(data_dir, file_basename, ".tsv", sep='')
  data_out <- data_outs[[file]][, c('ID', all_samples)]
  write.table(data_out,
              file = file_path,
              row.names=FALSE, sep="\t")
}

# write sample annotations to tsv
sample_anno_out <- sample_anno_out[sample_anno_out$sample %in% all_samples, ]
write.table(sample_anno_out,
            file = "sample_list.tsv",
            row.names=FALSE, sep="\t")

# write rna
if (!is.null(rna_file)) {
  rna_out <- rna_out[, c('ID', all_samples)]
  write.table(rna_out,
              file = 'x2.tsv',
              row.names=FALSE, sep="\t")
}

