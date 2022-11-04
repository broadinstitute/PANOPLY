################################################################################
# FUNCTION: preprocessing for OmicsEV on PANOPLY
# AUTHOR: Stephanie Vartany
################################################################################

library(optparse)
library(cmapR)
library(dplyr)
library(yaml)
library(utils)

################################################################################
## handle command line arguments
cat("\n\nExtracting command line arguments...\n")

parser <- OptionParser()
parser <- add_option(parser, c("--data_files"), type = 'character', dest = 'data_files')
parser <- add_option(parser, c("--rna_file"), type = 'character', dest = 'rna_file')
parser <- add_option(parser, c("--sample_anno_file"), type = 'character', dest = 'sample_anno_file')
parser <- add_option(parser, c("--yaml_file"), type = 'character', dest = "yaml_file")
parser <- add_option(parser, c("--STANDALONE"), type = 'logical', dest = "STANDALONE")
options <- parse_args(parser)

data_files <- options$data_files
rna_file <- options$rna_file
sample_anno_file <- options$sample_anno_file
yaml_file <- options$yaml_file
STANDALONE <- options$STANDALONE

## extract from yaml file
yaml_out <- read_yaml(yaml_file)
gene.id.col <- yaml_out$global_parameters$gene_mapping$gene_id_col
class_colname <- yaml_out$panoply_omicsev$class_column_name
batch_colname <- yaml_out$panoply_omicsev$batch_column_name
data_log_transformed <- yaml_out$panoply_omicsev$data_log_transformed
rna_log_transformed <- yaml_out$panoply_omicsev$rna_log_transformed
do_function_prediction <- yaml_out$panoply_omicsev$do_function_prediction

## make sure required inputs are present
if (is.null(yaml_file) | is.null(data_files) | is.null(sample_anno_file) | is.null(STANDALONE)) {
  stop("Required input missing")
}


cat("STANDALONE:", STANDALONE, '\n')
cat("data_files:", data_files, '\n')
cat("sample_anno_file:", sample_anno_file, '\n')
cat("rna_file:", rna_file, '\n')
cat("yaml_file:", yaml_file, '\n')
cat("class_colname:", class_colname, '\n')
cat('batch_colname:', batch_colname, '\n')
cat('data_log_transformed:', data_log_transformed, '\n')
cat('rna_log_transformed:', rna_log_transformed, '\n')
cat('gene.id.col:', gene.id.col, '\n\n')

################################################################################
## preprocessing functions

# for standalone inputs
preprocessing_STANDALONE <- function(data_files,
                                     sample_anno_file,
                                     rna_file,
                                     class_colname,
                                     batch_colname,
                                     data_log_transformed,
                                     rna_log_transformed,
                                     gene.id.col) {
  
  # define directory for dataset
  data_dir <- "dataset/"
  dir.create(data_dir)
  
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
}

# for harmonized inputs
preprocessing_harmonized <- function(data_file,
                                     sample_anno_file,
                                     rna_file,
                                     class_colname,
                                     batch_colname,
                                     data_log_transformed,
                                     rna_log_transformed,
                                     gene.id.col) {
  
  ## define directory for dataset
  data_dir <- "dataset/"
  dir.create(data_dir)
  
  ## load protein data file
  data_table <- read.csv(data_file)
  
  # change geneSymbol column name to ID
  stopifnot("First column of protein data not geneSymbol" = 
              names(data_table)[1] %in% c('geneSymbol', gene.id.col))
  names(data_table)[1] <- 'ID'
  
  # check that data is collapsed to gene level
  geneSymbols <- data_table$ID
  if (length(geneSymbols) != length(unique(geneSymbols))) {
    stop("There are gene symbol duplicates. PANOPLY harmonize not run correctly")
  }
  
  # undo log transformation 
  if (any(data_table < 0, na.rm=T) | data_log_transformed) { #log2-transform was performed
    data_table[,-1] <- 2^data_table[,-1] # undo log2-transform
  }
  
  
  ## read and process sample annotations
  sample_anno <- read.csv(sample_anno_file)
  
  # find the batch column
  if (!(batch_colname %in% names(sample_anno))) {
    warning("Cannot find batch column in sample annotation file. Defaulting to all same batch.")
    sample_anno$batch <- rep(1L, dim(sample_anno)[1])
  } else {
    sample_anno$batch <- sample_anno[, batch_colname]
  }
  
  # check that batch is all integers
  if (!all(is.integer(sample_anno$batch))) {
    stop("Something is wrong with the batch column. Make sure all values are integers.")
  }
  
  if ("Experiment" %in% names(sample_anno) & "Channel" %in% names(sample_anno)) {
    sample_anno$order <- order(sample_anno$Experiment, sample_anno$Channel)
  } else {
    warning("Not reordering data based on experiment/channel. This should only matter for how the data is displayed")
    sample_anno$order <- 1:nrow(sample_anno)
  }
  
  if (!(class_colname %in% names(sample_anno))) {
    stop("Class column name is not in the sample info csv")
  }
  sample_anno_out <- sample_anno[, c('Sample.ID', class_colname, 'batch', 'order')]
  names(sample_anno_out) <- c('sample', 'class', 'batch', 'order')
  
  
  ## read and process rna file (x2)
  rna_table <- read.csv(rna_file)
  
  # change geneSymbol column name to ID
  stopifnot("First column of rna data not geneSymbol" = 
              names(rna_table)[1] %in% c('geneSymbol', gene.id.col))
  names(rna_table)[1] <- 'ID'
  
  # check that rna is collapsed to gene level
  geneSymbols <- rna_table$ID
  if (length(geneSymbols) != length(unique(geneSymbols))) {
    stop("There are gene symbol duplicates in RNA. PANOPLY harmonize not run correctly")
  }
  
  # undo log transformation 
  if (any(rna_table < 0, na.rm=T) | rna_log_transformed) { #log2-transform was performed
    rna_table[,-1] <- 2^rna_table[,-1] # undo log2-transform
  }
  
  
  ## use only samples present in all files
  data_samples <- setdiff(names(data_table), 'ID')
  rna_samples <- setdiff(names(rna_table), 'ID')
  anno_samples <- sample_anno_out$sample
  intersect_samples <- intersect(intersect(data_samples, rna_samples), anno_samples)
  
  if (length(intersect_samples) == 0) {
    stop("No samples present in all three inputs (data, rna, sample annotation")
  } else if (length(intersect_samples) < max(sapply(list(data_samples, rna_samples, anno_samples), length))) {
    warning('not all samples present in data, rna, and sample annotation file)')
  }
  
  data_table <- data_table[, names(data_table) %in% c('ID', intersect_samples)]
  rna_table <- rna_table[, names(rna_table) %in% c('ID', intersect_samples)]
  sample_anno_out <- sample_anno_out[sample_anno_out$sample %in% intersect_samples, ]
  
  
  ## write output files
  
  # write data
  file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(data_file))
  file_path <- paste(data_dir, file_basename, ".tsv", sep='')
  write.table(data_table,
              file = file_path,
              row.names=FALSE, sep="\t")
  
  # write sample annotations to tsv
  write.table(sample_anno_out,
              file = "sample_list.tsv",
              row.names=FALSE, sep="\t")
  
  # write rna
  write.table(rna_table,
              file = 'x2.tsv',
              row.names=FALSE, sep="\t")
}

################################################################################
# main

# write do_function_pred to .txt file so it can be read later
write(do_function_prediction, "do_function_prediction.txt")

if (STANDALONE) {
  preprocessing_STANDALONE(data_files = data_files,
                           sample_anno_file = sample_anno_file,
                           rna_file = rna_file,
                           class_colname = class_colname,
                           batch_colname = batch_colname,
                           data_log_transformed = data_log_transformed,
                           rna_log_transformed = rna_log_transformed,
                           gene.id.col = gene.id.col)
} else {
  preprocessing_harmonized(data_file = data_files,
                           sample_anno_file = sample_anno_file,
                           rna_file = rna_file,
                           class_colname = class_colname,
                           batch_colname = batch_colname,
                           rna_log_transformed = rna_log_transformed,
                           data_log_transformed = data_log_transformed,
                           gene.id.col = gene.id.col)
}