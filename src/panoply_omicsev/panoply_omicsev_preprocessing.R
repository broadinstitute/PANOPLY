###############################################################################
# FUNCTION: preprocessing for OmicsEV on PANOPLY
# AUTHOR: Stephanie Vartany
###############################################################################



################################################################################
## handle command line arguments
print("extracting command line arguments")

args = commandArgs(trailingOnly = T)

if (length(args) == 4) {
  
  # all inputs provided
  data_files <- args[1]
  sample_anno_file <- args[2]
  
  if (args[3] == 'no_rna') {
    rna_file <- NULL
  } else {
    rna_file <- args[3]
  }
  
  class_colname <- args[4]
  
} else {
  stop("Incorrect number of inputs")
}

print(paste("data_files:", data_files, sep=' '))
print(paste("sample_anno_file:", sample_anno_file, sep=' '))
print(paste("rna_file:", rna_file, sep=' '))
print(paste("class_colname:", class_colname))

################################################################################

library(cmapR)
library(dplyr)

# define directory for datasets
data_dir <- "datasets/"
dir.create(data_dir)

data_files <- strsplit(data_files, ',')[[1]]

sample_names <- list()
column_descriptions <- list()

for (file in data_files) {

  ## read data
  data_gct <- parse_gctx(file)
  data_ids <- data.frame(ID = meta(data_gct, dimension='row')$geneSymbol)
  
  sample_names[[file]] <- colnames(mat(data_gct))
  column_descriptions[[file]] <- meta(data_gct, dimension='column')
  
  if (any(mat(data_gct) < 0, na.rm=T)) { #log2-transform was performed
    data_out <- cbind(data_ids, data.frame(2^mat(data_gct))) # undo log2-transform
  } else {
    data_out <- cbind(data_ids, data.frame(mat(data_gct)))
  }
  
  ## collapse to the gene level
  warning("collapsing protein data by geneSymbol")
  data_out <- aggregate(data_out[,-(1)], list(data_out$ID),
                                  function(x) mean(x, na.rm=T))
  names(data_out)[1] <- 'ID'
  
  # write data
  file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file))
  file_path <- paste(data_dir, file_basename, ".tsv", sep='')
  write.table(data_out,
              file = file_path,
              row.names=FALSE, sep="\t")
}

# # zip the proteome data into dataset directory
# zip("dataset.zip", files = file_path)

## read and process sample annotations
sample_anno <- read.csv(sample_anno_file)

## subset to only samples that are in the datasets files
all_samples <- unique(unlist(sample_names))
sample_anno <- filter(sample_anno, Sample.ID %in% all_samples)

# find the experiment column, make sure it is valid
if ('Experiment' %in% names(sample_anno)) {
  experiment <- sample_anno$Experiment
  names(experiment) <- sample_anno$Sample.ID
} else if (any(sapply(column_descriptions, function(x) 'Experiment' %in% names(x)))) {
  # see which omes have experiment cdesc
  ome_idx_with_experiment <- 
    which(sapply(column_descriptions, function(x) 'Experiment' %in% names(x)))
  
  # get all the samples with experiment info available
  all_samples_with_experiment <- unique(unlist(sample_names[ome_idx_with_experiment]))
  
  # find which -ome has the full experiment info
  ome_with_all_samples <- which(sapply(column_descriptions,
                                       function(x) setequal(x$Sample.ID, all_samples_with_experiment)))
  
  stopifnot("There is no data file with experiment info for all samples" = 
              length(ome_with_all_samples) > 0)
  
  # extract experiment
  experiment <- column_descriptions[[ome_with_all_samples[1]]]$Experiment
  names(experiment) <- column_descriptions[[ome_with_all_samples[1]]]$Sample.ID
  
  # make sure the experiment for all input data files are the same
  for (i in 1:length(column_descriptions)) {
    df1 <- column_descriptions[[1]]
    df2 <- column_descriptions[[i]]
    
    # subset to samples that overlap
    samples_intersect <- intersect(df1$Sample.ID, df2$Sample.ID)
    df1 <- filter(df1, Sample.ID %in% samples_intersect)
    df2 <- filter(df2, Sample.ID %in% samples_intersect)
    
    stopifnot('Experiment numbers do not correspond across different omes. Try runing each ome independently' = 
                setequal(df1$Experiment, df2$Experiment))
    
    # check that experiment names align with sample_anno file
    stopifnot(setequal(names(experiment), sample_anno$Sample.ID))
    
    # add to sample_anno
    sample_anno <- sample_anno[order(sample_anno$Sample.ID), ]
    experiment <- experiment[order(names(experiment))]
    sample_anno$Experiment <- experiment
  }
} else {
  stop("Cannot find experiment in sample annotation file or any of the dataset files")
}

# check that experiment is all integers
if (!all(is.integer(sample_anno$Experiment))) {
  stop("Something is wrong with the 'Experiment' column. Make sure all  values are integers.")
}

if ("Experiment" %in% names(sample_anno) & "Channel" %in% names(sample_anno)) {
  sample_anno$order <- order(sample_anno$Experiment, sample_anno$Channel)
} else {
  warning("Not reordering data based on experiment/channel. This should only matter for how the data is displayed")
  sample_anno$order <- 1:nrow(sample_anno)
}
sample_anno_out <- sample_anno[, c('Sample.ID', class_colname, 'Experiment', 'order')]
names(sample_anno_out) <- c('sample', 'class', 'batch', 'order')

# write sample annotations to tsv
write.table(sample_anno_out,
            file = "sample_list.tsv",
            row.names=FALSE, sep="\t")

## read and process rna file (x2)
if (!is.null(rna_file)) {
  rna_gct <- parse_gctx(rna_file)
  rna_ids <- data.frame(ID = meta(rna_gct, dimension='row')$geneSymbol)
  if (any(mat(rna_gct) < 0, na.rm=T)) { #log2-transform was performed
    rna_out <- cbind(rna_ids, data.frame(2^mat(rna_gct))) # undo log2-transform
  } else {
    rna_out <- cbind(rna_ids, data.frame(mat(rna_gct)))
  }
  
  # write rna to tsv
  write.table(rna_out,
              file = "x2.tsv",
              row.names=FALSE, sep="\t")
}
