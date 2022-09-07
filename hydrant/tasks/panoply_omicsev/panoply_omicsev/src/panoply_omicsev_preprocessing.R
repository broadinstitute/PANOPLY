###############################################################################
# FUNCTION: preprocessing for OmicsEV on PANOPLY
# AUTHOR: Stephanie Vartany
###############################################################################



################################################################################
## handle command line arguments
cat("Extracting command line arguments....\n")

args = commandArgs(trailingOnly = T)

if (length(args) == 3) {
  
  # all inputs provided
  harmonize_tar_file <- args[1]
  yaml_file <- args[2]
  class_colname <- args[3]
  
} else {
  stop("Incorrect number of inputs")
}

cat('harmonize_tar_file:', harmonize_tar_file, '\n')
cat('yaml_file:', yaml_file, '\n')
cat('class_colname:', class_colname, '\n')

################################################################################

library(cmapR)
library(dplyr)
library(yaml)
library(utils)

################################################################################
## untar harmonize inputs and extract files

ex_dir <- 'panoply_harmonize_output'
untar(harmonize_tar_file, exdir = ex_dir)

# get name of base directory that all files from tar are located in
base_dir <- list.files(ex_dir)
if (length(base_dir) != 1) {
  stop("Cannot find base directory in tar. Please provide tar file input from panoply_harmonize.")
}

# check that harmonize subdirectory is present
harmonize_dir <- 'harmonized-data'
if (!(harmonize_dir %in% list.files(file.path(ex_dir, base_dir)))) {
  stop("Cannot find harmonize directory in tar file. Please provide tar file input from panoply_harmonize.")
}

# check that correct matrix files exist from harmonize output
rna_file <- file.path(ex_dir, base_dir, harmonize_dir, 'rna-matrix.csv')
if (!file.exists(rna_file)) {
  stop("RNA file cannot be found. Please provide tar file input from panoply_harmonize.")
}
sample_anno_file <- file.path(ex_dir, base_dir, harmonize_dir, 'sample-info.csv')
if (!file.exists(sample_anno_file)) {
  stop("Sample annotation file cannot be found. Please provide tar file input from panoply_harmonize.")
}

# find protein file
data_file <- setdiff(list.files(file.path(ex_dir, base_dir, harmonize_dir), 
                                pattern = '*-matrix.csv'),
                     c('cna-matrix.csv', 'rna-matrix.csv'))
if (length(data_file) != 1) {
  stop("Protein data file cannot be found. Please provide tar file input from panoply_harmonize.")
}
data_file <- file.path(ex_dir, base_dir, harmonize_dir, data_file)




################################################################################


## define directory for dataset
data_dir <- "dataset/"
dir.create(data_dir)


## extract from yaml file
yaml_out <- read_yaml(yaml_file)
gene.id.col <- yaml_out$global_parameters$gene_mapping$gene_id_col


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
if (any(data_table < 0, na.rm=T)) { #log2-transform was performed
  data_table[,-1] <- 2^data_table[,-1] # undo log2-transform
}


## read and process sample annotations
sample_anno <- read.csv(sample_anno_file)

# find the experiment column
if (!('Experiment' %in% names(sample_anno))) {
  warning("Cannot find experiment column in sample annotation file. Defaulting to all same experiment")
  sample_anno$Experiment <- rep(1, dim(sample_anno)[1])
}

# check that experiment is all integers
if (!all(is.integer(sample_anno$Experiment))) {
  stop("Something is wrong with the 'Experiment' column. Make sure all values are integers.")
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
sample_anno_out <- sample_anno[, c('Sample.ID', class_colname, 'Experiment', 'order')]
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
if (any(rna_table < 0, na.rm=T)) { #log2-transform was performed
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
