# preprocess PANOPLY gct files for input to COSMO

###############################################################################
# extract command line inputs

args = commandArgs(trailingOnly = T)

if (length(args) == 3) {

  harmonize_tar_file <- args[1]
  yaml_file <- args[2]
  sample_label <- args[3]
  
  if (sample_label == "none") {
    sample_label <- NULL
  }
  
} else {
  stop("Incorrect number of inputs")
}

cat('harmonize_tar_file:', harmonize_tar_file, '\n')
cat('yaml_file:', yaml_file, '\n')
cat('sample_label:', sample_label, '\n\n\n')


###############################################################################
library(dplyr)
library(yaml)
library(cmapR)
library(utils)

###############################################################################
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



###############################################################################

data_dir <- "cosmo_preprocessed_data/"
dir.create(data_dir)

## extract from yaml file
yaml_out <- read_yaml(yaml_file)
gene.id.col <- yaml_out$global_parameters$gene_mapping$gene_id_col


## preprocess data file
data_table <- read.csv(data_file)

# change geneSymbol column name to ID
stopifnot("First column of data not geneSymbol" = 
            names(data_table)[1] %in% c('geneSymbol', gene.id.col))
names(data_table)[1] <- 'ID'

# check that data is collapsed to gene level
geneSymbols <- data_table$ID
if (length(geneSymbols) != length(unique(geneSymbols))) {
  stop("There are gene symbol duplicates. PANOPLY harmonize not run correctly")
}

# convert gene symbols to row names
rownames(data_table) <- data_table[,1]
data_table <- data_table[,-1]


## preprocess rna file
rna_table <- read.csv(rna_file)

# change geneSymbol column name to ID
stopifnot("First column of rna not geneSymbol" = 
            names(rna_table)[1] %in% c('geneSymbol', gene.id.col))
names(rna_table)[1] <- 'ID'

# check that rna is collapsed to gene level
geneSymbols <- rna_table$ID
if (length(geneSymbols) != length(unique(geneSymbols))) {
  stop("There are gene symbol duplicates. PANOPLY harmonize not run correctly")
}

# convert gene symbols to row names
rownames(rna_table) <- rna_table[,1]
rna_table <- rna_table[,-1]


## read and process sample annotations
sample_anno <- read.csv(sample_anno_file)

## use only samples present in all files
data_samples <- names(data_table)
rna_samples <- names(rna_table)
anno_samples <- sample_anno$Sample.ID
intersect_samples <- intersect(intersect(data_samples, rna_samples), anno_samples)

if (length(intersect_samples) == 0) {
  stop("No samples present in all three inputs (data, rna, sample annotation")
} else if (length(intersect_samples) < max(sapply(list(data_samples, rna_samples, anno_samples), length))) {
  warning('Not all samples are present in data, rna, and sample annotation file. Subsetting to samples present in all three files.')
}

data_table <- data_table[, names(data_table) %in% intersect_samples]
rna_table <- rna_table[, names(rna_table) %in% intersect_samples]
sample_anno <- sample_anno[sample_anno$Sample.ID %in% intersect_samples, ]


## check sample labels, extract from sample annotation file
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
  cat("Sample labels used:", good_sample_labels, '\n')
}

sample_anno_out <- sample_anno[, c('Sample.ID', good_sample_labels)]
names(sample_anno_out)[1] <- 'sample'



# save good sample labels
writeLines(paste(good_sample_labels, collapse=','), 
           con = paste(data_dir, "sample_label.txt", sep=''))



## write to tsv files

# write data
write.table(data_table,
            file = paste(data_dir, 'data.tsv', sep=''),
            row.names=TRUE, sep="\t")

# write rna
write.table(rna_table,
            file = paste(data_dir, 'rna.tsv', sep=''),
            row.names=TRUE, sep="\t")

# write sample annotations to tsv
write.table(sample_anno_out,
            file = paste(data_dir, 'sample-list.tsv', sep=''),
            row.names=FALSE, sep="\t")


