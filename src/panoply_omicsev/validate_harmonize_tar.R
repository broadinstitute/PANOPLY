#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

################################################################################
# FUNCTION: validate harmonize tar file for omicsev
# AUTHOR: Stephanie Vartany
################################################################################

# extract command line inputs
tar_dir = commandArgs(trailingOnly=TRUE)[1]
ome_type = commandArgs(trailingOnly=TRUE)[2]

# remove the last '/' from tar_dir if it is there, messes up the file paths
last_char <- substr(tar_dir, nchar(tar_dir), nchar(tar_dir))
if (last_char == '/') {
  tar_dir <- substr(tar_dir, 1, nchar(tar_dir)-1)
}

# check that harmonize subdirectory is present
harmonize_dir <- 'harmonized-data'
if (!(harmonize_dir %in% list.files(tar_dir))) {
  stop("Cannot find harmonize directory in tar file. Please provide tar file input from panoply_harmonize.")
}

# check that correct matrix files exist from harmonize output
rna_file <- file.path(tar_dir, harmonize_dir, 'rna-matrix.csv')
if (!file.exists(rna_file)) {
  stop("RNA file cannot be found. Please provide tar file input from panoply_harmonize.")
}
sample_anno_file <- file.path(tar_dir, harmonize_dir, 'sample-info.csv')
if (!file.exists(sample_anno_file)) {
  stop("Sample annotation file cannot be found. Please provide tar file input from panoply_harmonize.")
}

# find protein file
data_file <- file.path(tar_dir, harmonize_dir, paste0(ome_type, '-matrix.csv'))
if (!file.exists(data_file)) {
  stop("Protein data file cannot be found. Please provide tar file input from panoply_harmonize.")
}
