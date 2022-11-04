################################################################################
# FUNCTION: validate harmonize tar file for omicsev
# AUTHOR: Stephanie Vartany
################################################################################

# extract command line inputs
ex_dir = commandArgs(trailingOnly=TRUE)[1]


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
write(data_file, "data_file.txt")
