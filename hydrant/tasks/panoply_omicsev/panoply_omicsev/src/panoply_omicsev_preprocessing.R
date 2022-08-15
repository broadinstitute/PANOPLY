###############################################################################
# FUNCTION: preprocessing for OmicsEV on PANOPLY
# AUTHOR: Stephanie Vartany
###############################################################################



################################################################################
## handle command line arguments
print("extracting command line arguments")

args = commandArgs(trailingOnly = T)

if (length(args) == 3) {
  
  # all inputs provided
  data_files <- args[1]
  sample_anno_file <- args[2]
  
  if (args[3] == 'no_rna') {
    rna_file <- NULL
  } else {
    rna_file <- args[3]
  }
  
} else {
  stop("Incorrect number of inputs")
}

print(paste("data_files:", data_files, sep=' '))
print(paste("sample_anno_file:", sample_anno_file, sep=' '))
print(paste("rna_file:", rna_file, sep=' '))

################################################################################

library(cmapR)
library(dplyr)

# define directory for datasets
data_dir <- "datasets/"
dir.create(data_dir)

data_files <- strsplit(data_files, ',')[[1]]

for (file in data_files) {

  ## read data
  data_gct <- parse_gctx(file)
  data_ids <- data.frame(ID = meta(data_gct, dimension='row')$geneSymbol)
  
  if (any(mat(data_gct) < 0)) { #log2-transform was performed
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
sample_anno$order <- order(sample_anno$Experiment, sample_anno$Channel)
sample_anno_out <- sample_anno[, c('Sample.ID', 'PAM50', 'Experiment', 'order')]
names(sample_anno_out) <- c('sample', 'class', 'batch', 'order')

# write sample annotations to tsv
write.table(sample_anno_out,
            file = "sample_list.tsv",
            row.names=FALSE, sep="\t")

## read and process rna file (x2)
if (!is.null(rna_file)) {
  rna_gct <- parse_gctx(rna_file)
  rna_ids <- data.frame(ID = meta(rna_gct, dimension='row')$geneSymbol)
  if (any(mat(rna_gct) < 0)) { #log2-transform was performed
    rna_out <- cbind(rna_ids, data.frame(2^mat(rna_gct))) # undo log2-transform
  }
  
  # write rna to tsv
  write.table(rna_out,
              file = "x2.tsv",
              row.names=FALSE, sep="\t")
}
