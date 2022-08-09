###############################################################################
# FUNCTION: preprocessing for OmicsEV on PANOPLY
# AUTHOR: Stephanie Vartany
###############################################################################



################################################################################
## handle command line arguments
print("extracting command line arguments")

args = commandArgs(trailingOnly = T)

if (length(args) == 2) {
  
  # all inputs provided
  data_file <- args[1]
  sample_anno_file <- args[2]
  rna_file <- NULL
  
} else if (length(args) == 3) {
  
  # rna correlation input not provided
  data_file <- args[1]
  sample_anno_file <- args[2]
  rna_file <- args[3]
} else {
  stop("Incorrect number of inputs")
}

print(paste("data_file:", data_file, sep=' '))
print(paste("sample_anno_file:", sample_anno_file, sep=' '))
print(paste("rna_file:", rna_file, sep=' '))

################################################################################

library(cmapR)
library(dplyr)

# other important file paths
data_dir <- 'dataset/'
sample_anno_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(sample_anno_file))
sample_anno_path <- paste(sample_anno_basename, ".tsv", sep='')

## read proteome data
data_gct <- parse_gctx(data_file)
data_ids <- data.frame(ID = meta(data_gct, dimension='row')$geneSymbol)
data_out <- cbind(data_ids, data.frame(2^mat(data_gct))) # undo log2-transform

## GET RID OF ARTIFACT IN PROTEOME FILE
warning("getting rid of weird artifact in proteome_aggregate.gct file!")
weird.gene.ids <- sort(unique(data_out$ID))[1:17]
data_out <- data_out[!(data_out$ID %in% weird.gene.ids), ]

## collapse to the gene level
if (!is.null(rna_file)) {
  data_out <- aggregate(data_out[,-(1)], list(data_out$ID),
                                  function(x) mean(x, na.rm=T))
  names(data_out)[1] <- 'ID'
}


# write proteome data to 'dataset' directory
dir.create(data_dir)
file_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(data_file))
write.table(data_out,
            file = paste(data_dir, file_basename, ".tsv", sep=''),
            row.names=FALSE, sep="\t")

## read and process sample annotations
sample_anno <- read.csv(sample_anno_file)
sample_anno$order <- order(sample_anno$Experiment, sample_anno$Channel)
sample_anno_out <- sample_anno[, c('Sample.ID', 'PAM50', 'Experiment', 'order')]
names(sample_anno_out) <- c('sample', 'class', 'batch', 'order')

# write sample annotations to tsv
write.table(sample_anno_out,
            file = sample_anno_path,
            row.names=FALSE, sep="\t")

## read and process rna file (x2)
if (!is.null(rna_file)) {
  # make path for RNA data
  rna_basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(rna_file))
  rna_path <- paste(rna_basename, ".tsv", sep='')
  
  rna_gct <- parse_gctx(rna_file)
  rna_ids <- data.frame(ID = meta(rna_gct, dimension='row')$geneSymbol)
  if (any(mat(rna_gct) < 0)) { #log2-transform was performed
    rna_out <- cbind(rna_ids, data.frame(2^mat(rna_gct))) # undo log2-transform
  }
  
  # write rna to tsv
  write.table(rna_out,
              file = rna_path,
              row.names=FALSE, sep="\t")
}







