#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

# check to make sure libraries needed by cna-analysis.r are already installed
# (else, many jobs can fail before one of the jobs ends up installing the libraries)
if (!require("pacman")) install.packages ("pacman")
pacman::p_load (psych)


source ('config.r')
min.cna.N <- 5

# parameters
if (! exists ("pe.max")) pe.max <- pe.max.default
if (! exists ("cna.subgroups")) cna.subgroups <- NULL



create.subset <- function (data, sub, fname, na.max.subset=0.5) {
  # create data subset; ensure there aren't too many NAs after subsetting
  data.sub <- data [, c (TRUE, sub)]
  nas <- apply (data.sub[-1], 1, function (x) sum (is.na (x))) / (ncol (data.sub) - 1)
  data.sub <- data.sub [nas < na.max.subset, ]
  write.csv (data.sub, fname, row.names=FALSE)
  return (fname)
}


setup.corr.calcs <- function (groups=NULL, rna.file=file.path (harmonize.dir, 'rna-matrix.csv'),
                              cna.file=file.path (harmonize.dir, 'cna-matrix.csv'),
                              pome.file=file.path (harmonize.dir, paste (type, '-matrix.csv', sep='')),
                              matrix.files.table="file_table.tsv", subgroups.file="subgroups.txt",
                              force=FALSE) {
  
  rna.data <- read.csv (rna.file)
  cna.data <- read.csv (cna.file)
  pome.data <- read.csv (pome.file)
  
  # defaults for cls and groups
  if (is.null (groups)) {
    subgroup.table <- data.frame (all=rep ('all', ncol(pome.data)-1)) # by default, look at all samples together
    cls.list <- c ('all')
  } else {
    # groups file specified
    # (format similar to expt-design-file with Sample.ID and additional columns)
    # CNA analysis will be run for each additional column, excluding samples marked 'ignore'
    # (different columns cannot have the same subgroup name)
    subgroup.table <- read.csv (groups)
    rownames (subgroup.table) <- subgroup.table [, 'Sample.ID']
    subgroup.table <- subgroup.table [colnames (pome.data)[-1], ]  # reorder to match harmonzied data files
    cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
  }
  
  # create data subsets and populate matrix files tables
  mfile <- file (matrix.files.table, "w")
  sfile <- file (subgroups.file, "w")
  for (i in 1:length (cls.list)) {
    cls <- as.character (subgroup.table [, cls.list[i]])
    groups <- setdiff (unique (cls), 'ignore')
    for (x in groups) {
      subsamp <- cls == x
      subsamp[is.na(subsamp)] <- FALSE # overwrite NA with FALSE
      if (sum (subsamp, na.rm=TRUE) >= min.cna.N) {
        # run only if there are enough samples
        create.subset (rna.data, subsamp, paste (x, '-rna-matrix.csv', sep=''))
        create.subset (cna.data, subsamp, paste (x, '-cna-matrix.csv', sep=''))
        create.subset (pome.data, subsamp, paste (x, '-pome-matrix.csv', sep=''))
        cat (x, file=sfile, sep='\n')
        cat (paste (x, c ('-rna-matrix.csv', '-cna-matrix.csv', '-pome-matrix.csv'), sep=''), '\n',
             file=mfile, sep='\t')
      } else {
        cat ('Not enough samples for CNA analysis: Skipping', x, '... \n')
      }
    }
  }
  close (mfile)
  close (sfile)
  
}


setup.corr.calcs (groups=cna.subgroups)
