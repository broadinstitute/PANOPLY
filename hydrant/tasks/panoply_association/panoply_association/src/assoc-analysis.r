#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')
Source ('markersel-classify.r')
Source ('stats.r')


gct.file <- file.path (filt.dir, paste (master.prefix, '.gct', sep=''))


run.marker.selection <- function (input.gct.file, input.cls.file, prefix, run.1vAll=FALSE) {
  # runs marker selection on input data for given class vector in input.cls
  tryCatch (marker.selection.and.classification (input.gct.file, input.cls.file, paste (prefix, '-analysis', sep=''), 
                                                 gene.id.col=gene.id.col, id.col=id.col, desc.col=desc.col, gsea=FALSE,
                                                 id.to.gene.map=NULL,   # GeneSymbol already present in GCT v1.3 input
                                                 duplicate.gene.policy=duplicate.gene.policy,
                                                 impute.colmax=sample.na.max,
                                                 official.genenames=file.path ('..', 'data', 'gene-symbol-map.csv'),
                                                 fdr=assoc.fdr,
                                                 models=c("pls","rf","glmnet")),
            error = function(cond) {
              message(paste("Failed to complete marker selection for ", prefix))
              message(cond)
            })
 
  if (run.1vAll) {
    # if class vector has > 2 classes, and sufficient numbers per class,
    # run 1 vs. all marker selection for each class
    cls <- read.cls (input.cls.file)
    if (nlevels (factor (cls)) > 2 && min (summary (factor (cls))) >= 10) {
      cls.1vAll <- classes.1vAll (cls)
      for (i in 1:ncol(cls.1vAll)) {
        prefix.1vA <- paste (prefix, colnames(cls.1vAll)[i], sep='-')
        new.clsf <- paste (prefix.1vA, '.cls', sep='')
        write.cls (cls.1vAll[,i], new.clsf)
        marker.selection.and.classification (input.gct.file, new.clsf, paste (prefix.1vA, '-analysis', sep=''),
                                             gene.id.col=gene.id.col, id.col=id.col, desc.col=desc.col, gsea=FALSE,
                                             id.to.gene.map=NULL,   # GeneSymbol already present in GCT v1.3 input
                                             duplicate.gene.policy=duplicate.gene.policy,
                                             impute.colmax=sample.na.max,
                                             official.genenames=file.path ('..', 'data', 'gene-symbol-map.csv'),
                                             fdr=assoc.fdr,
                                             models=c("pls","rf","glmnet"))
      }
    }
  }
}

# this module ASSUMES that a groups file (assoc.subgroups) is specified

if (! exists ("assoc.subgroups")) {
  stop("A groups file must be provided. This file should include Sample.ID,
       additional row-description columns to be analyzed for enrichent.")
} else {
  # (groups file format is similar to expt-design-file, with Sample.ID and additional columns)
  # association analysis will be run for each additional column, excluding samples marked 'ignore'
  # (different columns cannot have the same subgroup name)
  subgroup.table <- read.csv (assoc.subgroups)
  rownames (subgroup.table) <- subgroup.table [, 'Sample.ID']
  cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
  
  ds <- parse.gctx (gct.file)
  sample.order <- ds@cid
  
  if (length (cls.list) > 0) {
    for (g in cls.list) {
      #replace "" with NA, otherwise it gets kept
      group <- subgroup.table[sample.order, g]
      group[group==""]=NA
      group <- make.names (group)  # use make.names to convert text to proper class labels
      #REMOVE NAs!
      subsamp <- !group%in%c('ignore',"NA.")
      # write subsamples dataset and class labels
      f <- paste (type, '-', g, sep='')
      ds.g <- col.subset.gct (ds, subsamp)
      gct.g <-  sprintf ("%s.gct", f)
      cls.g <- sprintf ("%s.cls", f)
      write.gct (ds.g, gct.g, appenddim=FALSE)
      write.cls (group [subsamp], cls.g)
      
      if (min (summary (factor (read.cls (cls.g)))) < 3) {
        # too few members in some class(es)
        warning ( paste (g, "has classes with < 3 members ... skipping") )
        next
      }
      run.marker.selection (gct.g, cls.g, f)
    }
    
  } else {
    warning ("No cls files to run association analysis on ... ignoring")
  }
}
