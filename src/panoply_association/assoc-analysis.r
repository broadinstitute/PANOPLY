
source ('config.r')
Source ('markersel-classify.r')
Source ('stats.r')


gct.file <- file.path (norm.dir, paste (master.prefix, '.gct', sep=''))


run.marker.selection <- function (input.gct.file, input.cls.file, prefix) {
  # runs marker selection on input data for given class vector in input.cls
  marker.selection.and.classification (input.gct.file, input.cls.file, paste (prefix, '-analysis', sep=''), 
                                       id.col=id.col, desc.col=desc.col, gsea=TRUE,
                                       id.to.gene.map=NULL,   # GeneSymbol already present in GCT v1.3 input
                                       duplicate.gene.policy=duplicate.gene.policy,
                                       impute.colmax=sample.na.max,
                                       official.genenames=file.path ('..', 'data', 'gene-symbol-map.csv'))
  
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
                                           id.col=id.col, desc.col=desc.col, gsea=TRUE,
                                           id.to.gene.map=NULL,   # GeneSymbol already present in GCT v1.3 input
                                           duplicate.gene.policy=duplicate.gene.policy,
                                           impute.colmax=sample.na.max,
                                           official.genenames=file.path ('..', 'data', 'gene-symbol-map.csv'))
    }
  }
}


# obtain cls's
if (! exists ("assoc.subgroups")) {
  # no class vectors specified -- use all cls files from norm.dir (these have between 2-5 classes)
  cls.files <- list.files (norm.dir, pattern=".*\\.cls", full.names=TRUE)
  if (length (cls.files) > 0) {
    # marker selection and classification for each file in cls.files
    for (cls in cls.files) {   
      out.prefix <- gsub ("\\.cls$", "", basename (cls))
      if (min (summary (factor (read.cls (cls)))) < 3)  {
        # too few members in some class(es)
        warning ( paste (cls, "has classes with < 3 members ... skipping") )
        next
      }
      run.marker.selection (gct.file, cls, out.prefix)
    } 
  } else {
    warning ("No cls files to run association analysis on ... ignoring")
  }
} else {
  # groups file specified
  # (format similar to expt-design-file with Sample.ID and additional columns)
  # association analysis will be run for each additional column, excluding samples marked 'ignore'
  # (different columns cannot have the same subgroup name)
  subgroup.table <- read.csv (assoc.subgroups)
  rownames (subgroup.table) <- subgroup.table [, 'Sample.ID']
  cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
  
  ds <- parse.gctx (gct.file)
  sample.order <- ds@cid
  
  if (length (cls.list) > 0) {
    for (g in cls.list) {
      group <- make.names (subgroup.table [ sample.order, g ])  # use make.names to convert text to proper class labels
      subsamp <- group != 'ignore'
      # write subsamples dataset and class labels
      f <- paste (type, '-', g, sep='')
      ds.g <- col.subset.gct (ds, subsamp)
      gct.g <-  sprintf ("%s.gct", f)
      cls.g <- sprintf ("%s.cls", f)
      write.gct (ds.g, gct.g)
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
