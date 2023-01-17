#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
library(optparse)

options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-i", "--input"), action="store", dest='tar.file' , type="character" ,help="input tar-file or (extracted) folder."),
  make_option(c("-l", "--label"), action="store", dest='label', type="character", help=""),
  make_option(c("-t", "--type"), action="store", dest='type', type="character", help="Data type, i.e. proteome, phospho, rna, ..."),
  make_option(c("-d", "--tempdir"), action="store", dest='tmp.dir', type="character", help="temp-folder", default="tmp"),
  make_option(c("-f", "--filtdir"), action="store", dest='filt.dir', type="character", help="filtered-folder", default="filtered-data"),
  make_option(c("-c", "--clustdir"), action="store", dest='clust.dir', type="character", help="clustering-folder", default="clustering"),
  make_option(c("-u", "--k_min"), action="store", dest='k_min', type="numeric", help="Minimal cluster number.", default=2),
  make_option(c("-v", "--k_max"), action="store", dest='k_max', type="numeric", help="Maximal cluster number.", default=10),
  make_option(c("-s", "--sdfilter"), action="store", dest='clustering.sd.threshold', type="numeric", help="Minmal standard deviation acreoss samples required for clustering.", default=1.5),
  make_option(c("-e", "--min_features_sd"), action="store", dest='min.feat.sd', type="numeric", help="Minmal standard deviation across samples required for clustering.", default=500),
  make_option(c("-m", "--method"), action="store", dest='method', type="character", help="clustering method. hclust, kmeans", default='kmeans'),
  make_option(c("-b", "--bootstrapiter"), action="store", dest='bs.nrun', type="numeric", help="Number of bootstrap iterations.", default=20),
  make_option(c("-a", "--annotation"), action="store", dest='class.var', type="character", help="Name of column annotation fields in GCT file used to plot as consensus heatmap tracks.", default=''),
  #make_option(c("-x", "--firecloud"), action="store", dest='firecloud', type="logical", help="Running in FireCloud?", default=TRUE)
  make_option(c("-z", "--libdir"), action="store", dest='lib.dir', type="character", help="the src directory.", default='.')
  )
opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=TRUE)

source(file.path(opt$lib.dir, 'consensus_clustering.R'))


## ###########################################################################################
## main
main <- function(opt) {

  wd <- getwd()
  
  ## prepare log file
  logfile <- paste(opt$label, '_consensus_clustering.log', sep='')
  start.time <- Sys.time()
  cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'consensus_clustering\'--\n\n', file=logfile ) 
  cat('## parameters\ntar file:', opt$tar.file, '\ntmp dir:', opt$tmp.dir, '\nlabel:', opt$label, '\nlog file:', logfile, '\n', file=logfile, append=T)
  
  cat('current folder:', getwd(), '\n')
  
  ## ###############################################
  ## check whether input is a .tar archive
  suffix <- sub('.*\\.(.*)$', '\\1', opt$tar.file)
  if(suffix == 'tar'){
    
    ## extract tar ball
    if(!dir.exists(opt$tmp.dir))
      dir.create(opt$tmp.dir)
    cat('\n## Extracting tar file to', opt$tmp.dir, '\n', file=logfile, append=T)
    if(!dir.exists( file.path(opt$tmp.dir, opt$label) ))
      untar(opt$tar.file, exdir=opt$tmp.dir)
    
    ## cluster folder
    cluster.path.full <- file.path( opt$tmp.dir, opt$label, opt$clust.dir)
    
    ## path to GCT file
    gct.str <- file.path(opt$tmp.dir, opt$label, opt$filt.dir, glue('{opt$type}-ratio-norm-NArm.gct'))
    
    ## ################################################
    ## if the input is not a .tar file assume that theÍ
    ## function is called from the PGDAC-main pipeline
    } else if(suffix == 'gct'){
      
      gct.str <- opt$tar.file
      cluster.path.full <- opt$clust.dir
      
    } else {
      cat(glue('Using minimal sd of: {opt$clustering.sd.threshold}\n\n'))
      
        cluster.path.full <- getwd()
        gct.str <- file.path('..', opt$filt.dir, glue('{opt$type}-ratio-norm-NArm.gct'))
    }
  
  if(!dir.exists(cluster.path.full))
    dir.create(cluster.path.full)
  
 
  ## import data
  cat('\n## reading data from', gct.str, '\n', file=logfile, append=T)
  gct <- parse.gctx(gct.str)
  mat <- gct@mat
  rdesc <- gct@rdesc
  cdesc <- gct@cdesc
  
  ## #########################################
  ## eliminate features with not enough variation
  cat('\n## eliminating features with not enough variation ( standard deviation <', opt$clustering.sd.threshold, ')\n')
  cat('\n## eliminating features with not enough variation ( standard deviation <', opt$clustering.sd.threshold, ')\n', file=logfile, append=T)
  
  feature.sd <- apply (mat, 1, sd, na.rm=TRUE)
  keep <- which( feature.sd > opt$clustering.sd.threshold )
  
  cat(glue( '\nRemaining features used for clustering: {length(keep)}\n\n'))
  cat(glue( '\nRemaining features used for clustering: {length(keep)}\n\n'), file=logfile, append=T)
  
  ## #########################################
  ## require at least 200 features for clustering
  ## if fewer passing the SD threshold use the top 10%
  if(length(keep) < opt$min.feat.sd){
    
    upper.decile <- quantile(feature.sd, c(.9))
    cat('\n## WARNING! SD filtering returned to few features for clustering. Switching to upper decile mode keeping the top 10 percent (equivalent to SD of', upper.decile, ')\n')
    cat('\n## WARNING! SD filtering returned to few features for clustering. Switching to upper decile mode keeping the top 10 percent (equivalent to SD of', upper.decile, ')\n', file=logfile, append=T)
    
    keep <- which( feature.sd >  upper.decile)
    cat(glue( '\nRemaining features used for clustering: {length(keep)}\n\n'))
    cat(glue( '\nRemaining features used for clustering: {length(keep)}\n\n'), file=logfile, append=T)
  } 
  mat <- mat[keep, ]
  rdesc <- rdesc[keep, ]
 

  ## #########################################
  ## class variables of interest
  ## - used as annotation tracks in consensus heatmaps
  if( nchar(opt$class.var) == 0){ # | !(opt$class.var %in% colnames(cdesc))){
    cdesc.plot <- NULL
  } else {
    class.var <- strsplit( opt$class.var, ';') %>% unlist
    class.var <- class.var[which(class.var %in% colnames(cdesc))]
    
    if(length(class.var) > 0){
      
      cdesc.plot <- data.frame(cdesc[, class.var])
      colnames(cdesc.plot) <- class.var
    } else {
      cdesc.plot <- NULL
    }
  }
  
  ## #########################################
  ## run clustering
  setwd(cluster.path.full  )
  cons.clust <- consensus_clustering(mat,
                                     method=opt$method,
                                     nmf.nrun=2,
                                     bs.nrun=opt$bs.nrun,
                                     nmf.method='brunet',
                                     seed='random',      
                                     k.min=opt$k_min,            
                                     k.max=opt$k_max,
                                     make.nn=FALSE,
                                     ncore=detectCores(),
                                     plot=T,
                                     logfile=logfile,
                                     prefix=opt$type,
                                     cdesc=cdesc.plot
  )
  
  setwd(wd)
  invisible(file.copy(logfile, cluster.path.full) )
  
  return(0)
}

## ######################################################################
## run the pipeline
main(opt)

