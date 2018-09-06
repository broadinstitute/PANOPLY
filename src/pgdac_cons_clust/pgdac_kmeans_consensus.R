# #!/usr/bin/env Rscript
#options( warn = -1 )
library(optparse)

options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-i", "--input"), action="store", dest='tar.file' , type="character" ,help="input tar-file or (extracted) folder."),
  make_option(c("-l", "--label"), action="store", dest='label', type="character", help=""),
  make_option(c("-t", "--type"), action="store", dest='type', type="character", help="Data type, i.e. proteome, phospho, rna, ..."),
  make_option(c("-d", "--tempdir"), action="store", dest='tmp.dir', type="character", help="temp-folder", default="tmp"),
  make_option(c("-n", "--normdir"), action="store", dest='norm.dir', type="character", help="normalization-folder", default="normalized-data"),
  make_option(c("-c", "--clustdir"), action="store", dest='clust.dir', type="character", help="clustering-folder", default="clustering"),
  make_option(c("-u", "--k_min"), action="store", dest='k_min', type="numeric", help="Minimal cluster number.", default=2),
  make_option(c("-v", "--k_max"), action="store", dest='k_max', type="numeric", help="Maximal cluster number.", default=10),
  make_option(c("-s", "--sdfilter"), action="store", dest='clustering.sd.threshold', type="numeric", help="Minmal standard deviation acreoss samples required for clustering.", default=1.5),
  make_option(c("-m", "--method"), action="store", dest='method', type="character", help="clustering method. hclust, kmeans, nmf", default='kmeans'),
  make_option(c("-b", "--bootstrapiter"), action="store", dest='bs.nrun', type="numeric", help="Number of bootstrap iterations.", default=20),
  make_option(c("-a", "--annotation"), action="store", dest='class.var', type="character", help="Name of column annotation field of GCT file used to plat as heatmap track.", default=''),
  make_option(c("-z", "--libdir"), action="store", dest='lib.dir', type="character", help="the src directory.", default='.')
  
  )
opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=TRUE)

source(file.path(opt$lib.dir, 'consensus_clustering.R'))


# opt <- c()
# opt$tar.file = 'cluster-test.tar'
# opt$label = 'cluster-test'
# opt$k_min =2
# opt$k_max =4
# opt$type = 'proteome'
# opt$clust.dir <- 'cons-clust'
# opt$bs.nrun <- 100
# opt$tmp.dir = 'tmp'
# opt$norm.dir = 'normalized-data' 
# opt$clustering.sd.threshold <-0
#opt$class.var='class'


## ###########################################################################################
## main
main <- function(opt) {
  
  ## analysis_dir = opt$tmp.dir/opt$label
  
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
    gct.str <- file.path(opt$tmp.dir, opt$label, opt$norm.dir, glue('{opt$type}-ratio-norm-NArm.gct'))
    
    ## ################################################
    ## if the input is not a .tar file assume that the
    ## function is called from the PGDAC-main pipeline
    } else {
      ## source 'config.r' to get all parameters  
      #system(glue('R CMD BATCH --vanilla "--args {opt$type}" config.r;'))
      #source(glue('config.r'))
      
      #opt$clustering.sd.threshold <- clustering.sd.threshold
      cat(glue('test: {opt$clustering.sd.threshold}\n\n'))
      
      ## assume the current folder is the clustering-folder
      cluster.path.full <- getwd()
    
      ## path to GCT file
      gct.str <- file.path('..', opt$norm.dir, glue('{opt$type}-ratio-norm-NArm.gct'))
  }
  
  if(!dir.exists(cluster.path.full))
    dir.create(cluster.path.full)
  
 
  ## import data
  cat('\n## reading data from', gct.str, '\n', file=logfile, append=T)
  gct <- parse.gctx(gct.str)
  mat <- gct@mat
  rdesc <- gct@rdesc
  cdesc <- gct@cdesc
  
  ## eliminate features with not enough variation
  cat('\n## eliminating features with not enough variation ( standard deviation <', opt$clustering.sd.threshold, ')\n', file=logfile, append=T)
  feature.sd <- apply (mat, 1, sd, na.rm=TRUE)
  keep <- which( feature.sd > opt$clustering.sd.threshold )
  mat <- mat[keep, ]
  rdesc <- rdesc[keep, ]
  
  cat(glue( '\nRemaining features used for clustering: {length(keep)}\n'), file=logfile, append=T)
  
  
  
  ## class variable of interest
  if( nchar(opt$class.var) == 0 | !(opt$class.var %in% colnames(cdesc))){
    cdesc.plot <- NULL
  } else {
    
    class.var <- opt$class.var
    cdesc.plot <- data.frame(cdesc[, class.var])
    colnames(cdesc.plot) <- class.var
  }
  
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
  file.copy(logfile, cluster.path.full) 
}

## ######################################################################
## run the pipeline
main(opt)

