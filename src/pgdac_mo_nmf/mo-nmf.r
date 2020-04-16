#!/usr/bin/env Rscript
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("pacman"))

# specify command line arguments
option_list <- list(
  make_option( c("-t", "--tar"), action='store', type='character',  dest='tar.file', help='tar file containing data tables in GCT v1.3 format.'),
  make_option( c("-l", "--lowrank"), action='store', type='numeric',  dest='rank.low', help='Minimal factorization rank.', default = 2),
  make_option( c("-m", "--maxrank"), action='store', type='numeric',  dest='rank.high', help='Maximal factorization rank.', default = 4),
  make_option( c("-n", "--nrun"), action='store', type='numeric',  dest='nrun', help='Number of NMF runs with different starting seeds.', default = 5),
  #make_option( c("-b", "--bootstrap"), action='store', type='numeric',  dest='nrun.bs', help='Number of bootstrap iterations for consensus clustering.', default = 100),
  make_option( c("-b", "--bayesian"), action='store', type='logical',  dest='bnmf', help='If TRUE Bayesian NMF to determine optimal rank will be used (package ccfindR).', default = FALSE),
  make_option( c("-s", "--seed"), action='store', type='character',  dest='seed', help='Seed method for NMF factorization.', default = 'random'),
  #make_option( c("-c", "--classvariable"), action='store', type='character', dest='class.variable', help='Class variable of interest, NMF results will be mapped to levels of that variable. Must be present in the cdesc object of the GCT 1.3 file.'),
  make_option( c("-o", "--categorial_annotation"), action='store', type='character', dest='cat.anno', help='Categorial annotations used in enrichemnt analysis and annotation tracks. If multiple variables are of interest, separate by semicolon, e.g. -o "var1;var2;varN". First entry will be used as class variable. NMF results will be mapped to levels of that variable and thus is required.'),
  make_option( c("-c", "--continuous_annotation"), action='store', type='character', dest='cont.anno', help='Continuous annotations used in enrichment analysis and annotation tracks (optional).', default=NA),
  make_option( c("-d", "--define_colors"), action='store', type='character', dest='class.colors', help='Specifiy colors for each level in "categorial_annotation". Example: "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"', default = NA),
  make_option( c("-r", "--runonly"), action='store', type='logical', dest='no.plots', help='If TRUE no plots will be generated and the R-object after NMF clustering will be returned.', default=FALSE),
 # make_option( c("-g", "--genesofinterest"), action='store', type='character', dest='genes.of.interest', help='Character vector of gene symbols separated by semicolon. Example: "^PTPN1$;ERBB2;CDK5"', default=''),
  make_option( c("-f", "--sdfilter"), action='store', type='numeric', dest='sd.perc.rm', help='Lowest standard deviation across columns percentile to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest sd will be removed. Will be applied before z-scoring (optional)', default=0),
  make_option( c("-g", "--filt_mode"), action='store', type='character', dest='filt.mode', help='determines how the dataset will be filtered. XX', default='global'),
  make_option( c("-u", "--z_score"), action='store', type='logical', dest='zscore.all', help='If TRUE, rows in th matrix will be z-scored.', default=TRUE),
  make_option( c("-i", "--impute"), action='store', type='logical', dest='impute', help='If TRUE, the matrix will be first filtered for features present in >70% of samples. The remaining missing values will be imputed via KNN (R packe impute).', default=TRUE),
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene.column', help='Column name in rdesc in the GCT that contains gene names.', default='geneSymbol'),
  make_option( c("-e", "--exclude_2"), action="store", dest='exclude.2', type="logical", help="If TRUE, a number of clusters of will be excluded from calculation of the optimal rank.", default=TRUE),
  make_option( c("-z", "--libdir"), action="store", dest='lib.dir', type="character", help="the src directory.", default='/home/pgdac/src/'),
  make_option( c("-y", "--yaml"), action="store", dest='yaml.file', type="character", help="Path to .yaml file with parameters.", default=NA)
)

## #########################################################
# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

################################################
## funtion to parse yaml file  if specified 
## and update parameters
parse_param_yaml <- function(opt){
  
  if(!is.na(opt$yaml.file) > 0 & file.exists(opt$yaml.file)){
    require(pacman)
    p_load(yaml)
    p_load(magrittr)
    
    cat('\n\nFound yaml paramater file.\n - Updating parameters....')
    
    ##################################
    ## parse parameter file
    param <- try(read_yaml(opt$yaml.file))
    if(class(param) == 'try-error'){
      warning(paste0("Error parsing '",opt$yaml.file,"':\n", param, '\n\n'))
      return(opt)
    }
    param <- param$mo_nmf
    
    ##################################
    ## update parameters
    if(!is.null(param$tar_file)) opt$tar.file <- param$tar_file
    if(!is.null(param$kmin)) opt$rank.low <- param$kmin
    if(!is.null(param$kmax)) opt$rank.high <- param$kmax
    if(!is.null(param$exclude_2)) opt$exclude.2 <- param$exclude_2
    if(!is.null(param$core_membership)) opt$core.membership <- param$core_membership
    ## NMF specific
    if(!is.null(param$nrun)) opt$nrun <- param$nrun
    if(!is.null(param$method)) opt$nmf.method <- param$method
    if(!is.null(param$bnmf)) opt$bnmf <- param$bnmf
    if(!is.null(param$seed)) opt$seed <- param$seed
    ## annotations
    if(!is.null(param$cat_anno))  opt$cat.anno <- param$cat_anno %>% unlist %>% paste(., collapse=';')
    if(!is.null(param$cont_anno))  opt$cont.anno <- param$cont_anno %>% unlist %>% paste(., collapse=';')
    ## colors
    if(!is.null(param$colors)){
      colors.tmp <- sapply(param$colors, function(x) paste( unlist(paste(names(x), x, sep  =':')), collapse = ';' ))
      class.colors <- paste(paste(names(colors.tmp), colors.tmp, sep='='), collapse='|')
      opt$class.colors <- class.colors
    }
    if(!is.null(param$sd_filt)) opt$sd.perc.rm <- param$sd_filt
    if(!is.null(param$filt_mode)) opt$filt.mode <- param$filt_mode
    if(!is.null(param$z_score)) opt$zscore.all <- param$z_score
    ## imputation
    if(!is.null(param$impute)) opt$impute <- param$impute
    if(!is.null(param$max_na_row)) opt$max.na.row <- param$max_na_row
    if(!is.null(param$max_na_col)) opt$max.na.col <- param$max_na_col
    ## heatmap parameters
    if(!is.null(param$hm_cw)) opt$hm.cw <- param$hm_cw
    if(!is.null(param$hm_ch)) opt$hm.ch <- param$hm_ch
    if(!is.null(param$hm_max_val)) opt$hm.max.val <- param$hm_max_val
    if(!is.null(param$hm_max_val_z)) opt$hm.max.val.z <- param$hm_max_val_z
    ## misc
    if(!is.null(param$gene_col)) opt$gene.column <- param$gene_col
    if(!is.null(param$nmf_only)) opt$no.plots <- param$nmf_only
    if(!is.null(param$lib_dir)) opt$lib.dir <- param$lib_dir
    cat('done\n\n')
  } 
  
  return(opt)
} 

############################################################
## parse yaml
p_load(yaml)
#opt$yaml.file <- 'param.yaml'
opt <- parse_param_yaml(opt)
save(opt, file='options.RData')

########################
## import libraries
p_load(NMF)
p_load(ccfindR)
p_load(RSQLite)
p_load(org.Hs.eg.db)
p_load(pheatmap)
p_load(WriteXLS)
p_load(RColorBrewer)
p_load(rhdf5)
if(!require(cmapR))devtools::install_github("cmap/cmapR")
p_load(dplyr)
p_load(ggpubr)
p_load(Rtsne)
p_load(plotly)
p_load(rmarkdown)
p_load(maptools)
p_load(UpSetR)
p_load(ComplexHeatmap)
p_load(glue)
p_load(circlize)
p_load(impute)

p_load(doParallel)
p_load(foreach)

p_load(readxl)
p_load(tibble)
p_load(magrittr)

p_load(limma)

#opt$lib.dir <- 'c:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_mo_nmf/'

################################################
## source required files
cat('Sourcing files...')
source(file.path(opt$lib.dir, 'my_plots.r'))
source(file.path(opt$lib.dir, 'my_io.r'))
source(file.path(opt$lib.dir, 'nmf_functions.R'))
source(file.path(opt$lib.dir, 'modT.r'))
cat('done!\n')

## #############################################
## main function
main <- function(opt){

    opt$sd.perc.rm <- as.numeric(opt$sd.perc.rm)
    
    ## ##################################################
    ##           parse parameters
    ## ###################################################
    tar.file <- opt$tar.file
  
    ## percentile of stddev to remove from the data
    sd.perc.rm <- opt$sd.perc.rm
    ## filter mode: global, sepaarte, equal
    filt.mode <- opt$filt.mode
      
     ## apply z-score to all data types?
    zscore.all <- opt$zscore.all
    
    ## #############################
    ## imputation parameters (KNN)
    impute <- opt$impute
    impute.k <- 5
    na.max.row <- 0.3
    na.max.col <- 0.9
     
    ## ################################
    ##            NMF parameters
    ranks <- sort(opt$rank.low:opt$rank.high)
    seed <- opt$seed
    nrun <- opt$nrun
    #nrun.bs <- opt$nrun.bs ## bootstrap iterations
    method <- ifelse(opt$bnmf, 'bnmf', 'nmf')
      
    cores <- detectCores()
    opts <- paste('vp', cores,'t', sep='')
  
          
    ## ###################################################
    ##      other parmaters
    ## ###################################################
   # tmp.dir <- 'tmp'
    tmp.dir <- tempdir()
    
    ## data filtering / normalization
    zscore.cnv <-  F     ## should CNVs be z-scored?
 

    ## #########################################################################################
    ##
    ##                               START
    ##
    ## #########################################################################################
   
    ## ##################################################
    ## import data files
    data.imported <- import.data.sets(tar.file, tmp.dir, zscore.cnv)
    
    expr <- data.imported$expr
    cdesc <- data.imported$cdesc
    rdesc <- data.imported$rdesc
    data.str <- data.imported$data.str
    
    ## data types
    data.type.dist <- table(sub('^(.*?)-.*','\\1', rdesc$Data.Type.ID))
    
    ## ###############################
    ##     prepare output folder
    label <- paste(names(data.str), collapse='_')
    
    date.str <- glue("{sub(' .*', '', Sys.time())}_{label}_{nrun}_{method}{ifelse(impute,'_knn_','_')}sd_{sd.perc.rm}{ifelse(zscore.all,'_zcore','')}")
    dir.create(date.str)
    setwd(date.str)
    
    ## ###############################
    ##     generate parameter file
    param <- c('## NMF parameters:', 
               paste('method:', method), 
               
               paste('rank:', paste(ranks, collapse=',')), 
               paste('seed:',seed), paste('nrun:', nrun), 
              
               paste('opts:', opts),
               '',
               '## data filter / normalization / imputation', 
               paste('filter mode:', filt.mode),
               paste('remove lowest StdDev perc.:', sd.perc.rm),
               
               paste('z-score rows:', zscore.all), 
               paste('imputation (KNN):', impute),
               paste('imputation: no. neighbors (K)', impute.k),
               paste('imputation: max.missing row', na.max.row),
               paste('imputation: max.missing col', na.max.col),
               '',
               paste('## data tables'), 
               paste(data.str) )
    
    writeLines(param, con='parameters.txt')
    
    ## ####################################################
    ## KNN imputation
    if(impute){
       keep.idx <- which(apply(expr, 1, function(x) ifelse( sum(is.na(x))/length(x) < na.max.row, T, F)))
       expr <- expr[keep.idx, ]
       rdesc <- rdesc[rownames(expr), ]
       ## impute
       expr.imp <- try(impute.knn(expr, k = impute.k, colmax = na.max.col ))
       expr <- expr.imp$data
       
    } else { ## keep only fully quantified features
        keep.idx <- which(apply(expr, 1, function(x) ifelse(sum(is.na(x)) > 0, F, T)))
        expr <- expr[keep.idx, ]
        rdesc <- rdesc[rownames(expr), ]
    }
    
    ## ############################################
    ##  export merged data used for nmf as gct
    gct.comb <- new('GCT')
    gct.comb@mat <- data.matrix(expr)
    gct.comb@cdesc <- cdesc
    gct.comb@rdesc <- rdesc
    gct.comb@cid <- colnames(gct.comb@mat)
    gct.comb@rid <- rownames(gct.comb@mat)
    write.gct(gct.comb, ofile = 'mo-data-matrix')
    
    ## ###################################################
    ##       calculate SD accross samples
    ## ###################################################
    ## sd.expr <- apply(gct.comb@mat, 1, sd, na.rm=T)
    
    ## ####################################################
    ##    sd filter
    gct.filt <- filter.datasets(gct.comb, sd.perc.rm, mode=filt.mode )
    write.gct(gct.filt, ofile = 'mo-data-matrix-sd-filt')
    
    ## update
    expr <- gct.filt@mat

    ## ###################################################
    ## Z-score
    ## ###################################################
    if(zscore.all)
      expr <- scale(expr)
   
    # ## #########################
    # ## density kernel
    # sd.d <- density( sd.expr )
    # sd.d.list <- list()
    # ymax <- 0
    # for(i in 1:length(data.str)){
    #     sd.d.list[[i]] <- density( sd.expr[ grep(names(data.str)[i], names(sd.expr)) ] )
    #     if(max(sd.d.list[[i]]$y) > ymax)
    #         ymax <- max(sd.d.list[[i]]$y)
    # }
    # ymax <- (ymax+(0.1*ymax))
    # 
    # ## #####################
    # ## plot
    # fn <- '0_density_StdDev.pdf'
    # pdf(fn)
    # plot(sd.d, ylim=c(0, ymax), main='StdDev accross samples', lwd=2)
    # for(i in 1:length(data.str))
    #     lines(sd.d.list[[i]], col=i+1, lwd=2)
    # legend('topright', legend=c('Aggregated', names(data.str)), lty='solid', lwd=3, col=1:(length(data.str)+1))
    # 
    # 
    # ####################################################
    # ## filter dataset
    # 
    # if(sd.perc.rm > 0){
    #   
    #     cat("\napplying SD-filter: ")
    #     sd.perc <- quantile(sd.expr, c(sd.perc.rm))
    #     sd.keep <- which( sd.expr > sd.perc)
    #     cat("removed", nrow(expr)-length(sd.keep), "features with SD<=", round(sd.perc, 2), "(", names(sd.perc),"-tile)\n\n")
    #     expr <- expr[sd.keep,]
    #     
    #     ## add line to plot
    #     abline(v=sd.perc, lty='dashed', lwd=2, col='grey')
    #     text(sd.perc, ymax, paste(100*sd.perc.rm, 'th percentile', sep=''), pos=4, col='grey')
    #     
    #     ## export filtered GCT
    #     gct.filt <- gct.comb
    #     gct.filt@mat <- data.matrix(expr)
    #     gct.filt@rdesc <- rdesc[sd.keep, ]
    #     gct.filt@rid <- rownames(gct.filt@mat)
    #     write.gct(gct.filt, ofile = 'mo-data-matrix-sd-filt')
    #     
    # }
    # dev.off()
    # #im.convert(fn, output = sub('\\.pdf','.png', fn), extra.opts="-density 150")
    
    
    
    ## ##############################################
    ##           make matrix non-negativ
    ## - separate up/down
    ## ##############################################
    expr.comb <- make.non.negative(expr)
    
    ## remove rows only containing zero values
    keep.idx <- which(apply(expr.comb, 1, function(x) sum(x != 0) ) > 0)
    expr.comb <- expr.comb[keep.idx, ]
    
    
    ## ###################################################
    ##                run NMF
    ## ###################################################
    res.rank <- vector('list', length(ranks))
    names(res.rank) <- ranks
    
    ##########################################
    ##               bayesian NMF
    if(opt$bnmf){
      
      #cl <- makeCluster(cores)
      #registerDoParallel(cl)
      
      ## loop over ranks
      for(i in 1:length(ranks)){
      #res.rank <- foreach(i = 1:length(ranks)) %dopar% {
      #  library(ccfindR)
        res.tmp <- vb_factorize(scNMFSet( count=expr.comb), ranks = ranks[i], nrun = nrun, 
                                verbose = 1, Tol=1e-5, progress.bar = T, ncores=cores,
                                initializer='seed')
        res.tmp
      }
      #on.exit(stopCluster(cl))
      #names(res.rank) <- ranks
      
        
      if('evidence' %in% names(res.rank[[1]]@measure))  
        rank.mev <- sapply(res.rank, function(x) x@measure$evidence)
      else
        rank.mev <- sapply(res.rank, function(x) x@measure$lml) ## names changed in newer verion of the package...
      
        names(rank.mev) <- as.character(ranks)
        #save(rank.mev, res.rank, file='debug.RData')
        ## ####################################################
        ##      Estimate factorization rank
        ## ####################################################
        
        ## extract top N ranks based on maximum evidence
        topn.rank <- 1
        exclude.2 <- F ## should rank=2 be excluded?
        
        rank.top <- GetBestRank(rank.mev, topn.rank=topn.rank, exclude.2=exclude.2)
        rank.sil <- rank.sil.avg <- rank.coph <- NULL
        
        ###########################################################   
        ## plot 'max. evidence'
        pdf('1_cluster_metrics.pdf', 4, 4)
        
        col.tmp <- rep('darkblue', length(rank.mev))
        names(col.tmp) <- names(rank.mev)
        col.tmp[rank.top] <- 'red'
        
        try(plot(ranks, rank.mev, xlab='rank', ylab='log(evidence)', type='b', axes=F, col=col.tmp))
        try(axis(1, at=ranks))
        try(axis(2))
        
      dev.off()
        
      ##################################################
      ##              non-bayesian NMF
    } else {
        
      ## loop over ranks
      for(i in 1:length(ranks))
          res.rank[[i]] <-  nmf(expr.comb, rank=ranks[i], method='brunet', seed=seed, nrun=nrun, .options = opts)
      
      ## ####################################################
      ##             calculate cluster metrics
      ## ####################################################
      
      ## SVD
      rank.svd <- svd(expr.comb)
      
      ## silhouette
      rank.sil <- lapply(res.rank, silhouette)
      rank.sil.avg <- lapply(rank.sil, function(x) tapply( x[,3], x[, 1], mean))
      
      ## cophenic
      rank.coph <- sapply(res.rank, cophcor)
      
      ## dispersion of sonsensus matrix
      rank.disp <- sapply(res.rank, dispersion)
      
      ## combine
      rank.coph.disp <- rank.coph * rank.disp
      
      ## ####################################################
      ##      Estimate factorization rank
      ## ####################################################
      
      ## extract top N ranks based on cophentic correlation
      topn.rank <- 1
      
      #rank.top <- GetBestRank(rank.coph, topn.rank=topn.rank, exclude.2=opt$exclude.2)
      rank.top <- GetBestRank(rank.coph.disp, topn.rank=topn.rank, exclude.2=opt$exclude.2)
      rank.mev  <- NULL
      
      #########################################
      ## plot cluster metrics
      if(length(ranks) > 1){
        
        col.tmp <- rep('darkblue', length(rank.coph))
        names(col.tmp) <- names(rank.coph)
        col.tmp[rank.top] <- 'red'
        
        
        ## plot
        pdf('1_cluster_metrics.pdf', 5, 5)
        
        ## SVD
        plot(rank.svd$d, cex=2, type='b', pch=20, main='SVD', xlab= 'Sample index', ylab='Singular values of x', col='darkblue')
        
        ## cophenetic
        plot(as.numeric(names(rank.coph)), rank.coph, type='b', pch=20, cex=2, col=col.tmp, xaxt='n', xlab='Factorization rank / No. of clusters', ylab='Cophenetic score', main='Cophenetic correlation coefficient', 
             ylim=c( min(rank.coph)-(min(rank.coph)*0.1), 1) )
        axis(1, at=as.numeric(names(rank.coph)))
        
        ## dispersion
        plot(as.numeric(names(rank.disp)), rank.disp, type='b', pch=20, cex=2, col=col.tmp, xaxt='n', xlab='Factorization rank / No. of clusters', ylab='Dispersion coefficient', main='Dispersion score of consensus matrix', 
             ylim=c( min(rank.disp)-(min(rank.disp)*0.1), 1) )
        axis(1, at=as.numeric(names(rank.disp)))
        
        ##################################################################
        ## combined plot
        ylim=c( range( c(rank.coph, rank.disp, rank.coph.disp)) )
        col <- c('coph'=my.col2rgb('black'), 'disp'=my.col2rgb('darkblue'), 'comb'=my.col2rgb('red'))#, 'sil'='grey')
        pch <- c('coph'=17, 'disp'=20, 'comb'=18)
        lty <- c('coph'='dashed', 'disp'='dashed', 'comb'='solid')
        cex <- 2.2
        
        plot(names(rank.coph), rank.coph, ylab='Score', xlab='Factorization rank', col=col['coph'], pch=pch['coph'], lty=lty['coph'], cex=cex, type='b', main='Cluster metrics', ylim=ylim,
             cex.axis=1.5, cex.lab=1.5)
        lines(names(rank.disp), rank.disp, col=col['disp'], pch=pch['disp'], lty=lty['disp'], cex=cex, type='b')
        lines(names(rank.coph.disp), rank.coph.disp, col=col['comb'], pch=pch['comb'], lty=lty['comb'], type='b', cex=cex)
        
        legend('bottomleft', legend=names(col), col=col, pch=pch, bty='n', cex=1.1)
        
        rt <- as.numeric(rank.top)
        rect( xleft = rt-0.2, xright = rt+0.2, ybottom = ylim[1]-.01, ytop = ylim[2]+0.01, border = 'red', lwd=1.6 )
        
        
        
        ## silhouette
        # par(mfrow=c(2,1))
        #fancyBoxplot( lapply(rank.sil, function(x)x[, 3] ), xlab='Factorization rank / No. of clusters',
        #              main='Silhouette scores', ylab='Silhouette scores', show.numb = 'mean', numb.cex = 1.5, numb.pos = 4, 
        #              col='white', box.border = 'darkblue', vio.alpha = 0, lwd=2.5)
        #fancyBoxplot( lapply(rank.sil.random, function(x)x[, 3] ), xlab='Factorization rank / No. of clusters', main='Randomized data', ylab='Silhouette scores', show.numb = 'mean', numb.cex = 1.5, numb.pos = 4)
        dev.off()
      }
      
         
    }
        
    ## #####################################################
    ## save results needed to 
    fn.ws <- 'workspace_after_NMF.RData'
    save(opt, res.rank, rank.top, expr.comb, expr, rdesc, cdesc, rank.sil, rank.sil.avg, 
         rank.coph, rank.mev, fn.ws, gct.comb, zscore.all, file=fn.ws)
    
    
    ## #######################################################
    ##         generate some plots
    if(opt$no.plots == FALSE){ 
      p <- try(nmf.post.processing(fn.ws))
    }
    
    return(0)
} ## end main()
    
   
## #########################################################
##               run pipeline
main(opt)


