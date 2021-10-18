#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))

# specify command line arguments
option_list <- list(
  make_option( c("-t", "--tar"), action='store', type='character',  dest='tar_file', help='tar file containing data tables in GCT v1.3 format.'),
  make_option( c("-l", "--lowrank"), action='store', type='numeric',  dest='kmin', help='Minimal factorization rank.'),  # default = 2),
  make_option( c("-m", "--maxrank"), action='store', type='numeric',  dest='kmax', help='Maximal factorization rank.'), # default = 4),
  make_option( c("-n", "--nrun"), action='store', type='numeric',  dest='nrun', help='Number of NMF runs with different starting seeds.'), # default = 5),
  make_option( c("-b", "--bayesian"), action='store', type='logical',  dest='bnmf', help='If TRUE Bayesian NMF to determine optimal rank will be used (package ccfindR).', default = FALSE),
  make_option( c("-s", "--seed"), action='store', type='character',  dest='seed', help='Seed method for NMF factorization.'), # default = 'random'),
  make_option( c("-o", "--categorial_annotation"), action='store', type='character', dest='cat_anno', help='Categorial annotations used in enrichemnt analysis and annotation tracks. If multiple variables are of interest, separate by semicolon, e.g. -o "var1;var2;varN". First entry will be used as class variable. NMF results will be mapped to levels of that variable and thus is required.'),
  make_option( c("-c", "--continuous_annotation"), action='store', type='character', dest='cont_anno', help='Continuous annotations used in enrichment analysis and annotation tracks (optional).') , # default=NA),
  make_option( c("-d", "--define_colors"), action='store', type='character', dest='cat_colors', help='Specifiy colors for each level in "categorial_annotation". Example: "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"'), # default = NA),
  make_option( c("-r", "--runonly"), action='store', type='logical', dest='nmf_only', help='If TRUE no plots will be generated and the R-object after NMF clustering will be returned.'), # default=FALSE),
  make_option( c("-f", "--sdfilter"), action='store', type='numeric', dest='sd_filt', help='Lowest standard deviation across columns percentile to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest sd will be removed. Will be applied before z-scoring (optional)') , ## default=0),
  make_option( c("-g", "--filt_mode"), action='store', type='character', dest='filt_mode', help='Determines how the low variable features (paramater "sd_filt") will be removed from the dataset. "global": the filter will be applied to the entrire datasets. "separate": the will be applied separately to each data type (e.g. proteome, RNA, etc.). "equal": all data tables of different data types (e.g. proteome, RNA) will be filtered down to have the same number of features and ths is limited by the data type with the smallest number of features "N_feat". The data tables of other data types are filtered to retain the "N_feat" highest variable features.  Only relevant for multi-omics data.'), # default='global'),
  make_option( c("-u", "--z_score"), action='store', type='logical', dest='z_score', help='If TRUE, the data matrix will be z-scored (see --z_score_mode).'), # default=TRUE),
  make_option( c("-v", "--z_score_mode"), action='store', type='character', dest='z_score_mode', help='z-scoring mode: row, col, rowcol'), # default=TRUE),
  make_option( c("-i", "--impute"), action='store', type='logical', dest='impute', help='If TRUE, the matrix will be first filtered for features present in >70% of samples. The remaining missing values will be imputed via KNN (R packe impute).'), # default=TRUE),
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  make_option( c("-e", "--exclude_2"), action="store", dest='exclude_2', type="logical", help="If TRUE, a '2' will be excluded from calculation of the optimal rank."), #default=TRUE),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/home/pgdac/src/'),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters.", default=NA)
)

## #########################################################
# parse command line parameters
# - do it here already to speed up the USAGE message
# - the actual parsing of parameters happens in 
#   function 'parse_yaml_mo_nmf()'
opt_cmd <- parse_args( OptionParser(option_list=option_list) )

###########################################################
## load libraries requires to parse and update parameters
library("pacman")
p_load("yaml")
p_load("magrittr")

################################################
## funtion to parse and update parameters
## - cmd line
## - yaml file
## parameters in yaml file will be updated with 
## parameters specified on cmd
parse_yaml_mo_nmf <- function(cmd_option_list, 
                              yaml_section='panoply_mo_nmf', 
                              yaml_colors='groups.colors', 
                              yaml_groups_cat='groups.cols', 
                              yaml_groups_cont='groups.cols.continuous'
                              ){
  
  ## #########################################################
  # parse command line parameters
  opt_cmd <- parse_args( OptionParser(option_list=cmd_option_list) )
  #opt_cmd$yaml_file <- "master-parameters.yaml"
  
  ############################################################
  ## parse yaml file
  if(!is.na(opt_cmd$yaml_file) & file.exists(opt_cmd$yaml_file)){
    
    ## import yaml
    opt_yaml_all <- read_yaml(opt_cmd$yaml_file)
    
    ###################################
    ## parse groups
    opt_yaml_groups <- opt_yaml_all[[yaml_groups_cat]] %>% unlist %>% paste(., collapse=';')
    names(opt_yaml_groups)[length(opt_yaml_groups)] <- 'cat_anno'
    
    ## optional continuous variables
    if(yaml_groups_cont %in% names(opt_yaml_all)){
      opt_yaml_groups[['cont_anno']] <- opt_yaml_all[[yaml_groups_cont]] %>% unlist %>% paste(., collapse=';')
    } else {
      opt_yaml_groups[['cont_anno']] <- NA
    }
    
    ###################################
    ## parse colors
    tmp <- opt_yaml_all[[yaml_colors]]
    #colors.tmp <- sapply(opt_yaml[["cat_colors"]], function(x) paste( unlist(paste(names(x), x, sep  =':')), collapse = ';' ))
    colors.tmp <- sapply(tmp, function(x) paste( unlist(paste(names(x), x, sep  =':')), collapse = ';' ))
    opt_yaml_colors <- list()
    opt_yaml_colors[['cat_colors']] <- paste(paste(names(colors.tmp), colors.tmp, sep='='), collapse='|')
    
    
    ###################################
    ## extract relevant section
    opt_yaml <- opt_yaml_all[[yaml_section]]
    
    ## parse cmd params
    cat('\n\nparsing command line parameters:\n')
    for(x in names(opt_cmd))
     cat('---', x, opt_cmd[[x]], '; prefer yaml?', opt_cmd[[x]] == 'NA','\n')
    
    ## update yaml with parameters specified on cmd line 
    cat('\n\nUpdating parameter file with command line parameters:\n')
    ## cmd parameters
    cmd_not_null <- which( !sapply(opt_cmd, function(x) x == 'NA' ) )
    cmd_to_update <- intersect( names(opt_cmd)[ cmd_not_null], names(opt_yaml) )
    cmd_to_add <- setdiff( names(opt_cmd), names(opt_yaml) )
    
    ## update yaml by cmd 
    opt_yaml[cmd_to_update] <- opt_cmd[cmd_to_update]
    sapply(cmd_to_update, function(x) cat(x, '\n'))
    
    ## add parameters only specified on cmd
    if(length(cmd_to_add) > 0){
      opt_cmd_to_add <- opt_cmd[cmd_to_add]
      opt_yaml <- append(opt_yaml, opt_cmd_to_add)
    }
    
    ## append grroups and colors
    opt_yaml <- append(opt_yaml, opt_yaml_groups)
    opt_yaml <- append(opt_yaml, opt_yaml_colors)
    
    
    ## updated params
    opt <- opt_yaml
    
  } else {
    ## no yaml file
    opt <- opt_cmd
  }
  
  ## force correct mode
  opt$sd_filt <- as.numeric(opt$sd_filt)
  opt$core_membership <- as.numeric(opt$core_membership)
  opt$core_membership_mode <- as.character(opt$core_membership_mode)
    
  opt$z_score <- as.logical(opt$z_score)
  opt$impute <- as.logical(opt$impute)
  opt$exclude_2 <- as.logical(opt$exclude_2)
  opt$bnmf  <- as.logical(opt$bnmf)
  opt$nmf_only  <- as.logical(opt$nmf_only)
  opt$method <- as.character(opt$method)
    
  opt$max_na_row <- as.numeric(opt$max_na_row)
  opt$max_na_col  <- as.numeric(opt$max_na_col)
  opt$hm_cw <- as.numeric(opt$hm_cw)
  opt$hm_ch <- as.numeric(opt$hm_ch)
  opt$hm_max_val <-  as.numeric(opt$hm_max_val)
  opt$hm_max_val_z <- as.numeric(opt$hm_max_val_z)
  
  opt$impute_k <- as.integer(opt$impute_k)
  opt$kmin <- as.integer(opt$kmin)
  opt$kmax <- as.integer(opt$kmax)
  opt$nrun <- as.integer(opt$nrun)
  
  opt$filt_mode <- as.character(opt$filt_mode)
  opt$gene_col <- as.character(opt$gene_col)
  opt$organism <- as.character(opt$organism)
  opt$z_score_mode <- as.character(opt$z_score_mode)
  
  opt$feature_fdr <- as.numeric(opt$feature_fdr)
  opt$ora_pval <- as.numeric(opt$ora_pval)
  opt$ora_max_categories <- as.integer(opt$ora_max_categories)
  
  return(opt)
} 

########################
## import libraries
p_load(NMF)
#p_load(ccfindR)
p_load(RSQLite)

## annotation packages
p_load(org.Hs.eg.db)
p_load(org.Rn.eg.db)
p_load(org.Mm.eg.db)

p_load(pheatmap)
p_load(WriteXLS)

p_load(RColorBrewer)
p_load(pals)

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
p_load(circlize)

p_load(glue)

p_load(impute)

p_load(doParallel)
p_load(foreach)

p_load(readxl)
p_load(tibble)
p_load(magrittr)

p_load(limma)

## parse parameters
opt <- parse_yaml_mo_nmf(option_list) 

################################################
## source required files
cat('\nSourcing files...')
source(file.path(opt$lib_dir, 'my_plots.r'))
source(file.path(opt$lib_dir, 'my_io.r'))
source(file.path(opt$lib_dir, 'nmf_functions.R'))
source(file.path(opt$lib_dir, 'modT.r'))
cat('done!\n')

## #############################################
## main function
main <- function(opt){

    ## set defaults for backwards compatibility
    if(is.null(opt$organism)) opt$organism <- 'human'
    if(is.null(opt$blank_anno)) opt$blank_anno='N/A'
    if(is.null(opt$blank_anno_col)) opt$blank_anno_col='white'
    if(is.null(opt$feature_fdr)) opt$feature_fdr <- 0.01
    if(is.null(opt$ora_pval)) opt$ora_pval <- 0.01
    if(is.null(opt$ora_max_categories)) opt$ora_max_categories <- 10
    if(is.null(opt$core_membership_mode)) opt$core_membership_mode <- 'legacy'
    if(is.null(opt$sample_order)) opt$sample_order <- ifelse(opt$core_membership_mode == 'mindiff', 'clust_then_core_then_score', 'clust_then_score')
    
    ## ##################################################
    ##           parse parameters
    ## ###################################################
    tar_file <- opt$tar_file
    
    ## percentile of stddev to remove from the data
    sd_filt <- opt$sd_filt

    ## filter mode: global, sepaarte, equal
    filt_mode <- opt$filt_mode
      
    ## apply z-score to all data types?
    z_score <- as.logical(opt$z_score)
    
    ## #############################
    ## imputation parameters (KNN)
    impute <- opt$impute
    impute.k <- opt$impute_k
    na.max.row <- opt$max_na_row ##0.3
    na.max.col <- opt$min_na_row ## 0.9
     
    ## ################################
    ##            NMF parameters
    ranks <- sort(opt$kmin:opt$kmax)
    seed <- opt$seed
    nrun <- opt$nrun
    #nrun.bs <- opt$nrun.bs ## bootstrap iterations
    #nmf_method <- ifelse(opt$bnmf, 'bnmf', 'nmf')
    nmf_method <- opt$method 
      
    cores <- detectCores()
    opts <- paste('vp', cores,'t', sep='')
  
          
    ## ###################################################
    ##      other parmaters
    ## ###################################################
    tmp.dir <- tempdir()
    
    ## data filtering / normalization
    #zscore.cnv <-  F     ## should CNVs be z-scored?
 
    ## #########################################################################################
    ##
    ##                               START
    ##
    ## #########################################################################################
   
    ## ##################################################
    ## import data files
    data.imported <- import.data.sets(tar_file, tmp.dir)
    
    expr <- data.imported$expr
    cdesc <- data.imported$cdesc
    rdesc <- data.imported$rdesc
    data.str <- data.imported$data.str
    
    ## data types
    data.type.dist <- table(sub('^(.*?)-.*','\\1', rdesc$Data.Type.ID))
    
    ## ###############################
    ##     prepare output folder
    label <- paste(names(data.str), collapse='_')
    
    date.str <- glue("{sub(' .*', '', Sys.time())}_{label}_{nrun}_{make.names(nmf_method)}{ifelse(impute,'_knn_','_')}sd_{sd_filt}{ifelse(z_score,paste0('_zcore_', opt$z_score_mode),'')}")
    
    dir.create(date.str)
    setwd(date.str)
    
    ## ###############################
    ##     generate parameter file
    param <- c('## NMF parameters:', 
               paste('nmf_method:', nmf_method), 
               
               paste('rank:', paste(ranks, collapse=',')), 
               paste('seed:',seed), paste('nrun:', nrun), 
              
               paste('opts:', opts),
               '',
               paste('core membership method:', opt$core_membership_mode),
               '',
               '## data filter / normalization / imputation', 
               paste('filter mode:', filt_mode),
               paste('remove lowest StdDev perc.:', sd_filt),
               '',
               '## normalization',
               paste('z-score:', z_score), 
               paste('z-score mode:', opt$z_score_mode),
               '',
               '## imputation', 
               paste('imputation (KNN):', impute),
               paste('imputation: no. neighbors (K)', impute.k),
               paste('imputation: max.missing row', na.max.row),
               paste('imputation: max.missing col', na.max.col),
               '',
               paste('## data tables'), 
               paste(data.str) )
    
    writeLines(param, con='parameters.txt')
    
    ## export yaml
    write_yaml(opt, 'parameters.yaml')
    
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
    ##    sd filter
    ######################################################
    gct.filt <- filter.datasets(gct.comb, sd_filt, mode=filt_mode )
    write.gct(gct.filt, ofile = 'mo-data-matrix-sd-filt')
    
    ## update
    expr <- gct.filt@mat

    ## ###################################################
    ## Z-score
    ## ###################################################
    
    pdf('0_boxplot-data-distribution.pdf')
    expr_org <- expr
    boxplot(expr_org, cex=0.5, main='before z-scoring')
    
    if(z_score){
      
      ## row
      expr_z_row <- t(scale(t(expr_org))) ## row
      boxplot(expr_z_row, cex=0.5, main='after z-scoring (row)', outline=F)
      
      ## column
      expr_z_col <- scale(expr_org)
      boxplot(expr_z_col, cex=0.5, main='after z-scoring (column)', outline=F)

      ## row + column
      expr_z_rowcol <- scale(expr_z_row) ## row + column
      boxplot(expr_z_rowcol, cex=0.5, main='after z-scoring (row + column)', outline=F)
     
      
      if(opt$z_score_mode == 'rowcol')
        expr <- scale(expr_z_rowcol)
      if(opt$z_score_mode == 'row')
        expr <- expr_z_row
      if(opt$z_score_mode == 'col')
        expr <- expr_z_col
      
      ###############################
      ## export
      gct.filt.zscore <- gct.filt
      gct.filt@mat <- expr
      write.gct(gct.filt.zscore, ofile = paste('mo-data-matrix-sd-filt-zscore', opt$z_score_mode, sep='-'))
    }
    dev.off()
    
    ## ##############################################
    ##           make matrix non-negative
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
    ## - only applicable to countr data
    ## - deactivated in PANOPLY/TERRA
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
        rank.mev <- sapply(res.rank, function(x) x@measure$lml) ## names changed in newer version of the package...
        names(rank.mev) <- as.character(ranks)
    
        ## ####################################################
        ##      Estimate factorization rank
        ## ####################################################
        
        ## extract top N ranks based on maximum evidence
        topn.rank <- 1
        exclude_2 <- F ## should rank=2 be excluded?
        
        rank.top <- GetBestRank(rank.mev, topn.rank=topn.rank, exclude_2=exclude_2)
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
      ##              non-Bayesian NMF
      ##################################################
    } else {
        
      ## loop over ranks
      for(i in 1:length(ranks))
          res.rank[[i]] <-  nmf(expr.comb, rank=ranks[i], method=nmf_method, seed=seed, nrun=nrun, .options = opts)
      
      ## ####################################################
      ##             calculate cluster metrics
      ## ####################################################
      
      ## SVD
      rank.svd <- svd(expr.comb)
      
      ## silhouette
      rank.sil <- lapply(res.rank, silhouette)
      rank.sil.avg <- lapply(rank.sil, function(x) tapply( x[,3], x[, 1], mean))
      
      ## cophenetic
      rank.coph <- sapply(res.rank, cophcor)
      
      ## dispersion of consensus matrix
      rank.disp <- sapply(res.rank, dispersion)
      
      ## combine
      #rank.coph.disp <- rank.coph * rank.disp
      #rank.coph.disp <- rank.coph^(1-rank.disp)
      rank.coph.disp <- rank.disp^(1-rank.coph)
      
      ## ####################################################
      ##      Estimate factorization rank
      ## ####################################################
      
      ## extract top N ranks based on cophenetic correlation
      topn.rank <- 1
      
      rank.top <- GetBestRank(rank.coph.disp, topn.rank=topn.rank, exclude_2=opt$exclude_2)
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
         rank.coph, rank.disp, rank.coph.disp, rank.mev, fn.ws, gct.comb, z_score, file=fn.ws)
    
    
    ## #######################################################
    ##         generate some plots
    if(opt$nmf_only == FALSE){ 
      p <- try(nmf.post.processing(fn.ws, 
                                   core_membership=opt$core_membership, 
                                   organism=opt$organism,
                                   blank.anno=opt$blank_anno,
                                   blank.anno.col=opt$blank_anno_col,
                                   feature.fdr=opt$feature_fdr,
                                   ora.pval=opt$ora_pval,
                                   ora.max.categories=opt$ora_max_categories))
    }
    
    return(0)
} ## end main()
    
   
## #########################################################
##               run pipeline
main(opt)


