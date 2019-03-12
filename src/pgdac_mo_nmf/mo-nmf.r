#!/usr/bin/env Rscript
options( warn = -1, stringsAsFactors = F )
#library(optparse)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("pacman"))

#Rscript c:\Users\karsten\Dropbox\Devel\PGDAC\src\pgdac_mo_nmf\mo-nmf.r -t ..\data\v2\ssGSEA\westbrook_baseline_ssGSEA_rna_psty_prot.tar -c Response -o Group -l 2 -m 3 -n 2
#Rscript c:\Users\karsten\Dropbox\Devel\PGDAC\src\pgdac_mo_nmf\mo-nmf.r -t .. -c PTPN12Class -o PTPN12IHCscores;Group -l 3 -m 3 -n 2 -b 10

# parse the directory this file is located in
#this.file.dir <- commandArgs()[4]
#this.file.dir <- '/home/pgdac/'

#this.file.dir <- 'c:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_mo_nmf/'
#this.file.dir <- sub('^(.*(/|\\\\)).*', '\\1', sub('.*?\\=','', this.file.dir))

#cat('------:',this.file.dir, '\n')
# source(paste(this.file.dir,  'my_plots.r', sep='/'))
# source(paste(this.file.dir,  'my_io.r', sep='/'))
#source(paste(this.file.dir,  'consensus_clustering.R', sep='/'))

#source(paste(this.file.dir, 'src', 'my_plots.r', sep='/'))
#source(paste(this.file.dir, 'src', 'my_io.r', sep='/'))

# specify command line arguments
option_list <- list(
  make_option( c("-t", "--tar"), action='store', type='character',  dest='tar.file', help='tar file containing data tables in GCT v1.3 format.'),
  make_option( c("-l", "--lowrank"), action='store', type='numeric',  dest='rank.low', help='Minimal factorization rank.', default = 2),
  make_option( c("-m", "--maxrank"), action='store', type='numeric',  dest='rank.high', help='Maximal factorization rank.', default = 4),
  make_option( c("-n", "--nrun"), action='store', type='numeric',  dest='nrun', help='Number of NMF runs with different starting seeds.', default = 5),
  #make_option( c("-b", "--bootstrap"), action='store', type='numeric',  dest='nrun.bs', help='Number of bootstrap iterations for consensus clustering.', default = 100),
  make_option( c("-s", "--seed"), action='store', type='character',  dest='seed', help='Seed method for NMF factorization.', default = 'random'),
  make_option( c("-c", "--classvariable"), action='store', type='character', dest='class.variable', help='Class variable of interest, NMF results will be mapped to levels of that variable. Must be present in the cdesc object of the GCT 1.3 file.'),
  make_option( c("-d", "--definecolors"), action='store', type='character', dest='class.colors', help='Specifiy colors for each level in "class variable". Example: "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"', default = "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"),
  make_option( c("-o", "--othervariable"), action='store', type='character', dest='variable.other', help='Other variable(s) in cdesc that should be included in the annotation column tracks of heatmaps. If multiple variables are of interest, separate by semicolon, e.g. -o "var1;var2;varN".', default=FALSE),
  make_option( c("-r", "--runonly"), action='store', type='logical', dest='no.plots', help='If TRUE no plots will be generated and the R-object after NMF clustering will be returned.', default=FALSE),
  make_option( c("-g", "--genesofinterest"), action='store', type='character', dest='genes.of.interest', help='Character vector of gene symbols separated by semicolon. Example: "^PTPN1$;ERBB2;CDK5"', default=''),
  make_option( c("-f", "--sdfilter"), action='store', type='numeric', dest='sd.perc.rm', help='Lowest standard deviation across columns percentile to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest sd will be removed. Will be applied before z-scoring (optional)', default=0),
  make_option( c("-b", "--z_score"), action='store', type='logical', dest='zscore.all', help='If TRUE, rows in th matrix will be z-scored.', default=TRUE),
  make_option( c("-i", "--impute"), action='store', type='logical', dest='impute', help='If TRUE, the matrix will be first filtered for features present in >70% of samples. The remaining missing values will be imputed via KNN (R packe impute).', default=TRUE),
  make_option( c("-a", "--geneColumn"), action='store', type='character', dest='gene.column', help='Column name in rdesc in the GCT that contains gene names.', default='geneSymbol'),
  make_option( c("-z", "--libdir"), action="store", dest='lib.dir', type="character", help="the src directory.", default='/home/pgdac/src/')
  
  )
# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

#save(opt, file='opt.RData')
## libraries
require('pacman')
p_load(NMF)
p_load(org.Hs.eg.db)
p_load(pheatmap)
p_load(WriteXLS)
p_load(RColorBrewer)
p_load(rhdf5)
if(!require(cmapR))devtools::install_github("cmap/cmapR")
p_load(dplyr)
p_load(Rtsne)
p_load(plotly)
p_load(rmarkdown)
p_load(maptools)
p_load(UpSetR)
p_load(ComplexHeatmap)
p_load(glue)
p_load(circlize)
p_load(impute)


#####################
## for local testing
dummy <- F
if(dummy){
    opt <- c()
    opt$tar.file <- 'data/prospBC_prot_phgcr_rna_hallmark.tar'
    opt$variable.other <- 'ER;PR;HER2'
    #opt$class.colors <- "PTPN12Class=high:red;low:blue;NotAvailable:white|Group=ExtraCluster:darkred;NonRespondingCluster:darkgreen;RespondingCluster:orange" 
    # "PAM50=Basal:red;Her2:magenta;LumA:blue;LumB:cyan;Normal:grey|ER=positive:black;negative:white;unknown:grey|PR=positive:black;negative:white;unknown:grey|HER2=positive:black;negative:white;unknown:grey"
    opt$class.colors <- "PAM50.ext=Basal:red;Her2:magenta;LumA:blue;LumB:cyan;Normal:grey|ER=positive:black;negative:white;unknown:grey|PR=positive:black;negative:white;unknown:grey|HER2=positive:black;negative:white;unknown:grey"
    #opt$class.variable <- "PTPN12Class" #"PAM50"
    opt$class.variable <- "PAM50.ext"
    opt$nrun <- 10
    opt$rank.low <- 2
    opt$rank.high <- 6
    opt$seed <- "random"
    opt$genes.of.interest <- 'PTPN;CDK'
    opt$lib.dir <- 'C:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_mo_nmf/'
    opt$sd.perc.rm <- 0.1
    opt$zscore.all <- T
    opt$gene.column <- 'geneSymbol'
    opt$no.plots <- F
   
}

## source required files
#cat('test1:',file.path(opt$lib.dir, 'my_plots.r'))
source(file.path(opt$lib.dir,  'my_plots.r'))
#cat('test2')
source(file.path(opt$lib.dir,  'my_io.r'))
#cat('test3')
source(file.path(opt$lib.dir, 'nmf_functions.R'))
#cat('test4')

## #############################################
## main function
main <- function(opt){

    ## ##################################################
    ##           parse parameters
    ## ###################################################
    tar.file <- opt$tar.file
  
    ## percentile to remove from the data
    sd.perc.rm <- opt$sd.perc.rm
    
    ## apply z-score to all data types
    zscore.all = opt$zscore.all
    
    ## ########################
    ## imputation (KNN)
    impute <- opt$impute
    impute.k <- 5
    na.max.row <- 0.7
    na.max.col <- 0.8
    # 
    # ## ###################################################
    # ## class variable/annotation tracks of interest
    # class.variable <- opt$class.variable
    # variable.other <- strsplit( opt$variable.other, ';' ) %>% unlist
    # class.variable <- gsub('\\\'', '', class.variable) 
    # variable.other <- gsub('\\\'', '', variable.other) 
    # 
    ## ################################
    ##            NMF paramters
    ranks <- opt$rank.low:opt$rank.high
    seed <- opt$seed
    nrun <- opt$nrun
    nrun.bs <- opt$nrun.bs ## bootstrap iterations
    method <- 'brunet'
    cores <- detectCores()
    opts <- paste('vp', cores,'t', sep='')
  
          
    ## ###################################################
    ##      other parmaters
    ## ###################################################
    tmp.dir <- 'tmp'
    
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
    #data.imported <- import.data.sets(data.str, zscore.cnv)
    expr <- data.imported$expr
    cdesc <- data.imported$cdesc
    rdesc <- data.imported$rdesc
    data.str <- data.imported$data.str
    
    ## ###############################
    ##     prepare output folder
    label <- paste(names(data.str), collapse='_')
    
    #date.str <- paste(paste( sub(" .*", "", Sys.time()), label ,sep='_'), nrun, method, 'Zscore', ifelse(zscore.all, 'T', 'F'), 'SD', sd.perc.rm, sep='-')
    date.str <- glue("{sub(' .*', '', Sys.time())}_{label}_{nrun}_{method}{ifelse(impute,'_knn_','_')}sd_{sd.perc.rm}{ifelse(zscore.all,'_zcore','')}")
    #paste(paste( sub(" .*", "", Sys.time()), label ,sep='_'), nrun, method, 'Zscore', ifelse(zscore.all, 'T', 'F'), 'SD', sd.perc.rm, sep='-')
    dir.create(date.str)
    setwd(date.str)
    
    ## ###############################
    ##     generate parameter file
    param <- c('## NMF parameters:', 
               paste('rank:', paste(ranks, collapse=',')), 
               paste('seed:',seed), paste('nrun:', nrun), 
               paste('method:', method), paste('opts:', opts),
               '',
               '## data filter / normalization / imputation', 
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
    #expr.full.org <- expr
    if(impute){
       keep.idx <- which(apply(expr, 1, function(x) ifelse( sum(is.na(x)) > na.max.row, F, T)))
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
    
    # ## ############################################
    # ## save a copy
    # expr.org <- expr
    # cdesc.org <- cdesc
    # rdesc.org <- rdesc
    # variable.other.org <- variable.other
    # class.variable.org <- class.variable
    # 
    # 
    ## ############################################
    ##  export merged data used for nmf as gct
    gct.comb <- new('GCT')
    gct.comb@mat <- data.matrix(expr)
    gct.comb@cdesc <- cdesc
    gct.comb@rdesc <- rdesc
    gct.comb@cid <- colnames(gct.comb@mat)
    gct.comb@rid <- rownames(gct.comb@mat)
    write.gct(gct.comb, ofile = 'mo-data_matrix')
    
    
    ## ###################################################
    ## Z-score
    ## ###################################################
    if(zscore.all)
      expr <- scale(expr)
    
    ## ###################################################
    ##       calculate SD accross samples
    ## ###################################################
    sd.expr <- apply(expr, 1, sd, na.rm=T)
    
 
    ## #########################
    ## density kernel
    sd.d <- density( sd.expr )
    sd.d.list <- list()
    ymax <- 0
    for(i in 1:length(data.str)){
        sd.d.list[[i]] <- density( sd.expr[ grep(names(data.str)[i], names(sd.expr)) ] )
        if(max(sd.d.list[[i]]$y) > ymax)
            ymax <- max(sd.d.list[[i]]$y)
    }
    ymax <- (ymax+(0.1*ymax))
    
    ## #####################
    ## plot
    pdf('0_density_StdDev.pdf')
    plot(sd.d, ylim=c(0, ymax), main='StdDev accross samples', lwd=2)
    for(i in 1:length(data.str))
        lines(sd.d.list[[i]], col=i+1, lwd=2)
    legend('topright', legend=c('Aggregated', names(data.str)), lty='solid', lwd=3, col=1:(length(data.str)+1))
    
    ## ####################################################
    ## sd filter
    if(sd.perc.rm > 0){
        sd.perc <- quantile(sd.expr, c(sd.perc.rm))
        sd.keep <- which( sd.expr > sd.perc)
        expr <- expr[sd.keep,]
    
        ## add line to plot
        abline(v=sd.perc, lty='dashed', lwd=2, col='grey')
        text(sd.perc, ymax, paste(100*sd.perc.rm, 'th percentile', sep=''), pos=4, col='grey')
    }
    dev.off()
    
    
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
    
    ## loop over ranks
    for(i in 1:length(ranks))
          res.rank[[i]] <-  nmf(expr.comb, rank=ranks[i], method=method, seed=seed, nrun=nrun, .options = opts)
      
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

    if(length(ranks) > 1){
      
      ## plot
      pdf('1_cluster_metrics.pdf')
      ## SVD
      plot(rank.svd$d, cex=2, type='b', pch=20, main='SVD', xlab= 'Sample index', ylab='Singular values of x', col='darkblue')
 
      ## cophenetic
      plot(as.numeric(names(rank.coph)), rank.coph, type='b', pch=20, cex=2, col='darkblue', xaxt='n', xlab='Factorization rank / No. of clusters', ylab='Cophenetic score', main='Cophenetic correlation coefficient', 
           ylim=c( min(rank.coph)-(min(rank.coph)*0.1), 1) )
      axis(1, at=as.numeric(names(rank.coph)))
       
      ## silhouette
      # par(mfrow=c(2,1))
      #fancyBoxplot( lapply(rank.sil, function(x)x[, 3] ), xlab='Factorization rank / No. of clusters',
      #              main='Silhouette scores', ylab='Silhouette scores', show.numb = 'mean', numb.cex = 1.5, numb.pos = 4, 
      #              col='white', box.border = 'darkblue', vio.alpha = 0, lwd=2.5)
      #fancyBoxplot( lapply(rank.sil.random, function(x)x[, 3] ), xlab='Factorization rank / No. of clusters', main='Randomized data', ylab='Silhouette scores', show.numb = 'mean', numb.cex = 1.5, numb.pos = 4)
      dev.off()
    }
    
    ## ####################################################
    ##      Estimate factorization rank
    ## ####################################################
    
    ## extract top N ranks based on cophentic correlation
    topn.rank <- 3
    exclude.2 <- T ## should rank=2 be excluded?
    
    rank.top <- names(rank.coph)[ order(rank.coph, decreasing=T)] [1:topn.rank]
    if(exclude.2){
      rank.top <- setdiff(rank.top, '2') 
    }
    
    ## #####################################################
    ## save results needed to 
    fn.ws <- 'workspace_after_NMF.RData'
    save(opt, res.rank, rank.top, expr.comb, expr, rdesc, cdesc, rank.sil, 
         rank.coph, fn.ws, gct.comb, file=fn.ws)
    
    
    ## #######################################################
    ##         generate some plots
    if(opt$no.plots == FALSE){ 
      p <- try(nmf.plots(fn.ws))
    }
    
    return(0)
} ## end main()
    
   
## #########################################################
##               run pipeline
main(opt)
#fn.ws <- 'workspace_after_NMF.RData'
#save.image(fn.ws)
#nmf.plots(fn.ws)

