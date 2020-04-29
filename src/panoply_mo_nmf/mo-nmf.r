#!/usr/bin/env Rscript
options( warn = -1, stringsAsFactors = F )
#library(optparse)
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
  make_option( c("-c", "--classvariable"), action='store', type='character', dest='class.variable', help='Class variable of interest, NMF results will be mapped to levels of that variable. Must be present in the cdesc object of the GCT 1.3 file.'),
  make_option( c("-d", "--definecolors"), action='store', type='character', dest='class.colors', help='Specifiy colors for each level in "class variable". Example: "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"', default = "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"),
  make_option( c("-o", "--othervariable"), action='store', type='character', dest='variable.other', help='Other variable(s) in cdesc that should be included in the annotation column tracks of heatmaps. If multiple variables are of interest, separate by semicolon, e.g. -o "var1;var2;varN".', default=FALSE),
  make_option( c("-r", "--runonly"), action='store', type='logical', dest='no.plots', help='If TRUE no plots will be generated and the R-object after NMF clustering will be returned.', default=FALSE),
 # make_option( c("-g", "--genesofinterest"), action='store', type='character', dest='genes.of.interest', help='Character vector of gene symbols separated by semicolon. Example: "^PTPN1$;ERBB2;CDK5"', default=''),
  make_option( c("-f", "--sdfilter"), action='store', type='numeric', dest='sd.perc.rm', help='Lowest standard deviation across columns percentile to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest sd will be removed. Will be applied before z-scoring (optional)', default=0),
  make_option( c("-u", "--z_score"), action='store', type='logical', dest='zscore.all', help='If TRUE, rows in th matrix will be z-scored.', default=TRUE),
  make_option( c("-i", "--impute"), action='store', type='logical', dest='impute', help='If TRUE, the matrix will be first filtered for features present in >70% of samples. The remaining missing values will be imputed via KNN (R packe impute).', default=TRUE),
  make_option( c("-a", "--geneColumn"), action='store', type='character', dest='gene.column', help='Column name in rdesc in the GCT that contains gene names.', default='geneSymbol'),
  make_option( c("-e", "--exclude_2"), action="store", dest='exclude.2', type="logical", help="If TRUE, a number of clusters of will be excluded from calculation of the optimal rank.", default=TRUE),
  make_option( c("-z", "--libdir"), action="store", dest='lib.dir', type="character", help="the src directory.", default='/home/pgdac/src/')
 
  
  )
# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

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

#####################
## for local testing
dummy <- F
if(dummy){
    opt <- c()
    #opt$tar.file <- 'prospBC_prot_data_freeze_v2_1.tar'
    #opt$tar.file <- 'prospBC_psty_data_freeze_v2_1.tar'
    #opt$tar.file <- 'prospBC_rna_prot_psty_data_freeze_v2_1.tar'
    #opt$tar.file <- '../luad-v2.0-prot.tar'
    #opt$tar.file <- 'brca-combined-prot.tar'
    #opt$tar.file <- "data/prosp-brca-v3.0-rna-tumor.tar"
    opt$tar.file <- 'data/luad-v2.1-interim-prot-psty-rna.tar'
    #opt$tar.file <- "'data/prosp-brca-v3.0-prot-psty-rna-tumor.tar'"
    #opt$variable.other <- 'ER;PR;HER2.status.Satpathy;CDH1.mutation.status;GATA3.mutation.status;MAP3K1.mutation.status;PIK3CA.mutation.status;PTEN.mutation.status;TP53.mutation.status;ESTIMATE.ImmuneScore;ESTIMATE.StromalScore;Stemness.Score'
    opt$variable.other <- 'Region.of.Origin;Country.of.Origin;Gender;Smoking.Status;TP53;KRAS;STK11;EGFR;KEAP1;RB1;IL21R;EGFL6;LMO2;C10orf62;DKK3;fusion.EML4.ALK'
    #opt$class.colors <- "PTPN12Class=high:red;low:blue;NotAvailable:white|Group=ExtraCluster:darkred;NonRespondingCluster:darkgreen;RespondingCluster:orange" 
    # "PAM50=Basal:red;Her2:magenta;LumA:blue;LumB:cyan;Normal:grey|ER=positive:black;negative:white;unknown:grey|PR=positive:black;negative:white;unknown:grey|HER2=positive:black;negative:white;unknown:grey"
    #opt$class.colors <- "PAM50=Basal:red;Her2:magenta;LumA:blue;LumB:cyan;Normal:grey|ER=positive:black;negative:white;unknown:grey|PR=positive:black;negative:white;unknown:grey|HER2.status.Satpathy=positive:black;negative:white;unknown:grey|TP53.mutation.status=0:white;1:darkblue|PTEN.mutation.status=0:white;1:darkblue|PIK3CA.mutation.status=0:white;1:darkblue|MAP3K1.mutation.status=0:white;1:darkblue|GATA3.mutation.status=0:white;1:darkblue|CDH1.mutation.status=0:white;1:darkblue"
    opt$class.colors <- 'Gender=male:lightblue;female:pink|Country.of.Origin=usa:blue;vietnam:yellow;china:red;other:grey;ukraine:orange;bulgaria:darkgreen;poland:azure;russia:black|Stage=1:grey70;1A:grey60;1B:grey50;2A:grey40;2B:grey30;3:grey20;3A:grey10|Region.of.Origin=western:red;asian:green|Smoking.Status=smoker:darkred;non-smoker:darkgreen|TP53=0:white;1:black|KRAS=0:white;1:black|EGFR=0:white;1:black|STK11=0:white;1:black|KEAP1=0:white;1:black|RB1=0:white;1:black|IL21R=0:white;1:black|EGFL6=0:white;1:black|LMO2=0:white;1:black|C10orf62=0:white;1:black|DKK3=0:white;1:black|BIRC6=0:white;1:black|fusion.EML4.ALK=0:white;1:darkblue'
    #opt$class.variable <- "PTPN12Class" #"PAM50"
    opt$class.variable <- "PAM50"
    #opt$class.variable <- "Stage"
    opt$nrun <- 10
    opt$rank.low <- 2
    opt$rank.high <- 6
    opt$seed <- "random"
   # opt$genes.of.interest <- 'PTPN;CDK'
    opt$lib.dir <- 'C:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_mo_nmf/'
    opt$sd.perc.rm <- 0.05
    opt$zscore.all <- T
    opt$gene.column <- 'geneSymbol'
    opt$no.plots <- F
    opt$impute <- T
    opt$bnmf <- F
}

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
    
    ## apply z-score to all data types?
    zscore.all = opt$zscore.all
    
    ## #############################
    ## imputation parameters (KNN)
    impute <- opt$impute
    impute.k <- 5
    na.max.row <- 0.3
    na.max.col <- 0.9
     
    ## ################################
    ##            NMF parameters
    ranks <- opt$rank.low:opt$rank.high
    seed <- opt$seed
    nrun <- opt$nrun
    nrun.bs <- opt$nrun.bs ## bootstrap iterations
    
    ##method <- 'brunet'
    method <- ifelse(opt$bnmf, 'bnmf', 'nmf')
      
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
    
    expr <- data.imported$expr
    cdesc <- data.imported$cdesc
    rdesc <- data.imported$rdesc
    data.str <- data.imported$data.str
    
    ## ###############################
    ##     prepare output folder
    label <- paste(names(data.str), collapse='_')
    
    date.str <- glue("{sub(' .*', '', Sys.time())}_{label}_{nrun}_{method}{ifelse(impute,'_knn_','_')}sd_{sd.perc.rm}{ifelse(zscore.all,'_zcore','')}")
    dir.create(date.str)
    setwd(date.str)
    
    ## ###############################
    ##     generate parameter file
    param <- c('## NMF parameters:', 
               #paste('Bayesian NMF:', opt$bnmf), 
               paste('method:', method), 
               
               paste('rank:', paste(ranks, collapse=',')), 
               paste('seed:',seed), paste('nrun:', nrun), 
              
               paste('opts:', opts),
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
    sd.expr <- apply(expr, 1, sd, na.rm=T)
    
    ## ###################################################
    ## Z-score
    ## ###################################################
    if(zscore.all)
      expr <- scale(expr)
    
   
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
      
        cat("\napplying SD-filter: ")
        sd.perc <- quantile(sd.expr, c(sd.perc.rm))
        sd.keep <- which( sd.expr > sd.perc)
        cat("removed", nrow(expr)-length(sd.keep), "features with SD<=", round(sd.perc, 2), "(", names(sd.perc),"-tile)\n\n")
        expr <- expr[sd.keep,]
        
        ## add line to plot
        abline(v=sd.perc, lty='dashed', lwd=2, col='grey')
        text(sd.perc, ymax, paste(100*sd.perc.rm, 'th percentile', sep=''), pos=4, col='grey')
        
        ## export filtered GCT
        gct.filt <- gct.comb
        gct.filt@mat <- data.matrix(expr)
        gct.filt@rdesc <- rdesc[sd.keep, ]
        gct.filt@rid <- rownames(gct.filt@mat)
        #write.gct(gct.filt, ofile = 'mo-data-matrix-sd-filt')
        
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
        pdf('1_cluster_metrics.pdf')
        
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
      
      ## ####################################################
      ##      Estimate factorization rank
      ## ####################################################
      
      ## extract top N ranks based on cophentic correlation
      topn.rank <- 1
      #exclude.2 <- T ## should rank=2 be excluded?
      
      rank.top <- GetBestRank(rank.coph, topn.rank=topn.rank, exclude.2=opt$exclude.2)
      rank.mev  <- NULL
      
      #########################################
      ## plot cluster metrics
      if(length(ranks) > 1){
        
        col.tmp <- rep('darkblue', length(rank.coph))
        names(col.tmp) <- names(rank.coph)
        col.tmp[rank.top] <- 'red'
        
        
        ## plot
        pdf('1_cluster_metrics.pdf')
        ## SVD
        plot(rank.svd$d, cex=2, type='b', pch=20, main='SVD', xlab= 'Sample index', ylab='Singular values of x', col='darkblue')
        
        ## cophenetic
        plot(as.numeric(names(rank.coph)), rank.coph, type='b', pch=20, cex=2, col=col.tmp, xaxt='n', xlab='Factorization rank / No. of clusters', ylab='Cophenetic score', main='Cophenetic correlation coefficient', 
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


