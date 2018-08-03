#!/usr/bin/env Rscript
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("pacman"))

#Rscript c:\Users\karsten\Dropbox\Devel\PGDAC\src\pgdac_mo_nmf\mo-nmf.r -t ..\data\v2\ssGSEA\westbrook_baseline_ssGSEA_rna_psty_prot.tar -c Response -o Group -l 2 -m 3 -n 2

# parse the directory this file is located in
this.file.dir <- commandArgs()[4]
#this.file.dir <- '/home/pgdac/'

#this.file.dir <- 'c:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_mo_nmf/'
this.file.dir <- sub('^(.*(/|\\\\)).*', '\\1', sub('.*?\\=','', this.file.dir))

#cat('------:',this.file.dir, '\n')
source(paste(this.file.dir,  'my_plots.r', sep='/'))
source(paste(this.file.dir,  'my_io.r', sep='/'))
source(paste(this.file.dir,  'nmf_consensus.R', sep='/'))

#source(paste(this.file.dir, 'src', 'my_plots.r', sep='/'))
#source(paste(this.file.dir, 'src', 'my_io.r', sep='/'))


# specify command line arguments
option_list <- list(
  make_option( c("-t", "--tar"), action='store', type='character',  dest='tar.file', help='tar file containing data tables in GCT v1.3 format.'),
  make_option( c("-l", "--lowrank"), action='store', type='numeric',  dest='rank.low', help='Minimal factorization rank.', default = 2),
  make_option( c("-m", "--maxrank"), action='store', type='numeric',  dest='rank.high', help='Maximal factorization rank.', default = 4),
  make_option( c("-n", "--nrun"), action='store', type='numeric',  dest='nrun', help='Number of NMF runs with different starting seeds.', default = 5),
  make_option( c("-b", "--bootstrap"), action='store', type='numeric',  dest='nrun.bs', help='Number of bootstrap iterations for consensus clustering.', default = 100),
  make_option( c("-s", "--seed"), action='store', type='character',  dest='seed', help='Seed method for NMF factorization.', default = 'random'),
  make_option( c("-c", "--classvariable"), action='store', type='character', dest='class.variable', help='Class variable of interest, NMF results will be mapped to levels of that variable. Must be present in the cdesc object of the GCT 1.3 file.'),
  make_option( c("-d", "--definecolors"), action='store', type='character', dest='class.colors', help='Specifiy colors for each level in "class variable". Example: "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"', default = "MedResponder:darkred;NonResponder:darkgreen;Responder:orange"),
  make_option( c("-o", "--othervariable"), action='store', type='character', dest='variable.other', help='Other variable(s) in cdesc that should be included in the annotation column tracks of heatmaps. If multiple variables are of interest, separate by semicolon, e.g. -o "var1;var2;varN".', default=FALSE),
  make_option( c("-r", "--runonly"), action='store', type='logical', dest='no.plots', help='If TRUE no plots will be generated and the R-object after NMF clustering will be returned.', default=FALSE),
  make_option( c("-g", "--genesofinterest"), action='store', type='character', dest='genes.of.interest', help='Character vector of gene symbols separated by semicolon. Example: "^PTPN1$;ERBB2;CDK5"', default='')
  
  )
# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

#####################
## for local testing
dummy <- F
if(dummy){
    opt$tar.file <- 'bc_prosp_pr_ph_rna-2comp-cnv_harmon.tar'
    opt$variable.other <- 'ER;PR;HER2'
    opt$class.colors <- "PAM50=Basal:red;Her2:magenta;LumA
    :blue;LumB:cyan;Normal:grey|ER=positive:black;negative:white;unknown:grey|PR=positive:black;negative:white;unknown:grey|HER2=positive:black;negative:white;unknown:grey"
    opt$class.variable <- "PAM50"
    opt$nrun <- 2
    opt$nrun.bs <- 10
    opt$rank.low <- 4
    opt$rank.high <- 4
    opt$seed <- "random"
    opt$genes.of.interest <- ''
}

tmp.dir ='tmp'

## #################################################
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

## #################################
## extract tar ball
if(!dir.exists(tmp.dir))
  dir.create(tmp.dir)
untar(opt$tar.file, exdir=tmp.dir)

##  import config file
conf <- read.delim(paste( tmp.dir, 'nmf.conf', sep='/'), row.names = NULL, stringsAsFactors = F, header=F)
data.str <- paste('..', tmp.dir, conf[, 2], sep='/')
names(data.str) <- conf[,1]


## ###################################################
## data filtering / normalization
## ###################################################
zscore.cnv = F     ## should CNVs be z-scored?
zscore.all = F     ## apply z-score to all data types
sd.perc.rm = 0     ## percentile to remove from the data

# variables of interest
class.variable <- opt$class.variable
variable.other <- strsplit( opt$variable.other, ';' ) %>% unlist

gene.column <- 'geneSymbol'


class.variable <- gsub('\\\'', '', class.variable) 
variable.other <- gsub('\\\'', '', variable.other) 

## ################################################
##            NMF paramters
## ###################################################
ranks <- opt$rank.low:opt$rank.high
seed <- opt$seed
nrun <- opt$nrun
nrun.bs <- opt$nrun.bs ## bootstrap iterations
method <- 'brunet'
cores <- detectCores()
opts <- paste('vp', cores,'t', sep='')

## ####################################################
## pheatmap
cw <- 15
ch <- 15
max.val <- 10


## ####################################################
## genes of interest
genes.of.interest <- opt$genes.of.interest
if(nchar(genes.of.interest) < 3){
    genes.of.interest <- NULL
} else{
  genes.of.interest <- strsplit( genes.of.interest, ';' ) %>% unlist
}
  #genes.of.interest <- c('^PTPN', 'ERBB2', 'CDK5', 'CDK12', 'PAK1', 'PTK2', 'RIPK2', 'TLK2', 'CETN3', 'SKP1', 'TP53', 'GATA3', 'ESR1', 'PGR')


## ##################################################
## colors for class vector
class.colors <- opt$class.colors
color.all <- class.colors %>% strsplit(., '\\|') %>% unlist

cdesc.color <- vector('list', length(color.all))
names(cdesc.color) <- sub('^(.*?)=.*', '\\1', color.all)

color.all <- sub('.*?=','', color.all)

for(i in 1:length(cdesc.color)){
  color.tmp <- color.all[i] %>% strsplit(., ';') %>% unlist %>% strsplit(. , ':')
  color.names <- sapply(color.tmp, function(x)x[1])
  color.tmp <- sapply(color.tmp, function(x)x[2])
  names(color.tmp) <- color.names
  color.tmp <- gsub('\\\'', '', color.tmp) 
  names(color.tmp) <- gsub('\\\'', '', names(color.tmp)) 
  cdesc.color[[i]] <- color.tmp
}
cdesc.color <- lapply(cdesc.color, function(x){names(x)=sub('\\n.*', '',names(x));x})
cdesc.color.org <- cdesc.color


## #########################################################################################
##
##                               START
##
## #########################################################################################

## ###############################
##     prepare output folder
label <- paste(names(data.str), collapse='_')

date.str <- paste(paste( sub(" .*", "", Sys.time()), label ,sep='_'), nrun, method, 'Zscore', ifelse(zscore.all, 'T', 'F'), 'SD', sd.perc.rm, sep='-')
dir.create(date.str)
setwd(date.str)

## ###############################
##     generate parameter file
param <- c('## NMF parameters:', paste('rank:', paste(ranks, collapse=',')), paste('seed:',seed), paste('nrun:', nrun), paste('method', method), paste('opts:', opts),'', '## data filter / normalization', paste('z-score ALL:', zscore.all), paste('remove lowest SD perc.:', sd.perc.rm), '',paste('## data tables'), paste(data.str) )
writeLines(param, con='parameters.txt')


## ###################################################
## import data tables
for(i in 1:length(data.str)){

    ##gct
    gct <-  parse.gctx2(data.str[i])  
    ## extract data matrix
    data.tmp <- gct@mat
    ## append data type to row names
    rownames(data.tmp) <- paste(names(data.str)[i], gct@rid, sep='-')
    colnames(data.tmp) <- gct@cid

    ## remove columns with no data values
    keep.idx <- which(apply(data.tmp, 2, function(x) ifelse(sum(is.na(x))/length(x) == 1, F, T)))
    data.tmp <- data.tmp[, keep.idx]
    
    ## row annotations    
    rdesc.tmp <- gct@rdesc
    rownames(rdesc.tmp) <- paste(names(data.str)[i], gct@rid, sep='-')
    rdesc.tmp <- data.frame( Data.Type.ID=rownames(rdesc.tmp), rdesc.tmp)
    
    ## column annotations
    cdesc.tmp <- gct@cdesc %>% data.frame
    
    #cdesc.tmp[] <- lapply(cdesc.tmp, as.character) %>% as.data.frame()
    if(nrow(cdesc.tmp) > 0){
      cdesc.tmp <- data.frame(ID=gct@cid, cdesc.tmp)
      cdesc.tmp <- cdesc.tmp[ keep.idx,  ]
    }
    
    ## Z-score CNV?
    if(names(data.str)[i] == 'CNV' & zscore.cnv)
        data.tmp <- apply(data.tmp, 2, function(x) (x-median(x, na.rm=T))/mad(x, na.rm=T))
    
    ## merge data tables, row and column annotations
    if(i == 1){
        data <- data.tmp
        samp <- colnames(data.tmp)
        rdesc <- rdesc.tmp
        cdesc <- cdesc.tmp
      
    } else {
        ## only sample columns present in all tables will be carried over
        samp <- intersect(samp, colnames(data.tmp))
        data <- data[, samp]
        cdesc <- cdesc[samp, ]
        cdesc[] <- lapply(cdesc, function(x){if(sum(is.na(as.numeric(x))) == 0){as.numeric(x);} else {x} })

        ## append to data matrix
        data <- rbind(data, data.tmp[, colnames(data)])
        
        ## row annotations
        rdesc <- full_join(rdesc, rdesc.tmp)
        ## add back rownames
        rownames(rdesc) <- rdesc$Data.Type.ID
        ## ensure same order as data
        rdesc <- rdesc[rownames(data), ]
        
        ## column annotations
        if(nrow(cdesc.tmp) > 0){
          cdesc.tmp[] <- lapply(cdesc.tmp, function(x){if(sum(is.na(as.numeric(x))) == 0){as.numeric(x);} else {x} })
          #cdesc.tmp[] <- lapply(cdesc.tmp, as.character) %>% as.data.frame()
          cdesc <- left_join(cdesc, cdesc.tmp)
          rownames(cdesc) <- cdesc$ID
          #cdesc2 <- full_join(data.frame(t(cdesc)), data.frame(t(cdesc.tmp)))
        }
      }
}


## ##################################################
## extract expression data
expr <- data.matrix(data)
expr.org <- expr
cdesc.org <- cdesc
variable.other.org <- variable.other

## ##################################################
## fully quantified
keep.idx <- which(apply(expr, 1, function(x) ifelse(sum(is.na(x)) > 0, F, T)))
expr <- expr[keep.idx, ]
expr.full.org <- expr

## ##################################################
## Z-score
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
expr.comb <- make.nn(expr)

## remove rows only containing zero values
keep.idx <- which(apply(expr.comb, 1, function(x) sum(x != 0) ) > 0)
expr.comb <- expr.comb[keep.idx, ]


## ###################################################
##     Estimate factorization rank
## ###################################################
res.rank <- res.rank.random <- vector('list', length(ranks))
names(res.rank) <- names(res.rank.random) <- ranks
res.cons <- res.cons.random <- res.rank

## shuffle data
expr.comb.random <- randomize(expr.comb)

## loop over clusters
for(i in 1:length(ranks)){
  
  ###############################
  # actual data
  res.rank[[i]] <-  nmf(expr.comb, rank=ranks[i], method=method, seed=seed, nrun=nrun, .options = opts)
  ## consensus matrix
  res.cons[[i]] <- consensus.clust( expr.comb, method='nmf', nmf.nrun = 10, bs.nrun = nrun.bs, nmf.opt = opts )
  
  #################################
  # shuffled data
  res.rank.random[[i]] <- nmf(expr.comb.random, rank=ranks[i], method=method, seed=seed, nrun=nrun, .options = opts)
  ## consensus matrix
  res.cons.random[[i]] <- consensus.clust( expr.comb.random, method='nmf', nmf.nrun = 10, bs.nrun = nrun.bs, nmf.opt = opts )
 
}

## ####################################################
## perform actual consensus NMF clustering
## ####################################################

## save results
save.image('workspace_after_NMF.RData')


## #######################################################
##         generate some plots
if(opt$no.plots == FALSE){ 
  
    ## ###################################################
    ##       calculate cluster metrics
    ## ####################################################
    
    ## SVD
    rank.svd <- svd(expr.comb)
    rank.svd.random <- svd(expr.comb.random)
    ## silhouette
    rank.sil <- lapply(res.rank, silhouette)
    rank.sil.random <- lapply(res.rank.random, silhouette)
    rank.sil.avg <- lapply(rank.sil, function(x) tapply(x[,3], x[, 1], mean)) 
    ## cophenic
    rank.coph <- sapply(res.rank, cophcor)
    rank.coph.random <- sapply(res.rank.random, cophcor)
    
    if(length(ranks) > 1){
      
          ## plot
          pdf('1_cluster_metrics.pdf')
          ## SVD
          plot(rank.svd$d, cex=2, type='b', pch=20, main='SVD', xlab= 'Sample index', ylab='Singular values of x', col='darkblue')
          lines(rank.svd.random$d, cex=2, type='b', col='grey', pch=20)
          legend('right', legend=c('Actual data', 'Randomized data'), col=c('darkblue', 'grey'), pch=20, cex=1.5, bty='n')
          
          ## cophenetic
          plot(as.numeric(names(rank.coph)), rank.coph, type='b', pch=20, cex=2, col='darkblue', xaxt='n', xlab='Factorization rank / No. of clusters', ylab='Cophenetic score', main='Cophenetic correlation coefficient', ylim=c(min(c(rank.coph.random, rank.coph)), 1) )
          axis(1, at=as.numeric(names(rank.coph)))
          lines(as.numeric(names(rank.coph.random)), rank.coph.random,  type='b', pch=20, cex=2, col='grey')
          legend('right', legend=c('Actual data', 'Randomized data'), col=c('darkblue', 'grey'), pch=20, cex=1.5, bty='n')
          
          ## silhouette
          par(mfrow=c(2,1))
          fancyBoxplot( lapply(rank.sil, function(x)x[, 3] ), xlab='Factorization rank / No. of clusters', main='Silhouette scores', ylab='Silhouette scores', show.numb = 'mean', numb.cex = 1.5, numb.pos = 4, col='darkblue')
          fancyBoxplot( lapply(rank.sil.random, function(x)x[, 3] ), xlab='Factorization rank / No. of clusters', main='Randomized data', ylab='Silhouette scores', show.numb = 'mean', numb.cex = 1.5, numb.pos = 4)
          dev.off()
    }
    
    
    
    ## #######################################################################
    ##
    ##              Generate some plots
    ##
    ## #######################################################################
    for(rank in as.character(ranks)){
      
      dir.create(paste('K', rank, sep='_'))
      setwd(paste('K', rank, sep='_'))
    
      
      ## extract NMF results
      res <- res.rank[[as.character( rank )]]
    
      
      ## ##########################################
      ## silhoutte plots
      pdf(paste('1.0_silhouette_K_', rank, '.pdf', sep=''), 10, 6)
      par(mfrow=c(1,2))
      plot(rank.sil[[rank]], main=paste('K=', rank, sep=''), col=palette()[1:as.numeric(rank)+1] )
      plot(rank.sil.random[[rank]], main=paste('Randomized data', sep=''), col=palette()[1:as.numeric(rank)+1])
      dev.off()
      
      
      ## colors
      cdesc.color <- cdesc.color.org
      cdesc <- cdesc.org
      variable.other <- variable.other.org
      
      
      ## ########################################
      ## add NMF classification to cdesc object
      NMF.basis <- predict(res)
      NMF.consensus <- predict(res, 'consensus')
      
      ## map consensuns to basis
      NMF.consensus.map <- NMF.consensus
      for(ii in 1:as.numeric(rank))
        NMF.consensus[ NMF.consensus == ii ] <-  NMF.basis[NMF.consensus == ii]
      
      ## add to cdesc
      cdesc <- data.frame(cdesc, NMF.basis=as.factor(NMF.basis), NMF.consensus=as.factor(NMF.consensus), stringsAsFactors=F)
      
      ## export
      write.table(cdesc, file=paste('clin_anno_nmf.txt'), sep='\t', quote=F, na='', col.names=NA)
        
      ## add to heatmap annotation tracks
      variable.other <- c('NMF.consensus', variable.other)
    
          
      ## #######################################################################################################
      ## add NMF to color list for pheatmap
      ## - map NMF classes to levels in 'class.variable'
      mb.nmf.map <- data.frame(mb=cdesc[, class.variable], nmf_basis=NMF.basis, nmf_cons=NMF.consensus)
        
      ## map consensus to 'class.variable'
      cons.map <- tapply(mb.nmf.map$mb, mb.nmf.map$nmf_cons, table)
      
      ## relative frequency
      cons.map.rel <- lapply(cons.map, function(x) x/sum(x))
      ## matrix
      cons.map.rel.mat <- matrix(0, nrow=as.numeric(rank), ncol=length(unique(cdesc[, class.variable])), dimnames = list(1:as.numeric(rank), unique(cdesc[, class.variable])))
      for(j in 1:nrow(cons.map.rel.mat))
        cons.map.rel.mat[j, names(cons.map.rel[[j]])] <- cons.map.rel[[j]]
      write.table(cons.map.rel.mat, col.names = NA, sep='\t', file=paste('nmf_vs_', class.variable, '.txt', sep=''), quote=F)
      
      ## pick the maximum
      #cons.map.max <- sapply(cons.map.rel, function(x) {m=max(x); ifelse( sum(x %in% m) > 1, paste(names(x)[which(x %in% m)] , collapse='|'), names(x)[which.max(x)])})
      cons.map.max <- sapply(cons.map.rel, function(x) {
        m=max(x); 
        ifelse( sum(x %in% m) > 1, 
                paste(names(x)[which(x %in% m)] , collapse='|'),
                paste(names(x)[x > 0.3], collapse='|' ) )
      })
    
      ## assign colors 
      NMF.consensus.col <- cdesc.color[[class.variable]][ cons.map.max ]
      names(NMF.consensus.col) <- cons.map.max
      
      ######################################################  
      ## check whether all NMF basis have color
      if( sum(is.na(NMF.consensus.col)) > 0 ){
        
        #idx.tmp <- which( !(levels(NMF.consensus) %in% names(NMF.consensus.col)) )
        idx.tmp <- which( is.na(NMF.consensus.col) )
        
        col.tmp <-  rev( palette() )[ (1:length(idx.tmp) ) +  1] 
        #names(col.tmp) <- idx.tmp
        #NMF.consensus.col <- c(NMF.consensus.col, col.tmp)
        NMF.consensus.col[idx.tmp] <- col.tmp
      }
      
      cdesc.color$NMF.consensus <- NMF.consensus.col
      names(cdesc.color$NMF.consensus) <- 1:length(NMF.consensus.col)
      

      ## ########################################
      ## heatmap coefficients
      #coefmap(res, annCol=cdesc[, c(class.variable, variable.other)], annColors=cdesc.color, filename='1_coefmap.pdf', tracks='consensus', main='Metagenes (consensus)')
      coefmap(res, annCol=cdesc[, c(class.variable, variable.other)], annColors=cdesc.color, filename='2.0_coefmap.pdf')
      
      # pheatmap version
      H <- res@fit@H
      H.norm <- apply(H, 2, function(x)x/max(x))
      pheatmap(H.norm, cluster_row=F, annotation_col=cdesc[, rev(c(class.variable, variable.other))],
                 annotation_colors = cdesc.color,
                 color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), 
                 filename='2.1_coefmap_pheatmap.pdf', cellwidth = cw, cellheight = ch)
      
      ## disable clustering
      pheatmap(H.norm[, order(cdesc[, 'NMF.consensus'])], cluster_row=F, cluster_cols = F, annotation_col=cdesc[, rev(c(class.variable, variable.other) )],
               annotation_colors = cdesc.color,
               color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), 
               filename='2.2_coefmap_pheatmap_sorted.pdf', cellwidth = cw, cellheight = ch)
      
        
        #pheatmap(H.norm, cluster_row=F, color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), filename='2_coefmap_pheatmap.pdf')
        
        
        ## #######################################
        ## consensus map
        cons <- res@consensus
        colnames(cons) <- rownames(cons) <- colnames(H)
        
        res2 <- res
        res2@consensus <- cons
        consensusmap(res2, filename=paste('3.0_consensusmap_nrun_', nrun, '.pdf', sep=''))
        
        cons.ord <- order(NMF.consensus)
        cons <- cons[cons.ord, cons.ord ]
        
        pheatmap(cons, cluster_row=F, cluster_col=F, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
                 annotation_colors=cdesc.color, color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), 
                 cellwidth = cw, cellheight = ch,
                 filename=paste('3.1_consensusmap_nrun_', nrun,'_pheatmap.pdf', sep=''))
       
        
        ## ##########################################
        ## consensus matrix from bootstrapping
        cons.bs <- res.cons[[as.character( rank )]]
        pheatmap(cons.bs, cluster_row=T, cluster_col=T, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
                 annotation_row = cdesc[ , rev(c(class.variable, variable.other))],
                 annotation_colors=cdesc.color, color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), 
                 cellwidth = cw-10, cellheight = ch-10, symm=T, fontsize_row = 5, fontsize_col = 5,
                 filename=paste('3.2_consensus_matrix_bootstrap_nrun_', nrun.bs,'.pdf', sep=''))
        
        
        cons.bs.rand <- res.cons.random[[as.character( rank )]]
        pheatmap(cons.bs.rand, cluster_row=T, cluster_col=T, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
                 annotation_colors=cdesc.color, color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), 
                 cellwidth = cw-10, cellheight = ch-10,
                 filename=paste('3.3_consensus_matrix_RANDOM_bootstrap_nrun_', nrun.bs,'.pdf', sep=''))
        
        
        
        ##V.hat = fitted(res)
        ##fit(res)
        
        ## ###############################################################
        ##
        ##                   extract NMF features
        ##
        ## ###############################################################
        feature.methods <- list('kim', 'max', 1, '10L')
        feature.methods.score <- list('kim', 'max', 'max', 'max')
        
        success <- F
        pdf('4.0_histogram_feature_scores.pdf')
        j <- 1
        while(!success){

          feature.method=feature.methods[[j]]
          feature.method.score=feature.methods.score[[j]]
          
          ## feature scores
          s.i <- featureScore(res, method=feature.method.score)
          
          hist( unlist(s.i), main=paste("Feature scores (method:", feature.method,")"), xlab='Feature score' )
      
          ## extract features
          s <- extractFeatures(res, method=feature.method)
        
          if( sum( sapply(s, function(x) sum(is.na(x))) ) == 0 ){
            success=T
          } else {
            warning(paste("\nExtracting NMF features using method '", feature.method, "' did not return features for all basis\n"))    
            j <- j + 1
            if(j > length(feature.methods)){
             break
            }  
          }
        }
        dev.off()
    
        ## ####################################
        ##  annotate the features
        
        ## get accession numbers
        #s.acc <- lapply(s, function(x) unique(sub('_up|_down','',rownames(expr.comb)[x])))
        s.acc <- lapply(s, function(x) data.frame(Accession=sub('_up|_down','',rownames(expr.comb)[x]) , id=rownames(expr.comb)[x], stringsAsFactors = F) ) 
        names(s.acc) <- paste('features_basis', 1:length(s.acc), sep='_')
    
        
        ## features plus scores
        ##s.acc.scores <- lapply(s, function(x) s.i[x])
        s.acc.scores <- lapply(s.acc, function(x) data.frame(x, Score=s.i[x$id], stringsAsFactors = F))
        
        
        
        ## upset plot
        upset.mat <- matrix(0, ncol=length(s.acc), nrow=length(unique(unlist(sapply(s.acc, function(x)x$Accession )))), dimnames=list(unique(unlist(sapply(s.acc, function(x)x$Accession ))), names(s.acc)))
        for(ii in names(s.acc))
          upset.mat[s.acc[[ii]]$Accession, ii] <- 1
        
        pdf('4.1_upset_NMF_markers.pdf')
        upset(data.frame(upset.mat), point.size = 4, text.scale = 1.5)
        dev.off()
        
        
        ## gene names
        if(gene.column %in% colnames(rdesc)){
          s.gn <- lapply(s.acc, function(x) rdesc[x$Accession, gene.column] )
          s.gn.red <- s.gn
          
          } else {
            s.gn <- lapply(s.acc, function(x) unique(sub('\\|.*','',sub('.*?\\|(.*).*', '\\1', sub('^CNV-(.*?)\\|.*', '\\1', x$Accession)))))
            ## gene names redundant
            s.gn.red <- lapply(s.acc, function(x) sub('\\|.*','',sub('.*?\\|(.*).*', '\\1', sub('^CNV-(.*?)\\|.*', '\\1', x$Accession))))
          }
        
        ## keep accession
        s.acc2 <- lapply(s.acc, function(x)sub('^(.*?)\\|.*$', '\\1', x$Accession))
        
        ## annotations as data frame
        s.acc.gn <- s.acc
        for(i in 1:length(s.acc)){
        
            tmp <- data.frame(
                Type=sub('-.*','', s.acc[[i]]$Accession),
                Accession=sub('.*?-','', s.acc2[[i]]),
                Direction=sub('.*_(up|down)$','\\1', s.acc[[i]]$id),
                SYMBOL=s.gn.red[[i]],
                Score=s.acc.scores[[i]]$Score,
                stringsAsFactors=F)
            s.acc.gn[[i]] <- tmp[order(tmp$SYMBOL),]
        }
        s.acc.gn.anno <- s.acc.gn
        
        ## ################################################################
        if(gene.column %in% colnames(rdesc)){
          ## add description and enzyme codes
            s.acc.gn.anno <- lapply( s.acc.gn, function(x)
              AnnotationDbi::select(org.Hs.eg.db, keys=x$SYMBOL , column=c( 'GENENAME',  'ENZYME'), keytype='SYMBOL', multiVals='first')
            )
        
          ## join in a single data frame
          for(i in 1:length(s.acc)){
            tmp <- merge(s.acc.gn[[i]], s.acc.gn.anno[[i]], 'SYMBOL' )
            tmp <- tmp[which(!duplicated(apply(tmp, 1, paste, collapse=' '))), ]
            s.acc.gn.anno[[i]] <- tmp
          }
            
            
            ## ##########################################################
            ## extract protein kinases and phosphatases
            s.acc.gn.kinase <- lapply(s.acc.gn.anno, function(x) x[ grep('2\\.7\\.1[0-2]|3\\.1\\.3\\.(16|48)', x$ENZYME), ])
            names(s.acc.gn.kinase) <- paste(names(s.acc.gn.kinase), 'KINASE', sep='_')
            
            s.acc.gn.kinase <- lapply(s.acc.gn.kinase, function(x) x[which(!duplicated(apply(x, 1, paste, collapse=' '))),] )
            
            ## sort
            s.acc.gn.kinase <- lapply(s.acc.gn.kinase, function(x)x[order(x$Score, decreasing=T), ])
            
            ## kinases
            list2env(s.acc.gn.kinase , envir=.GlobalEnv)
            WriteXLS(names(s.acc.gn.kinase), ExcelFileName=paste('NMF_KINASE_PHOSPHATASE_features_.xlsx', sep=''), FreezeRow=1, FreezeCol=1, SheetNames=make.unique(substr(names(s.acc.gn.kinase), 1, 20)), row.names=F, BoldHeaderRow=T, AutoFilter=T)
            
        }
        
        ## sort
        s.acc.gn.anno <- lapply(s.acc.gn.anno, function(x)x[order(x$Score, decreasing=T), ])
        
        
        ## ##########################################################
        ## export as txt
        ## all features
        list2env(s.acc.gn.anno , envir=.GlobalEnv)
        WriteXLS(names(s.acc.gn.anno), ExcelFileName=paste('NMF_features_N_', length(unique(unlist(s.acc2))),'.xlsx', sep=''), FreezeRow=1, FreezeCol=1, SheetNames=make.unique(substr(names(s.acc.gn.anno), 1, 20)), row.names=F, BoldHeaderRow=T, AutoFilter=T)
        
       
        ## #####################################################################################
        ##
        ##               contribution of each data type to NMF coefficients
        ##
        ## ######################################################################################
        data.types <- lapply(s.acc, function(x) table(sub('-.*','', x$Accession)))
        tmp <- lapply(data.types, names)
        tmp <- unique(unlist(tmp))
        data.types.mat <- matrix(0, nrow=length(data.types), ncol=length(tmp), dimnames=list(1:length(data.types), tmp) )
        for(i in 1:length(data.types))
            data.types.mat[i, names(data.types[[i]])] <- data.types[[i]]
        
        ## add MB subtype
        #idx.tmp <- match( cons.map.max, rownames(data.types.mat))
        idx.tmp <- match( names(cons.map), rownames(data.types.mat))
        rownames(data.types.mat)[idx.tmp] <- paste(rownames(data.types.mat)[idx.tmp], ' (', cons.map.max, ')', sep='' )
        
        
        col <- RColorBrewer::brewer.pal(length(tmp), "Set2")[1:length(tmp)]
        
        
        ## #############
        ## barplot
        names(col) <- tmp
        pdf('4.2_barplot_features_per_basis2.pdf')
        fancyBarplot(t(data.types.mat), srt=30, xlab='NMF basis', col=col, main='Proteogenomic features', ylab='No. features')
        legend('top', legend=names(col), fill=col, bty='n')
        legend('topright', legend=paste('total #:', length(unique(unlist( s.acc2 )))), bty='n')
        dev.off()
        
        ## ##############
        ## table
        m <- data.types.mat
        pheatmap(m, cluster_rows=F, cluster_cols=F, display_numbers=T, legend=F, number_format='%.0f', cellwidth=30, cellheight=30, fontsize=15, number_color='grey20', main='NMF features', scale='none', col='white', border_color='blue', filename='4.3_matrix_features_per_NMF.consensus.pdf')
        
        
        ## ###########################################################
        ## plot the actual expression values for the extracted features
        for(i in 1:length(s.gn)){
          
             m <- expr.org[ s.acc[[i]]$Accession, ]
            if(is.null(max.val)){
              m.max <- ceiling(max(abs(m), na.rm=T))
            } else {
              m.max <- max.val
              m[m > m.max] <- m.max
              m[m < -m.max] <- -m.max
            }
            pheatmap(m[, order(cdesc$NMF.consensus)], scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
                     annotation_colors=cdesc.color, filename=paste('5.',i-1,'_heatmap_expression_basis_', i, '.pdf', sep=''),  
                     breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cellwidth = cw,
                     main=names(NMF.consensus.col)[i], cluster_cols = F)
            
        }
        
        
        ## ############################################################
        ##
        ##                plot ALL markers
        ##
        feat.all <- unique(unlist(s.acc2))
        rdesc.feat <- matrix(0, nrow=length(feat.all), ncol=length(s.acc), dimnames=list( feat.all, names(s.acc) ))
        for(i in 1:length(s.acc))
          rdesc.feat[ s.acc[[i]]$Accession, i] <- 1
        rdesc.feat <- data.frame(rdesc.feat)
          
        m <- expr.org[ feat.all, ]
        
        if(is.null(max.val)){
          
          m.max <- ceiling(max(abs(m), na.rm=T))
        
          } else {
          
            m.max <- max.val
            m[m > m.max] <- m.max
            m[m < -m.max] <- -m.max
            
        }
        pheatmap(m, scale='none', fontsize_row=0, annotation_row = rdesc.feat, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, clustering_distance_col='euclidean', clustering_method='ward.D2', filename='6.0_heatmap_ALL_features.pdf', show_rownames=F, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cellwidth = cw)
        
        
        ## ############################################################
        ## plot ALL markersbut without CNV
        if('CNV' %in% names(data.str)){
            m.org <- m
            #m <- expr.org[ unique( unlist(s.acc)[-grep('CNV',  unlist(s.acc))] ), ]
            cnv.idx <- grep('CNV-', rownames(m.org))
            
            if(length(cnv.idx) > 0){
              m <- m[-cnv.idx, ]
              
              if(is.null(max.val)){
                m.max <- ceiling(max(abs(m), na.rm=T))
              } else {
                m.max <- max.val
                m[m > m.max] <- m.max
                m[m < -m.max] <- -m.max
              }
              pheatmap(m, scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, clustering_distance_col='euclidean', clustering_method='ward.D2', filename='6.1_heatmap_ALL_features_no_CNV.pdf', show_rownames=F,breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
            
        
              ## CNV only
              #m <- expr.org[ unlist(s.acc)[grep('CNV',  unlist(s.acc))] , ]
              m <-m.org[cnv.idx, ]
              m[m > 1] <- 1
              m[m < -1] <- -1
              m.max <- ceiling(max(abs(m), na.rm=T))
              pheatmap(m, scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, clustering_distance_col='euclidean', clustering_method='ward.D2', filename='6.2_heatmap_ALL_features_CNV_only.pdf', show_rownames=F,breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
            }
        }
        
        
        
        ## #############################################################################################
        ##
        ##                        positive conrols / genes of interest
        ##
        ## #############################################################################################
        if(!is.null(genes.of.interest))
        {
          
        
          pdf('barchart_positive_controls.pdf', 3, 3)
           for(i in genes.of.interest){
         
             ## extract features
             feat <- lapply(s.acc.gn, function(x){ 
               idx=grep(i, x$SYMBOL, value=F, ignore.case=F)
               paste(x$Type[idx], x$Accession[idx], sep='-' )
               })
       
             if(sum(sapply(feat, length)) > 0){
         
                 feat.expr <- expr.org[unlist(feat), ]
         
                 ## order: subtype + NMF
                 ord.idx <- with(cdesc, order( NMF.consensus))
         
         
                 ## if there is only a single row
                 if(is.null(dim(feat.expr))){
                     feat.expr <- feat.expr[ ord.idx ]
                     tmp <- feat.expr
         
                     dim(tmp) <- c(1, length(feat.expr))
                     colnames(tmp) <-  names(feat.expr)
                     rownames(tmp) <- unlist(feat)
                     m <- tmp
         
                     ## heatmap
                     m.max <- ceiling(max(abs(m), na.rm=T))
                     pheatmap( m,  annotation_col=cdesc[, rev(c(class.variable, variable.other) )], annotation_colors=cdesc.color, cluster_col=F, gaps_col=cumsum(table(cdesc[, 2])), main=paste('Positive ctrl:', i), filename=paste('7_heatmap_gene_of_interest_', sub('\\[.*|\\$','',i), '.pdf', sep=''), cellwidth=cw, cellheight=ch, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cluster_row=F)
         
                 } else {
                     m <- feat.expr[, ord.idx]
         
                     ## heatmap
                     m.max <-  ceiling(max(abs(m), na.rm=T))
                     pheatmap( m,  annotation_col=cdesc[ , rev(c(class.variable, variable.other) )],annotation_colors=cdesc.color, cluster_col=F, gaps_col=cumsum(table(cdesc[, 'NMF.consensus'])),main=paste('Positive ctrl:', i), filename=paste('7_heatmap_gene_of_interest_', sub('\\[.*|\\$','',i), '.pdf', sep=''), cellwidth=cw, cellheight=ch, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
                 }
         
                 feat.n <- unlist(lapply(feat, length))
         
                 ## barplot
                 pdf(paste('8_barplot_',sub('\\[.*|\\$','',i),'_features.pdf', sep=''), 3, 3)
                 fancyBarplot(feat.n, srt=45, xlab='NMF basis', names=1:length(s.acc), main=paste(i), ylab='No. features')
                 dev.off()
             }
          }
        # 
        }
        
        ## #######################################################################
        ##
        ##                      compare to PCA
        ##
        ## #######################################################################
        pch.vec <- c(15, 19, 17, 18, 8, 4, 3, 2, 1, 5:7, 9:14, 16, 20:25 )
        
        
        ## ############################################
        ## PCA on entire matrix
        pca <- prcomp( t(expr.full.org) )
     
        pc1 <- pca$x[,1]
        pc2 <- pca$x[,2]
        
        ## correct colors
        col.tmp <- cdesc.color[[class.variable]]
        col <- cdesc[,class.variable]
        #names(col) <- cdesc[,1]
        col <- sapply(col, function(i) col.tmp[i])
        
        ## plot
        pdf(paste('9.0_PCA.pdf'))
        plot(pc1, pc2, col=col,  cex=2, main=paste('PCA: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')), xlab='PC 1', ylab='PC2', pch=pch.vec[ kmeans(cbind(pc1, pc2), rank)$cluster] )
        pointLabel(pc1, pc2, labels=cdesc[,1], col=col, offset=50, method='SANN', cex=.5)
        legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
        
        
      #  fancyPlot(cdesc$CSresponseRank, pc1, xlab='CSresponseRank', ylab='Principle component 1', reg='linear', boxplots = 'n', cor = T, main=paste('PCA: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')))
      #  pointLabel(cdesc$CSresponseRank, pc1, labels=cdesc[,1], col=col, offset=50, method='SANN', cex=.5)
        
        dev.off()
        
    
        
        
        
        ## #############################################
        ## PCA on NMF features
        expr.nmf <- expr[feat.all, ]
        pca.nmf <- prcomp( t(expr.nmf) )
        
        
        pc1 <- pca.nmf$x[,1]
        pc2 <- pca.nmf$x[,2]
        
        col.tmp <- cdesc.color[[class.variable]]
        col <- cdesc[,class.variable]
        #names(col) <- cdesc[,1]
        col <- sapply(col, function(i) col.tmp[i])
        
        ## correct colors
        ## plot
        pdf(paste('9.1_PCA_NMF_features.pdf'))
        plot(pc1, pc2, col=col,  cex=2, main=paste('PCA on NMF features: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')), xlab='PC 1', ylab='PC2', pch=pch.vec[ kmeans(cbind(pc1, pc2), rank)$cluster])
        pointLabel(pc1, pc2, labels=cdesc[,1], col=col, offset=20, method='SANN', cex=.5)
        legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
        
       # fancyPlot(cdesc$CSresponseRank, pc1, xlab='CSresponseRank', ylab='Principle component 1', reg='linear', boxplots = 'n', cor = T, main=paste('PCA on NMF features: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')) )
      #  pointLabel(cdesc$CSresponseRank, pc1, labels=cdesc[,1], col=col, offset=50, method='SANN', cex=.5)
      
        dev.off()
        

        ## ################################################################
        ##
        ##   T-SNE plot of NMF markers
        ##
        ## ################################################################
        ## tsne
        tsne=Rtsne(expr.nmf, check_duplicates=F)
        
        ## tsne components
        tsneY = tsne$Y
      
        ## #################################
        ## color according to kmeans clustering
        tsne.clust.col=kmeans(tsneY, length(data.str))$cluster
        
        ## ##################################
        ## color according to data type
        tsne.dt.col=factor(sub('-.*','',rownames(expr.nmf)))
        
        ## ################################
        ## color according to NMF consensus
        tsne.nmf.col <- rownames(expr.nmf)
        names(tsne.nmf.col) <- tsne.nmf.col
        for(i in 1:length(s.acc))
            tsne.nmf.col[ s.acc[[i]]$Accession] <- NMF.consensus.col[ i ] %>% as.character #unlist(NMF.basis.col[as.character(i)])
        
        ## just the NMF class
        tsne.nmf.class <- rownames(expr.nmf)
        names(tsne.nmf.class) <- tsne.nmf.class
        for(i in 1:length(s.acc))
            tsne.nmf.class[s.acc[[i]]$Accession ] <- paste(i,  cons.map.max[i], sep=' ')
        
        
        
        d <- data.frame(
            tsne_1=tsneY[,1],
            tsne_2=tsneY[,2],
            Data_type=sub('-.*', '', rownames(expr.nmf)),
            Feature=rownames(expr.nmf),
            NMF_col=tsne.nmf.col,
            NMF_class=tsne.nmf.class
        )
        
        
        ##plot_ly(d, x=~tsne_1, y=~tsne_2, mode='markers', text=~Feature, color=~NMF_class)
        
        symbs <- c('circle', 'x', 'o', 'square', 'diamond', 'triangle-down', 'pentagon')
        
        #col.tmp <- unique(d$NMF_col)
        #names(col.tmp) <- unique(d$NMF_class)
        col.tmp <- NMF.consensus.col
        names(col.tmp) <- paste(1:length(col.tmp), names(col.tmp ))
            
        symb.tmp <- symbs[ 1:length(data.str) ]
        names(symb.tmp) <- unique(d$Data_type)
        
        ##
        dir.create('tsne')
        save(d, col.tmp, symb.tmp, file='tsne/tsne.RData')
        
        
        ## ####################################################
        ## Rmarkdown
        rmd <- c(paste("## T-SNE plot of NMF markers (n=",nrow(d),")\n
        ### Colored by NMF classification:\n
        ```{r echo=F}
        load('tsne.RData')
        library(plotly)
        plot_ly(d, x=~tsne_1, y=~tsne_2, mode='markers', type='scatter', text=~Feature, color=~NMF_class, colors=as.character(col.tmp))\n
        ```\n
        ### Colored by data type:\n
        ```{r echo=F}
        plot_ly(d, x=~tsne_1, y=~tsne_2, mode='markers', type='scatter', text=~Feature, color=~Data_type)\n
        ```
        \n
        ", sep=''))
        writeLines(rmd, con='tsne/tsne.rmd')
        try(rmarkdown::render('tsne/tsne.rmd'))
     
        setwd('..')   
    } 

}
    