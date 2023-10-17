#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))

# specify command line arguments
option_list <- list(
  make_option( c("-f", "--filenames"), action='store', type='character',  dest='ome_gcts_string', help='string of gct filenames, semicolon-separated'),
  make_option( c("-l", "--labels"), action='store', type='character',  dest='ome_labels_string', help='string of gct filenames, semicolon-separated'),
  make_option( c("-v", "--var"), action='store', type='numeric',  dest='var_pca_expl', help='Explained variance by PCA. Used to extract the number of PCs explaining the specified fraction of variance in the multiomics data matrix.', default = 0.9),
  make_option( c("-t", "--tol"), action='store', type='numeric',  dest='perc_tol', help='Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.', default = 0.05),
  make_option( c("-z", "--z_score_mode"), action='store', type='character', dest='z_score_mode', help='z-scoring mode: row, col, rowcol', default='rowcol')
  
)

###########################################################
## parse command line parameters'
opt_cmd <- parse_args( OptionParser(option_list=option_list) )

z_score_mode <- opt_cmd$z_score_mode

##################################################################
## explained variance by PCA
## - used to extract the n PCs explaining 'var_pca_expl' of
##   the variance in the multiomics data matrix
## - For these PCs the feature contributions
##   are calculated
## - cumulative contribution across all PCs
##   is used to balance features from different omes
var_pca_expl <- opt_cmd$var_pca_expl 

## stopping criterion: maxim
perc_tol <- opt_cmd$perc_tol 

library(pacman)
p_load(cmapR)
p_load(magrittr)

########################################
## import
## wdl-compatible
import_gct_from_wdl <- function( ome_gcts_string, ome_labels_string ){
  
  # parse strings into vectors
  ome_gct_str = unlist(strsplit(ome_gcts_string, ","))
  ome_labels = unlist(strsplit(ome_labels_string, ","))
  # ome_labels = sapply(ome_gct_str, USE.NAMES = FALSE, function(filepath) { basename(filepath) %>% gsub(".gct","", .) }) # assume that file-name is ome-label
  
  if (length(ome_gct_str)!=length(ome_labels)) stop("Different number of GCT files and ome-labels detected. Please check your inputs.")
  
  ## import
  ome_gcts <- lapply(ome_gct_str, parse_gctx)
  names(ome_gcts) <- ome_labels
  
  out <- list(ome_gcts=ome_gcts, ome_labels=ome_labels)
  
  return(out)
}

#####################################################
##
##                   fancyBarplot
##
## a better customizable version of 'barplot'
##
## 20110217 added support for matrices
## 20140201 'cex.numb'
#####################################################
fancyBarplot <- function(x, space=0.2, 
                         ylim=c(0, max(x, na.rm=T)+ 0.15*max(x, na.rm=T)), 
                         ndigits=3, 
                         add.numb=T,
                         numb=c('counts', 'dev.median.perc'),
                         cex.numb=.9, srt=90, border=NA, ...)
{
  numb <- match.arg(numb)
  
  #########################################
  #         vector
  #########################################
  if(length(dim(x)) <= 1){
    # make the plot
    barplot(x, space=space, ylim=ylim, border=border, ...)
    
    # add numbers
    if(add.numb){
      if(numb=='counts'){
        numb.to.add <- x
        ## add to plot
        text(  seq( 0.5+space, (0.5+space)+((length(numb.to.add)-1)*(1+space)), 1+space ),
               x+(0.05*max(x)), round(numb.to.add, ndigits), srt=srt, pos=4, offset=-0.05, cex=cex.numb )
      }
      
      if(numb == 'dev.median.perc'){
        numb.median <- median(x, na.rm=T)
        numb.to.add <- sapply(x, function(y) 100*(y-numb.median)/numb.median )
        abline(h=numb.median, col='darkblue', lty='dashed')
        #legend('topleft', legend=c(paste('Median:', round(numb.median, ndigits))), bty='n', text.col = 'darkblue')
        mtext(c(paste('Median:', round(numb.median, ndigits))),side = 2, at = numb.median, col = 'darkblue', cex = 0.8)
        ## add to plot
        text(  seq( 0.5+space, (0.5+space)+((length(numb.to.add)-1)*(1+space)), 1+space ),
               x+(0.05*max(x)), paste(round(numb.to.add, ndigits),'%', sep=''), srt=srt, pos=4, offset=-0.05, cex=cex.numb )
      }
    }
    
  } else{
    
    #########################################
    #        matrix
    #########################################
    offset = -0.1
    barplot(x, ylim=ylim, beside=T, border=border, ...)
    
    if(add.numb){
      p.tmp=1
      for(p in 1:dim(x)[2]){
        if( sum(is.na(x[, p])) < nrow(x) ){
          text( (p.tmp : (p.tmp + dim(x)[1]-1) ) + offset, max(x[, p], na.rm=T)+0.05*( max(x[, p], na.rm=T)), round(x[,p], ndigits), adj=c(1,1), pos=4, cex=cex.numb, srt=srt  )
        }
        p.tmp <- p.tmp+dim(x)[1]+1
      }
    }
    
  }
}

#################################################################################
## function to calculate feature importance across meaningful PCs, i.e. explaining
## 'var_pca_expl' variance
pca_feature_importance <- function(mat_z,              ## z-scored data matrix - features in rows, samples in columns            
                                   var_pca_expl=0.9    ## var explained in the entire dataset; used to extract meaningful PCs
) {
  
  ## PCA
  res.pca <- prcomp( t(mat_z), scale. = F, center = F)
  
  ## PCs to explain 'var_pca_expl' variance 
  res.pca.summ <- summary(res.pca)
  pcs_filt <- which(res.pca.summ$importance[3, ] <= var_pca_expl)
  
  ## Helper function 
  var_coord_func <- function(loadings, comp.sdev){
    loadings*comp.sdev
  }  
  
  ## Compute Coordinates
  loadings <- res.pca$rotation
  sdev <- res.pca$sdev
  var.coord <- t(apply(loadings, 1, var_coord_func, sdev))
  
  ## Compute Cos2
  var.cos2 <- var.coord^2
  
  # Compute contributions
  comp.cos2 <- apply(var.cos2, 2, sum)
  contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
  var.contrib <- t(apply(var.cos2, 1, contrib, comp.cos2))
  
  
  ## contributions for PCs explaining "var_pca_expl" of variance
  var.contrib.filt <- var.contrib[, pcs_filt]
  
  ## relative contributions
  var.contrib.filt.sum <- apply(var.contrib.filt, 1, sum)
  var.contrib.filt.perc <- var.contrib.filt.sum / sum(var.contrib.filt.sum)
  
  return(var.contrib.filt.perc)
}


##################################################
## function to balance contributions of each ome
balance_omes <- function(   contrib,           ## numeric, vector of contributions for each features 
                            ome,               ## character, vector of same length as 'contrib' mapping each feature to an ome
                            perc_tol=0.05,     ## stop criterion
                            plot=TRUE,         ## if TRUE barplots depicting the contributions and number of features before and after balancing
                            verbose=TRUE,       ## if TRUE some progress will be reported on the command line
                            z_score_mode=c('row', 'col', 'rowcol')
){    
  
  z_score_mode <- match.arg(z_score_mode)
  
  ## target: equal contributions (within 'perc_tol')
  ## of each ome
  ome_perc_target <- rep(1/length(unique(ome)), length(unique(ome)))
  
  ## split by ome
  ome_list <- tapply(contrib, ome, function(x) x) 
  
  ome_perc_org <- sapply(ome_list, sum) 
  ome_perc_org <- ome_perc_org / sum(ome_perc_org)
  
  cond <- TRUE
  iter <- 1
  
  while(cond) {
    if(verbose) cat('iter:', iter, '\n')
    
    ## total sum
    ome_sum <- sapply(ome_list, sum)
    
    ## percent variance with current features 
    ome_perc <- ome_sum / sum(ome_sum)
    
    ## deviation from target
    perc_deviation <- ome_perc - ome_perc_target
    
    ome_adjust_idx <- which(perc_deviation > perc_tol)
    
    if(length(ome_adjust_idx) == 0) cond <- FALSE
    
    if(length(ome_adjust_idx) > 0) {
      
      ome_adjust <- sort(perc_deviation[ ome_adjust_idx ], decreasing = T)
      
      ## only update the ome with maximal contribution
      for(ome_tmp in names(ome_adjust)[1] ){
        
        if(verbose) cat('__', ome_tmp, '___', ome_adjust[ome_tmp], '\n')        
        
        ome_tmp_var <- ome_list[[ ome_tmp ]] 
        ome_tmp_var <- ome_tmp_var[ -which.min(ome_tmp_var) ]
        ome_list[[ome_tmp]] <- ome_tmp_var
        
      }
    }
    iter <- iter + 1
  }
  
  if(plot){
    
    #source('https://raw.githubusercontent.com/karstenkrug/R-code/main/my_plots.r')
    
    col <- c('grey20', 'grey80')
    
    #pdf(paste0('balance-contrib-var-', var_pca_expl, '-tol-', perc_tol, '.pdf'), 10,5)
    pdf(paste0('balance-omes.pdf'), 10,5)
    par(mfrow=c(1,2))
    ## contribution
    fancyBarplot(rbind(ome_perc_org,
                       ome_perc),
                 srt=30,
                 main='Fractional contributions',
                 ylab='Fraction of variance', 
                 col=col
    )
    abline(h=unique(ome_perc_target), lty='dashed', lwd=2)
    
    legend('top', legend=c('all', 'balanced'), fill=col, bty='n', ncol = 2)
    ## number of features
    fancyBarplot(rbind(table(ome),
                       sapply(ome_list, length)),
                 srt=30,
                 main='Number of features',
                 ylab='Number',
                 col=col
    )
    legend('top', legend=c('all', 'balanced'), fill=col, bty='n', ncol = 2)
    
    dev.off()
  }
  return(ome_list)
}


#####################################
## import GCT files
##ome_gcts <- lapply(ome_labels, parse_gctx)
gct_imp <- import_gct_from_wdl( opt_cmd$ome_gcts_string, opt_cmd$ome_labels_string )
ome_gcts <- gct_imp$ome_gcts
ome_labels <- gct_imp$ome_labels

## extract data matrices
gct_mat <- lapply(ome_gcts, function(x) x@mat)

## remove features with missing values
gct_mat_fullquant <- lapply(gct_mat, function(m){
  na_idx <- apply(m, 1, function(x) sum(is.na(x)))
  keep_idx <- which(na_idx == 0)
  m[keep_idx, ]
})

## merge tables: 
## - mat
## - rdesc
for(ome in names(gct_mat_fullquant)){
  
  mat_tmp <- gct_mat_fullquant[[ome]]
  cid_tmp <- colnames(mat_tmp)
  rdesc_tmp <- data.frame(id=paste0(ome, '_', rownames(mat_tmp)), ome=rep(ome, nrow(mat_tmp)))
  
  if(ome == names(gct_mat_fullquant)[1]) {
    
    mat <- mat_tmp
    rdesc <- rdesc_tmp
    
  } else {
    
    cid <- colnames(mat)
    cid_common <- intersect(cid, cid_tmp)
    
    mat <- rbind( mat[, cid_common], mat_tmp[, cid_common] )
    rdesc <- rbind(rdesc, rdesc_tmp)
  }
}

################################
## z-score multi-omic data matrix

## - z-score by row
## - followed by z-score columns
if(z_score_mode == 'rowcol'){
  mat_z <- apply(mat, 1, function(x) (x-mean(x))/sd(x))
  mat_z <- apply(mat_z, 2, function(x) (x-mean(x))/sd(x))
  mat_z <- t(mat_z)
}
## - z-score by row
if(z_score_mode == 'row'){
  mat_z <- apply(mat, 1, function(x) (x-mean(x))/sd(x))
  mat_z <- t(mat_z)
}
## - z-score by column
if(z_score_mode == 'col'){
  mat_z <- apply(mat, 2, function(x) (x-mean(x))/sd(x))
  mat_z <- t(mat_z)
}


#################################
## calculate feature importance
var.contrib.filt.perc <- pca_feature_importance(mat_z, var_pca_expl = var_pca_expl)

#################################
## balance contributions from each ome
features_balanced <- balance_omes(var.contrib.filt.perc, rdesc$ome, perc_tol = perc_tol)

#################################
## filter tables
for(o in names(features_balanced)){
  index = which(ome_labels==o) # get the proper index for this file, so that WDL's glob() function pulls files in the right order
  gct_org <- ome_gcts[[o]]
  feat_filt <- features_balanced[[o]] %>% names
  
  ## filter gct
  gct_filt <- subset_gct(gct_org, rid=feat_filt)
  
  ## export
  fn <- paste0(sprintf("%02d", index), "-", # prepend padded-index to filename, to ensure that glob() pulls files in the correct order
               o, # use label in file-name
               '-balanced-contrib.gct') # append balanced-contrib.gct

  fn <- file.path(tempdir(), fn)
  write_gct(gct_filt, ofile = fn, appenddim = F)
  file.copy(fn, '.', overwrite = T)
  
}
