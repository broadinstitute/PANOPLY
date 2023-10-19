#!/usr/bin/env Rscript
#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### Input Parameters ####
  make_option( c("-n", "--nmf_results"), action='store', type='character',  dest='nmf_results', help='Rdata object with NMF Clustering results (res.rank object).'),
  make_option( c("-e", "--expr_comb"), action='store', type='character',  dest='expr_comb', help='GCT file with combined expression data.'),
  make_option( c("-f", "--expr_comb_nn"), action='store', type='character',  dest='expr_comb_nn', help='GCT file with combined (non-negative) expression data, used for NMF.'),
  make_option( c("-r", "--rank_top"), action='store', type='numeric',  dest='rank_top', help='Best number of clustering / rank.'),
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest. If not provided, all annotations in the cdesc will be analyzed.'),
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  #### Post-Processing Parameters ####
  make_option( c("-p", "--pval_signif"), action='store', type='numeric', dest='pval_signif', help='P-value threshold for significant enrichement.'), # default='geneSymbol'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters.", default=NA),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # # for testing arguments
                       # args = c('--nmf_results',"/opt/input/nmf_res.Rdata",
                       #          '--rank_top',"3",
                       #          '--expr_comb',"/opt/input/LSCC_GCP_combined_n155x101.gct",
                       #          '--expr_comb_nn',"/opt/input/LSCC_GCP_combinedNonNegative_n155x202.gct",
                       #          '-a',"geneSymbol",
                       #          '-p',"0.01",
                       #          '-x',"LSCC_GCP")
)

#### Parse YAML Arguments ####
library(yaml)
# todo: read in YAML, instead of just using opt_cmd
opt = opt_cmd


#### Load Libraries ####
library(cmapR)
library(NMF)
library(tidyverse)
library(ggplot2)

# source annotation-color-function from R-utilities
source("/prot/proteomics/Projects/R-utilities/color-mod-utils.r")





###########################################################
##         Data Import / Parameter Setup
###########################################################

prefix = paste0(opt$output_prefix,"_K", opt$rank_top, "_") # prefix for filenames
cw=10; ch=cw # heatmap cell size

#### Data Import / Annotations of Interest ####
gct_expr_comb = parse_gctx(opt$expr_comb)
gct_expr_comb_nn = parse_gctx(opt$expr_comb_nn)
if (!is.null(opt$groups_file)) {
  # todo: use groups file
} else {
  groups = gct_expr_comb@cdesc # otherwise use full cdesc
  rownames(groups) = gct_expr_comb@cid # set rownames to 'id' column / cid
}

#### Annotation Color Values ####
if (FALSE) {
  # todo: assign colors from YAML
  colors = NULL
  colors.is_discrete = NULL # determine which annotations are discrete
} else {
  # assign colors
  colors.full = set_annot_colors(groups)
  # pull out the color vectors, named with values
  # consider: this should probably happen in set_annot_color()
  colors = sapply(colors.full, function(annot) {
    annot_colors = annot$colors
    names(annot_colors) = annot$vals
    return(annot_colors)
  })
  colors.is_discrete = sapply(colors.full, function(annot) {annot$is_discrete})
}
# pull out discrete vs continuous annots
annots.is_discrete = names(colors.is_discrete[which(colors.is_discrete)])
annots.is_continuous = names(colors.is_discrete[which(!colors.is_discrete)])


#### Load NMF Results ####
load(opt$nmf_results) # loads in res.rank object

#### Retrieve Basis & Coefficient Matrix ####

res = res.rank[[opt$rank_top]]
nmf_fit = res@fit

NMF.basis = predict(res)
NMF.consensus = predict(res, 'consensus')
#consider: only valid for non-bayesian!! also, could this be made with the H-vector data?
consensusmap(res@consensus, filename = paste0(prefix, "consensusMap.pdf")) # idk what makes this meaningful compared to a different way of determining cluster membership but

basis.mat = nmf_fit@W
coef.mat = nmf_fit@H


###########################################################
##         H-Matrix Analysis / cluster membership
###########################################################

#### determine cluster membership ####
coef_norm  <- apply(coef.mat, 2, function(x) x/sum(x))
coef_consensus <- apply(coef_norm, 2, which.max) # which cluster had the maximum membership
coef_membership <- apply(coef_norm, 2, max) # maximum membership fraction
# mindiff core membership
coef_core_member <- apply(coef_norm, 2, function(x){
  min_diff = min(max(x) - setdiff(x, max(x))) # minimum difference between our largest coefficient, and our other coefficients
  min_diff > 1/as.numeric(opt$rank_top) # is this difference less than 1/rank?
})

NMF.annots = data.frame(row.names = colnames(coef.mat),
                        NMF.consensus = as.factor(coef_consensus),
                        NMF.cluster.membership = round(coef_membership, 3),
                        NMF.core.member = coef_core_member) %>%
  arrange(NMF.consensus) # order by cluster-membership
#### save NMF membership info to CSV ####
write.csv(NMF.annots, # sorted by cluster-membership
          paste0(prefix, "clusterMembership.csv"))

#### add NMF.annots to cdesc ####
cdesc.full = merge(gct_expr_comb@cdesc, NMF.annots, by='row.names') %>% # append NMF.annots
  arrange(NMF.consensus) %>%  # order by cluster-membership
  mutate_at(annots.is_continuous, as.numeric) %>% # mutate all continuous annotations to numerics
  column_to_rownames('Row.names') # re-add rownames

#### assign colors for NMF annotations ####
colors.full.NMF = set_annot_colors(NMF.annots)
colors.full.NMF$NMF.consensus$colors = c('#004488', '#DDAA33', '#BB5566') # manually overwrite to really clearly distinct colors...
colors.NMF = sapply(colors.full.NMF, function(annot) {
  annot_colors = annot$colors
  names(annot_colors) = annot$vals
  return(annot_colors)
})

#### filter cdesc.full to core members ####
cdesc.core = filter(cdesc.full, NMF.core.member) # get just core

#### write H-matrix & normalized H-matrix to a GCT file, with full cdesc ####
## unmodified H-matrix
gct.H <- new('GCT')
gct.H@mat <- coef.mat[,match(rownames(cdesc.full), colnames(coef.mat))] # un-normalized H-matrix, sorted to match cdesc
gct.H@cdesc <- cdesc.full # append full cdesc (matrix will match this order)
gct.H@rid <- paste0('C',1:nrow(gct.H@mat))
gct.H@cid <- rownames(gct.H@cdesc)
write_gct(gct.H, ofile=paste0(prefix, 'H')) # write H to GCT
## normalized H-matrix
gct.H.norm = gct.H # initialize
gct.H.norm@mat <- apply(coef.mat[,match(rownames(cdesc.full), colnames(coef.mat))], # un-normalized H-matrix, sorted to match cdesc
                        2, function(x) x/max(x)) # normalized to the max of each column
write_gct(gct.H.norm, ofile=paste0(prefix, 'H_normToMax')) # write H-norm to GCT

#### plot H.norm heatmap, with and without clustering ####
for (cluster_cols in c(TRUE,FALSE)) {
  pheatmap::pheatmap(gct.H.norm@mat, color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), # plot normalized coefficient matrix
                     annotation_col = gct.H.norm@cdesc[c(annots.is_discrete, "NMF.consensus")], # add discrete annotations
                     annotation_colors = c(colors[annots.is_discrete], NMF.consensus=colors.NMF$NMF.consensus), # use annotation color-scheme
                     show_rownames = TRUE,
                     cluster_cols = cluster_cols, cluster_rows = FALSE,
                     cellwidth=cw, cellheight=ch,
                     gaps_col = cumsum(table(gct.H.norm@cdesc$NMF.consensus)),
                     filename = paste0(prefix, 'heatmap_coeficientMatrixNorm', # save to PDF
                                       ifelse(cluster_cols, '_clustered', ''),'.pdf')) # with or without clustering suffix
}



###########################################################
##         enrichment analysis
###########################################################

## function for calculating enrichement on a cluster vector & class vector; returns table with cluster, annot_value, and p_value from comparing in-group vs out-group
CalcClustEnrich <- function(clust.vec,  ## vector of cluster labels
                            class.vec   ## vector of class labels, same order and length as 'clust'
                            
){
  clust.vec <- as.character(clust.vec)
  class.vec <- as.character(class.vec)
  
  res.per.clust <- data.frame(cluster=numeric(),
                              class.level=character(),
                              pval=double())
  
  ## loop over clusters
  for(clust in sort(unique(clust.vec))){
    
    ## in cluster
    class.vec.clust <- class.vec[ which(clust.vec == clust)]
    
    ## not in cluster
    class.vec.not.clust <- class.vec[ which(clust.vec != clust)]
    
    ## loop over class levels
    for(class.level in na.omit(unique(class.vec))){
      
      x11 <- sum(class.vec.clust %in% class.level)
      x12 <- sum(!(class.vec.clust %in% class.level))
      x21 <- sum( class.vec.not.clust %in%  class.level )
      x22 <- sum(!(class.vec.not.clust %in% class.level ))
      
      # Fisher's p
      pval = fisher.test( rbind(c(x11, x12), c(x21, x22)), alternative = 'greater')$p.value
      
      # append row
      res.per.clust = rbind(res.per.clust,
                            data.frame(cluster=clust,
                                       class.level=class.level,
                                       pval=pval))
    }
  }
  return(res.per.clust)
}

## stop if nothing to test #consider: is this necessary?
if(length(annots.is_discrete) == 0) stop('No categorial variables left to perform enrichment analysis.')

#### Calculate Enrichement on each annotation ####
clust.enrich.df = data.frame(cluster=numeric(), # initialize empty df
                             class=character(),
                             class.level=character(),
                             pval=double())
for (annot in annots.is_discrete) {
  # compute enrichement of each annot.value (excluding NAs)
  tmp = CalcClustEnrich(cdesc.core$NMF.consensus, cdesc.core[[annot]]) %>%
    mutate(class = annot) %>% relocate(class, .before=class.level) %>% # append class as column
    mutate(!!paste0("signif_at_",opt$pval_signif) := pval < opt$pval_signif) # boolean if significant
  # append to clust.df
  clust.enrich.df = rbind(clust.enrich.df, tmp) # append each annotation's enrichement
}

#### FWER correction (Bonferroni) ####
n_tests = dim(clust.enrich.df)[1] # number of fisher tests performed
clust.enrich.fwer = mutate(clust.enrich.df,
                           fwer.pval = sapply(clust.enrich.df$pval*n_tests, min, 1), # apply BF correction, cap at 1
                           !!paste0("fwer.signif_at_",opt$pval_signif) := fwer.pval < opt$pval_signif) # boolean if significant

#### Write Cluster Enrichement to CSV ####
write.csv(clust.enrich.fwer,
          paste0(prefix, "clusterEnrichement.csv"))



###########################################################
##         Silhouette Plots
###########################################################

#consider: only valid for non-bayesian!? also, could this be made with the H-vector data?
rank.sil <- silhouette(res.rank[[opt$rank_top]])
pdf(paste0(prefix, 'silhouettePlot_clusterMembership.pdf')) # open PDF
plot(rank.sil[[opt$rank_top]], main=paste('K=', opt$rank_top, sep=''),
     col=colors.full.NMF$NMF.consensus$colors)
dev.off() # close PDF



###########################################################
##         Boxplots of Continuous Features
###########################################################

pdf(paste0(prefix, 'boxplots_continuousAnnots.pdf')) # open PDF
for (annot in annots.is_continuous) { # for each continuous annotation
  p = ggplot(cdesc.full, aes(x=NMF.consensus, y=!!sym(annot), # compare the values across clusters
                             color = NMF.consensus)) + scale_color_manual(values = colors.full.NMF$NMF.consensus$colors) +
    geom_boxplot() + # make a boxplot
    ggtitle(annot) +
    ggpubr::stat_compare_means(comparisons=combn(1:opt$rank_top,2,simplify = FALSE)) # add wilcox.test comparison between each cluster
  plot(p) # plot as a PDF page
}
dev.off() # close PDF




###########################################################
##         W-Matrix Analysis / driver features 
###########################################################

## row-normalize
basis_norm <- apply(basis.mat, 1,  function(x)x/sum(x)) %>% # normalize rowwise
  t() # transform to get back into original orientation

# consider: combine up/down, get highest reading for each feature

#### write W-matrix & normalized W-matrix to a GCT file ####
## unmodified H-matrix
gct.W <- new('GCT')
gct.W@mat <- basis.mat # un-normalized H-matrix, sorted to match rdesc
# gct.W@rdesc <- NULL #rdesc[match(rownames(rdesc), rownames(basis.mat))] # append full cdesc (matrix will match this order) # consider: fix the rdesc bullshit
gct.W@rid <- rownames(basis.mat)
gct.W@cid <- paste0('C',1:ncol(gct.W@mat))
write_gct(gct.W, ofile=paste0(prefix, 'W')) # write H to GCT
## normalized W-matrix
gct.W.norm = gct.W # initialize
gct.W.norm@mat <- basis_norm
write_gct(gct.W.norm, ofile=paste0(prefix, 'W_rowNormalized')) # write H-norm to GCT


#### Identify Significant Features ####

#### Calculate Feature Scores to identify driver features ####
# consider: idk what this is
feature.scores <- featureScore(basis.mat, method='max')
driver.features <- extractFeatures(basis.mat, method='max')
driver.features.byName = lapply(driver.features, function(index_vec) {rownames(basis.mat)[index_vec]})

#### T-test on Driver Features, on Expression Data in/out of cluster ####
driver.features.df = data.frame(cluster=numeric(), # initialize empty df
                                feature=character(),
                                pval=double())
for (clust in 1:as.numeric(opt$rank_top)) { # for each clusters' driver features
  ft.vec = driver.features.byName[[clust]]
  for (ft in ft.vec) { # for each feature in that cluster
    # consider: account for positive / negative expression. should we use expr_comb or expr_comb_nn?
    ft.expr = gct_expr_comb_nn@mat[ft,] # get corresponding feature expression
    class.vec = colnames(gct_expr_comb_nn@mat) %in% rownames(filter(NMF.annots, NMF.consensus==clust)) # get TRUE/FALSE for whether a sample is in this cluster
    
    # consider: modified T test?? moderated T test????
    ttest_out = t.test(ft.expr[class.vec], ft.expr[!class.vec]) # t-test on expression values in the cluster, and out of the cluster
    
    driver.features.df = rbind(driver.features.df,
                               data.frame(cluster = clust,
                                          feature = ft,
                                          pval=ttest_out$p.value))
    
  }
}
driver.features.df = driver.features.df %>%  # add boolean for significance
  mutate(!!paste0("signif_at_",opt$pval_signif) := pval < opt$pval_signif)

#### FWER correction (Bonferroni) ####
n_tests = dim(driver.features.df)[1] # number of fisher tests performed
driver.features.fwer = mutate(driver.features.df,
                              fwer.pval = sapply(driver.features.df$pval*n_tests, min, 1), # apply BF correction, cap at 1
                              !!paste0("fwer.signif_at_",opt$pval_signif) := fwer.pval < opt$pval_signif) # boolean if significant
#### Write Driver Feature + P-value to CSV ####
## all features
write.csv(driver.features.fwer,
          paste0(prefix, "driverFeatures.csv"))
## significant only
driver.features.sigFeatOnly = filter(driver.features.fwer, #filter just to significant features
                                     !!sym(paste0("fwer.signif_at_",opt$pval_signif))) 
write.csv(driver.features.sigFeatOnly, # filter to 
          paste0(prefix, "driverFeatures_sigFeatOnly.csv"))
sign_driver_features = driver.features.sigFeatOnly$feature


#### Heatmap of Significant Features ####
sign_driver_features_unsigned = unique(gct_expr_comb_nn@rdesc[sign_driver_features,'id_unsigned'])
expr.drivers = gct_expr_comb@mat[sign_driver_features_unsigned, # filter to significant features
                                 match(rownames(cdesc.full), colnames(gct_expr_comb@mat))] # sort in order of cdesc.full
for (cluster_cols in c(TRUE,FALSE)) {
  pheatmap::pheatmap(expr.drivers, # plot expression of driver features
                     color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100)), # red-blue color palette
                     breaks = seq(-max(abs(expr.drivers)), max(abs(expr.drivers)), length.out=101),
                     annotation_col = cdesc.full[c(annots.is_discrete, "NMF.consensus")], # add discrete annotations
                     annotation_colors = c(colors[annots.is_discrete], NMF.consensus=colors.NMF$NMF.consensus), # use annotation color-scheme
                     show_rownames = TRUE,
                     cluster_cols = cluster_cols, cluster_rows = FALSE,
                     cellwidth=cw, cellheight=ch,
                     gaps_col = cumsum(table(gct.H.norm@cdesc$NMF.consensus)),
                     filename = paste0(prefix, 'heatmap_signFeatures', # save to PDF
                                       ifelse(cluster_cols, '_clustered', ''),'.pdf')) # with or without clustering suffix
}
#### Barplot of Significant Features per Ome ####

