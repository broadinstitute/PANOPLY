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
  make_option( c("-l", "--max_annot_levels"), action='store', type='numeric', dest='max_annot_levels', help='Maximum number of levels an annotation can have and be considered discrete.'), # default='geneSymbol'),
  make_option( c("-p", "--pval_signif"), action='store', type='numeric', dest='pval_signif', help='P-value threshold for significant enrichement.'), # default='geneSymbol'),
  make_option( c("-t", "--top_n_features"), action='store', type='numeric', dest='top_n_features', help='Number of driver-features to create expression-boxplots of, per cluster.'), # default='geneSymbol'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters.", default=NA),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # # for testing arguments
                       # args = c('--nmf_results',"/Users/wcorinne/Git/panoply-sandbox/hydrant/tasks/panoply_nmf/latest/outputs/panoply_nmf_workflow/4d729ec7-03f2-4abd-a01a-ef0c6c782e76/call-panoply_nmf/execution/nmf_res.Rdata",
                       #          '--rank_top',"6",
                       #          '--expr_comb',"/Users/wcorinne/Git/panoply-sandbox/hydrant/tasks/panoply_nmf/latest/outputs/panoply_nmf_workflow/4d729ec7-03f2-4abd-a01a-ef0c6c782e76/call-panoply_nmf/execution/glob-6bb1c1bc65aa92bfed4e9e23a70e6646/odg_test_combined_n61x26574.gct",
                       #          '--expr_comb_nn',"/Users/wcorinne/Git/panoply-sandbox/hydrant/tasks/panoply_nmf/latest/outputs/panoply_nmf_workflow/4d729ec7-03f2-4abd-a01a-ef0c6c782e76/call-panoply_nmf/execution/glob-5d28d7ac052751ffcd4ef752ba1f0fbf/odg_test_combinedNonNegative_n61x53147.gct",
                       #          '-g',"/Users/wcorinne/Downloads/groups-subset.csv",
                       #          '-a',"geneSymbol",
                       #          '-p',"0.01",
                       #          '-x',"odg_test")
)

#### Parse YAML Arguments ####
library(yaml)
# todo: read in YAML, instead of just using opt_cmd
opt = opt_cmd


#### Load Libraries ####
library(cmapR)
library(NMF)
library(limma)
# library(statmod)

library(tidyverse)
library(doParallel)
library(glue)
library(WriteXLS)

library(ggplot2)
library(ComplexHeatmap)
library(UpSetR)
library(Rtsne)


source("/prot/proteomics/Projects/PGDAC/src/modT.r")

# source modT.test.2class() and set_annot_colors() from R-utilities
# wd = getwd(); setwd("/prot/proteomics/Projects/R-utilities/"); # briefly change wd to r-utilities, so modT.r doesn't complain
# source("/prot/proteomics/Projects/R-utilities/modT.r"); setwd(wd); # source, and unset wd
source("/prot/proteomics/Projects/R-utilities/color-mod-utils.r")
# source("~/Git/proteomics-Rutil/color-mod-utils.r") # local testing

#### complex heatmap helper-function ####
# the only reason we're not just using Heatmap() directly is to:
# a.) easily toggle show / hide annotations
# b.) have defaults automatically set for certain common NMF plotting preferences (e.g. column-split by NMF.consensus)
MyComplexHeatmap <- function(matrix, heatmap_col, annot_df, annot_colors, na_color="#f7f7f7",
                             cluster_columns = FALSE, cluster_rows=TRUE,
                             show_column_names = FALSE, show_row_names = FALSE,
                             column_title_rot = 0, row_title_rot = 0,
                             show_heatmap_legend = FALSE,
                             row_split=NULL, 
                             column_split_annotCol="NMF.consensus",
                             show_annot = TRUE, ...){
  # create column-annotations separately
  if (show_annot) {
    ha <- HeatmapAnnotation(df=annot_df, col=annot_colors, na_col=na_color)
  } else {ha=NULL}
  
  # Heatmap arguments
  args_list = list(matrix = matrix,
                   top_annotation = ha,
                   row_split = row_split, # class vector with groupings for rows
                   column_split = annot_df[[column_split_annotCol]], # class vector with grouping from columns; taken from provided annot_df
                   row_dend_side = 'right',
                   cluster_columns = cluster_columns, cluster_rows=cluster_rows,
                   show_column_names = show_column_names, show_row_names = show_row_names,
                   column_title_rot = column_title_rot, row_title_rot = column_title_rot,
                   show_heatmap_legend = show_heatmap_legend, ...)
  if (!missing(heatmap_col) && !is.null(heatmap_col)) { # if we have custom heatmap-colors
    args_list <- modifyList(args_list, list(col = heatmap_col)) # add to arguments list
  }
  
  do.call(Heatmap, args_list) %>% #plot heatmap
    return() # return heatmap
}



###########################################################
##         Data Import / Parameter Setup
###########################################################

prefix = paste0(opt$output_prefix,"_K", opt$rank_top, "_") # prefix for filenames
cw=10; ch=cw # heatmap cell size
id_col = "Sample.ID" # id column used in groups file

#### NMF-Input Dataset Import ####
gct_expr_comb = parse_gctx(opt$expr_comb)
gct_expr_comb_nn = parse_gctx(opt$expr_comb_nn)

#### Annotations / Groups File Pre-processing ####
if (!is.null(opt$groups_file)) { # if we have a groups file
  #### import groups file ####
  groups = read_csv(opt$groups_file) %>% # import groups file
    column_to_rownames(id_col) # set id column to rownames
  #### ensure that groups file contains all of our samples ####
  if ( length(setdiff(gct_expr_comb_nn@cid, rownames(groups))) > 0  ) { # if there are Sample.IDs that are not in our groups file
    stop(paste0("Invalid groups file; the following ",id_col," are missing from groups file:\n", # throw an error
                setdiff(gct_expr_comb_nn@cid, rownames(groups)))) # and print the missing samples
  }
} else { # otherwise
  groups = gct_expr_comb@cdesc # otherwise use full cdesc
  rownames(groups) = gct_expr_comb@cid # set rownames to 'id' column / cid
}

#### Annotation Color Values ####
if (FALSE) {
  # todo: assign colors from YAML
  colors = NULL
  idx.is_discrete = NULL # determine which annotations are discrete'
  # continuous may be separate
} else {
  # assign colors
  colors.full = set_annot_colors(groups, continuous.return_function = TRUE)
  # pull out the color vectors, named with values
  # consider: this should probably happen in set_annot_color()
  colors = sapply(colors.full, function(annot) { # wrangle list of lists into a list of:
    annot_colors = annot$colors # vectors with colors
    names(annot_colors) = annot$vals # where each color is named for its corresponding annotation value
    return(annot_colors)
  })
  # create vector with discrete/non-discrete TRUE/FALSE values, for each annotation
  idx.is_discrete = sapply(colors.full, function(annot) {annot$is_discrete}) 
}
# pull out discrete vs continuous annots
annots.is_discrete.full = names(idx.is_discrete[which(idx.is_discrete)])
annots.is_continuous = names(idx.is_discrete[which(!idx.is_discrete)])
# exclude "discrete" variables with too many categories
annots.excluded = sapply(annots.is_discrete.full, function(annot) {
  length(unique(groups[[annot]])) > opt$max_annot_levels # determine which discrete annotaitons have too many values
}) %>% which() %>% names() # convert named TRUE/FALSE vector to vector of names
if (length(annots.excluded) > 0) {
  cat(paste(glue("The following {length(annots.excluded)} annotations have more than {opt$max_annot_levels} discrete levels, and will be excluded from analysis & figures:"), "\n",
            paste(annots.excluded, collapse = ", "), "\n"))
}
annots.is_discrete = setdiff(annots.is_discrete.full, annots.excluded) # remove that annotation from annot.is_discrete


#### Load NMF Results ####
load(opt$nmf_results) # loads in res.rank object

#### Retrieve Basis & Coefficient Matrix ####

res = res.rank[[opt$rank_top]]
nmf_fit = res@fit

NMF.basis = predict(res)
NMF.consensus = predict(res, 'consensus') %>% # predict results
  names() %>% NMF.basis[.] # map back onto basis matrix
#consider: only valid for non-bayesian!! also, could this be made with the H-vector data?
consensusmap(res, filename = paste0(prefix, "consensusMap.pdf")) # idk what makes this meaningful compared to a different way of determining cluster membership but

basis.mat = nmf_fit@W
coef.mat = nmf_fit@H


###########################################################
##         H-Matrix Analysis / cluster membership
###########################################################

#### determine cluster membership ####
coef_norm  <- apply(coef.mat, 2, function(x) x/sum(x)) # normalize each column to sum of column
coef_consensus <- apply(coef_norm, 2, which.max) # which cluster had the maximum membership
coef_membership <- apply(coef_norm, 2, max) # maximum membership fraction
# mindiff core membership
coef_core_member <- apply(coef_norm, 2, function(x){
  if ( sum(x == max(x)) > 1 ) return(FALSE) # if more than one value is the maximum, return FALSE
  min_diff = min(max(x) - setdiff(x, max(x))) # minimum difference between our largest coefficient, and our other coefficients
  min_diff > 1/as.numeric(opt$rank_top) # is this difference less than 1/rank?
})

NMF.annots = data.frame(row.names = colnames(coef.mat),
                        NMF.consensus = NMF.consensus, # use predict results, rather than max-membership
                        NMF.cluster.membership = round(coef_membership, 3),
                        NMF.core.member = coef_core_member) %>%
  arrange(NMF.consensus, desc(NMF.cluster.membership)) # order by cluster, then cluster-membership
#### save NMF membership info to CSV ####
write.csv(NMF.annots, # sorted by cluster-membership
          paste0(prefix, "clusterMembership.csv"))

#### add NMF.annots to groups-file ####
groups.full = merge(groups, NMF.annots, by='row.names') %>% # append NMF.annots to groups file / cdesc
  arrange(NMF.consensus) %>%  # order by cluster-membership
  mutate_at(annots.is_continuous, as.numeric) %>% # mutate all continuous annotations to numerics
  column_to_rownames('Row.names') # re-add rownames

#### assign colors for NMF annotations ####
colors.full.NMF = set_annot_colors(NMF.annots, continuous.return_function = TRUE)
# colors.full.NMF$NMF.consensus$colors = c('#004488', '#DDAA33', '#BB5566') # manually overwrite to really clearly distinct colors...
colors.NMF = sapply(colors.full.NMF, function(annot) {
  annot_colors = annot$colors
  names(annot_colors) = annot$vals
  return(annot_colors)
})

#### filter groups.full to core members ####
groups.core = filter(groups.full, NMF.core.member) # get just core

#### write H-matrix & normalized H-matrix to a GCT file, with full groups-file cdesc ####
## unmodified H-matrix
gct.H <- new('GCT')
gct.H@mat <- coef.mat[, match(rownames(NMF.annots), colnames(coef.mat))] # un-normalized H-matrix (sorted by cluster-membership, i.e. as they appear in NMF.annots object))
gct.H@cdesc <- groups.full[match(colnames(gct.H@mat), rownames(groups.full)), ] # append full cdesc (sorted to match matrix)
gct.H@rid <- paste0('C',1:nrow(gct.H@mat))
gct.H@cid <- colnames(gct.H@mat)
write_gct(gct.H, ofile=paste0(prefix, 'H')) # write H to GCT
## normalized H-matrix
gct.H.norm = gct.H # initialize
gct.H.norm@mat <- apply(gct.H@mat, 2, function(x) x/max(x)) # normalized to the max of each column
write_gct(gct.H.norm, ofile=paste0(prefix, 'H_normToMax')) # write H-norm to GCT

#### plot H.norm heatmap, with and without clustering ####
for (cluster_cols in c(TRUE,FALSE)) {
  pdf(paste0(prefix, 'heatmap_coeficientMatrixNorm', # save to PDF
             ifelse(cluster_cols, '_clustered', ''),'.pdf'), 20,8)
  MyComplexHeatmap(gct.H.norm@mat, # plot normalized heatmap
                   heatmap_col = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), # colorscale for heatmap
                   show_heatmap_legend = T, show_column_names = T, # show heatmap & sample IDs
                   select(gct.H.norm@cdesc, -c(annots.excluded)), # with all annotations
                   c(colors, colors.NMF), show_annot=T, # using appropriate annotation colors
                   cluster_columns = cluster_cols, #optionally cluster columns
                   cluster_rows = F) %>%
    draw(., annotation_legend_side='right')
  dev.off()
  
  # obselete pheatmap code
  # pheatmap::pheatmap(gct.H.norm@mat, color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), # plot normalized coefficient matrix
  #                    annotation_col = select(gct.H.norm@cdesc, -c(annots.excluded, "NMF.core.member")), # add all annotations
  #                    # consider: binary NMF.core.member column seems to break pheatmap...
  #                    annotation_colors = c(colors, colors.NMF), # add all annotation colors
  #                    show_rownames = TRUE,
  #                    cluster_cols = cluster_cols, cluster_rows = FALSE,
  #                    cellwidth=cw, cellheight=ch,
  #                    gaps_col = cumsum(table(gct.H.norm@cdesc$NMF.consensus)),
  #                    filename = paste0(prefix, 'heatmap_coeficientMatrixNorm', # save to PDF
  #                                      ifelse(cluster_cols, '_clustered', ''),'.pdf')) # with or without clustering suffix
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

#### Calculate Enrichement on each annotation ####
clust.enrich.df = data.frame(cluster=numeric(), # initialize empty df
                             class=character(),
                             class.level=character(),
                             pval=double())
for (annot in annots.is_discrete) {
  # compute enrichement of each annot.value (excluding NAs)
  tmp = CalcClustEnrich(groups.core$NMF.consensus, groups.core[[annot]]) %>%
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

#consider: only valid for non-bayesian?
rank.sil <- silhouette(res)
pdf(paste0(prefix, 'silhouettePlot_clusterMembership.pdf')) # open PDF
plot(rank.sil, main=paste('K=', opt$rank_top, sep=''),
     col=colors.full.NMF$NMF.consensus$colors)
dev.off() # close PDF



###########################################################
##         Boxplots of Continuous Features
###########################################################

if ( length(annots.is_continuous) > 0 ) {
  pdf(paste0(prefix, 'boxplots_continuousAnnots.pdf')) # open PDF
  for (annot in annots.is_continuous) { # for each continuous annotation
    p = ggplot(groups.core, aes(x=NMF.consensus, y=!!sym(annot), # compare the values across clusters
                                color = NMF.consensus)) + scale_color_manual(values = colors.full.NMF$NMF.consensus$colors) +
      geom_boxplot() + # make a boxplot
      ggtitle(annot) +
      ggpubr::stat_compare_means(comparisons=combn(1:opt$rank_top,2,simplify = FALSE)) # add wilcox.test comparison between each cluster
    plot(p) # plot as a PDF page
  }
  dev.off() # close PDF
} else {
  print("No continuous annotations found in groups-file. Skipping boxplots")
}





###########################################################
##         W-Matrix Analysis / driver features 
###########################################################
#### Calculate Feature Scores to identify driver features ####
method="max" # try kim first, and then use max if it fails
feature.scores <- featureScore(basis.mat, method=method) %>% # calculate feature-score for features in basis-matrix
  sort(decreasing=TRUE) # sort features from highest to lowest score
# todo: add in the histogram for feature scores

#### create W-matrix GCT ####
## unmodified W-matrix
gct.W <- new('GCT')
gct.W@mat <- basis.mat[names(feature.scores),] # basis.mat, sorted by feature scores
gct.W@rdesc <- gct_expr_comb_nn@rdesc[match(rownames(gct.W@mat), # append full rdesc, sorted to match W matrix
                                            rownames(gct_expr_comb_nn@rdesc)),] 
gct.W@rdesc$feature.scores <- feature.scores # will be in the correct order since gct.W@mat was sorted to match feature.score
gct.W@rid <- rownames(gct.W@mat)
gct.W@cid <- paste0('C',1:ncol(gct.W@mat))
write_gct(gct.W, ofile=paste0(prefix, 'W')) # write W to GCT

#### create normalized W-matrix GCT ####
## normalized W-matrix
basis_norm <- apply(basis.mat, 1,  function(x)x/sum(x)) %>% # normalize rowwise
  t() # transform to get back into original orientation
gct.W.norm = gct.W # initialize
gct.W.norm@mat <- basis_norm[match(gct.W.norm@rid, # overwrite data-matrix with normalized data, sorted to match rid
                                   rownames(basis_norm)), ]
write_gct(gct.W.norm, ofile=paste0(prefix, 'W_rowNorm')) # write W-norm to GCT

#### create directional/combined W-matrix GCT ####
W.norm.dir = sweep(gct.W.norm@mat, 1, ifelse(gct.W.norm@rdesc$sign=="neg", -1, 1), "*") # sign each normalized row
## combine pos/neg features; pick abs.max for each feature, sample-wise
cl <- makeCluster(detectCores()) # set up cluster, so this can be done in parallel
registerDoParallel(cl) # start cluster
features = unique(gct.W.norm@rdesc$id_unsigned) # get all unique unsiged-feature IDs
W.tmp.list <- foreach(feat = features, # for each unsigned feature id
                      .final = function(x) setNames(x, features) # set names of each output to the feature
) %dopar% { # do in prallel
  W.tmp = W.norm.dir[dplyr::filter(gct.W.norm@rdesc, id_unsigned==feat)$id, , drop = F] # select pos / neg features 
  apply(W.tmp, 2, function(x) x[ which.max(abs(x)) ]) # get abs.max for each feature
}
stopCluster(cl) # stop cluster

W.norm.comb = do.call(rbind, W.tmp.list) # combine list into data.frame


## save to GCT
gct.W.norm.comb = gct.W.norm # initialize w/ W.norm
gct.W.norm.comb@mat <- W.norm.comb # overwrite matrix with normalized & directional values, combined across pos/neg
# reformat rdesc
rdesc.tmp <- dplyr::select(gct.W@rdesc, -c("id","sign")) %>% # drop signed ID & sign
  distinct() %>% # subset to distinct rows
  dplyr::rename(id = id_unsigned) # set id_unsigned to id
gct.W.norm.comb@rdesc <- rdesc.tmp[match(rownames(gct.W.norm.comb@mat), # overwrite rdesc with new rdesc, in order of new matrix
                                         rdesc.tmp$id),] 
rownames(gct.W.norm.comb@rdesc) <- gct.W.norm.comb@rdesc$id # set rownames to unsigned ID
# overwrite rid
gct.W.norm.comb@rid <- rownames(gct.W.norm.comb@mat) # overwrite rid
write_gct(gct.W.norm.comb, ofile=paste0(prefix, 'W_rowNorm_combined_signed')) # write to GCT



###########################################################
#### Identify Driver Features ####
driver.features <- extractFeatures(basis.mat, method=method)
driver.features.byName = lapply(driver.features, function(index_vec) {rownames(basis.mat)[index_vec]}) # signed feature ID
driver.features.byName.unsigned = lapply(driver.features.byName, function(feat_vec) {
  gct_expr_comb_nn@rdesc[feat_vec,"id_unsigned"] %>% unique()  }) # unsigned feature IDs (unique)

#### T-test on Driver Features, on Expression Data in/out of cluster ####
driver.features.list = list() # initialize list for driver-feature t-test-dataframes for each cluster
for (clust in 1:as.numeric(opt$rank_top)) { # for each clusters' driver features
  #### set up t-test inputs ####
  # ft.vec = driver.features.byName[[clust]] # non-negative feature IDs
  # expr_mat = gct_expr_comb_nn@mat # non-negative feature expression
  ft.vec = driver.features.byName.unsigned[[clust]] # get driver featuress for this cluster
  if ( sum(is.na(ft.vec))!=0 ) {driver.features.list[[clust]]=NA; next} # if we have no driver features in this cluster, return NA val and skip test
  
  feat.cl <- gct.W.norm.comb@rdesc[ft.vec,] %>% # initialize a dataframe to pass back, if t-test fails
    mutate(cluster = clust) # add cluster 
  expr_mat = gct_expr_comb@mat # original feature expression
  ## subset expression matrix to driver-features, and add id column
  ft.expr = expr_mat[ft.vec, , drop=FALSE] %>% # get corresponding SIGNED feature expression
    as.data.frame() %>% rownames_to_column("feature.id") # convert to dataframe and turn rownames into a column
  ## create class vector assigning sample-columns to ingroup vs outgroup
  class.vec = colnames(expr_mat) %in% rownames(filter(NMF.annots, NMF.consensus==clust)) %>% # get TRUE/FALSE for whether a sample is in this cluster
    ifelse(., paste0("C",clust), "allother") # convert TRUE to cluster-label and FALSE to "allother" (using expr_mat since ft.expr has a feature-id column)
  names(class.vec) <- colnames(expr_mat) # add sample_id names back in (using expr_mat since ft.expr has a feature-id column)
  # note: it's important that the class.vec and expr_mat have columns IN THE SAME ORDER; names() are not actually taken into account...

  #### perform moderated test ####
  res.tmp <-  try({
    label = glue("C{clust}.over.rest")
    modT.test.2class( ft.expr, groups=class.vec,
                      id.col='feature.id', label=label )$output %>%
      .[, c(1, grep(label, colnames(.)))] # subset to JUST the id-column & t-test-columns
  })
  if(class(res.tmp) == 'try-error') { # if t-test failed
    driver.features.list[[clust]] <- feat.cl # just return the feature annotation df
  } else { # otherwise
    driver.features.list[[clust]] <- full_join(feat.cl, res.tmp, by="id") %>% # merge feature-annotation df with t-test df
      dplyr::arrange(!!sym(paste0("P.Value.", label))) # arrange by p-value column
  }
}

# todo: add the excel sheet for the driver.features.list


#### simplify t-test results into dataframe ####
driver.features.df = lapply(driver.features.list, function(df) {
  if (is.na(df)) return(NULL)
  df %>% # subset each data-frame to
    dplyr::select(c(!ends_with(".over.rest"), # annotation columns
                    starts_with("P.Value."))) %>% # and p-value column
    dplyr::rename("pval" = starts_with("P.value.")) # rename pvalue column for Easy Stacking
}) %>% do.call(rbind, .) %>% # rbind list into dataframe
  dplyr::mutate(!!paste0("signif_at_",opt$pval_signif) := pval < opt$pval_signif) # boolean for significance

#### Benjamini-Hochberg correction ####
driver.features.bh = mutate(driver.features.df,
                            bh.pval = p.adjust(driver.features.df$pval, method = 'BH'), # apply BH correction
                            !!paste0("bh.signif_at_",opt$pval_signif) := bh.pval < opt$pval_signif) # boolean if significant




#### Write Driver Feature + P-value to CSV ####
## all features
write.csv(driver.features.bh,
          paste0(prefix, "driverFeatures.csv"))
## significant only
driver.features.sigFeatOnly = filter(driver.features.bh, #filter just to significant features
                                     !!sym(paste0("bh.signif_at_",opt$pval_signif))) 
write.csv(driver.features.sigFeatOnly, # filter to 
          paste0(prefix, "driverFeatures_sigFeatOnly.csv"))

#### Write Excel File with full test results + BH correction ####
# add bh correction to full driver-feature df 
driver.features.list.bh <- driver.features.list # initialize with original list
names(driver.features.list.bh) = paste0('C', 1:length(driver.features.list.bh))
for (clust in 1:length(driver.features.list.bh)) { # for each cluster
  if (is.na(driver.features.list.bh[[clust]])) { # if we find a cluster with no features
    driver.features.list.bh[[clust]] <- NULL # set that entry to NULL
    next # and skip to next cluster
  } 
  # otherwise, append two new columns
  driver.features.list.bh[[clust]]$BH.adj.P.Value <- dplyr::filter(driver.features.bh, cluster==clust)$bh.pval # append a column with the BH corrected pvalues
  driver.features.list.bh[[clust]]$BH.adj.P.Value_signif <- dplyr::filter(driver.features.bh, cluster==clust)$bh.signif_at # append a column with whether the value was significant
}
WriteXLS(driver.features.list.bh,
         ExcelFileName=paste0(prefix, "driverFeatures_byCluster.xlsx"),
         FreezeRow=1, FreezeCol=1,
         SheetNames=make.unique(substr(names(driver.features.list.bh), 1, 20)),
         row.names=F, BoldHeaderRow=T, AutoFilter=T)

#### Write Excel File, subset to Kinase Only ####
driver.features.list.kinase = lapply(driver.features.list.bh, function(df) {
  filter(df, grepl('2\\.7\\.1[0-2]|3\\.1\\.3\\.(16|48)', ENZYME)) # subset each DF of driver-features to kinases only
})
WriteXLS(driver.features.list.kinase,
         ExcelFileName=paste0(prefix, "driverFeatures-kinase-phosphotase_byCluster.xlsx"),
         FreezeRow=1, FreezeCol=1,
         SheetNames=make.unique(substr(names(driver.features.list.kinase), 1, 20)),
         row.names=F, BoldHeaderRow=T, AutoFilter=T)


###########################################################
#### create Significant-Driver-Feature Heatmap ####
for (clust in unique(driver.features.sigFeatOnly$cluster)) { # for each cluster that had significant features
  # subset expression matrix to drivers for this cluster
  expr.drivers = gct_expr_comb@mat[filter(driver.features.sigFeatOnly, cluster==clust)$id, # filter to significant features for this cluster
                                   match(rownames(groups.full), colnames(gct_expr_comb@mat)),  # sort in order of groups.full
                                   drop = FALSE] # do not simplify matrix to array if there's one row
  # append signficant features to heatmap
  if (clust == unique(driver.features.sigFeatOnly$cluster)[1]) { # if we're on the first cluster-heatmap 
    # initialize hm.concat
    hm.concat <- MyComplexHeatmap(expr.drivers, NULL, # plot expression of drivers, using default color-palette
                                  select(groups.full, -c(annots.excluded)), # with all annotations
                                  c(colors, colors.NMF), # with custom annot-colors
                                  show_row_names = T, cluster_rows = F, # show feature_ids, but don't cluster
                                  row_names_side = "left") # put row names on left side
  } else { # otherwise
    hm.concat <- hm.concat %v% # append heatmaps vertically
      MyComplexHeatmap(expr.drivers, NULL, # plot expression of drivers, using default color-palette
                       groups.full, colors[annots.is_discrete], # plot next heatmap
                       show_annot=FALSE, # don't re-plot annotations
                       show_row_names = T, cluster_rows = F, # don't cluster rows
                       row_names_side = "left") # put row names on left side
  }
}

#### draw Heatmap to Files ####
pdf(paste0(prefix, "driverFeatures_complexHeatmap.pdf"), 20, 16)
draw(hm.concat, annotation_legend_side='bottom')
dev.off()

png(paste0(prefix, "driverFeatures_complexHeatmap.png"), width=20, height=16, units = 'in', res=100)
draw(hm.concat, annotation_legend_side='bottom')
dev.off()


###########################################################
#### create UpSet Plot of driver-feature-overlap b/t clusters ####
if ( length(unique(driver.features.bh$cluster)) > 1 ) { # if more than one cluster has driver features
  upset.mat <- driver.features.sigFeatOnly %>% # use only significant features
    dplyr::select(c(id, cluster)) %>% # select relevant columns
    dplyr::mutate(observed = 1) %>% # make placeholder "observed" column, to use as the value in the pivot_wider() function
    pivot_wider(id_col = 'id', # use feature-id as id
                names_from='cluster', names_prefix='C', # create column for each cluster with driver features
                values_from='observed', values_fill=0) %>% # fill with observed (1) if driver-feature was seen in cluster & fill wtih missing with 0 instead of NA if driver was NOT seen in cluster
    column_to_rownames('id') # convert feature-id to rowname
  
  p <- upset(upset.mat, # create upset-plot
             point.size = 4, text.scale = 1.5, nintersects = NA,
             nsets = dim(upset.mat)[2])
  
  pdf(paste0(prefix, "driverFeatures_UpSet.pdf")); p; dev.off(); # plot as PDF
  png(paste0(prefix, "driverFeatures_UpSet.png")); p; dev.off(); # plot as PNG
}

###########################################################
#### calculate Contributions of Each Ome ####
ft_countByOme <- driver.features.sigFeatOnly %>% # use only significant features
  dplyr::select(c(ome_type, cluster)) %>% # select relevant columns
  group_by(cluster, ome_type) %>% summarise(n_feat = n()) %>% # get number of features per-ome per-cluster
  ungroup() %>% complete(ome_type, cluster, fill = list(n_feat = 0)) %>% # add rows for ome_type/cluster combinations with n_feat = 0
  mutate(cluster = paste0('C',cluster))
# colnames(ft_countByOme) <- paste0('C', colnames(ft_countByOme)) # append C to cluster-number
#### write CSV of contributions per ome ####
write.csv(ft_countByOme,
          paste0(prefix, "driverFeatures_contributionsByOme.csv"))
#### plot Barplot of Significant Features per Ome ####
p <- ggplot(ft_countByOme, aes(x = as.factor(cluster), y = n_feat, fill=ome_type)) +
  geom_bar(position="dodge", stat="identity") + # grouped barplot with n_feat per clut per ome
  geom_text(aes(label = n_feat), # add labels with number of features
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) + # position label above each column
  labs(x = "NMF Basis", y = "Number of Features") + # x axis and y axis labels
  ggtitle('Proteogenomic Features by Cluster') + # add title
  theme_minimal() + theme(panel.grid.major.x = element_blank()) # adjust formatting

ggsave(paste0(prefix, "driverFeatures_contributionsByOme_barplot.pdf"),
       plot = p, device="pdf") # plot as PDF
ggsave(paste0(prefix, "driverFeatures_contributionsByOme_barplot.png"),
       plot = p, device="png") # plot as PNG








###########################################################
#### create Barplots of Top Driver Features by cluster ####

driver.features.topNFeat <- driver.features.sigFeatOnly %>% # take significant driver-features
  group_by(cluster) %>% # group by cluster
  slice_min(order_by = bh.pval, # order by p-value
            n = opt$top_n_features) %>% # select n features with lowest p-values
  ungroup() # ungroup



pdf(paste0(prefix, "driverFeatures_expressionBoxplots.pdf"))
driverFeats = driver.features.topNFeat$id # get vector of driver features we wanna plot 
for (ft in driverFeats) {
  # collect the info we need for our boxplot
  expr.df = data.frame(ft.expression = gct_expr_comb@mat[ft,], # get expression data
                       # sample.id = gct_expr_comb@cid, # get corresponding sample name (unnecessary but convenient)
                       cluster = NMF.annots[gct_expr_comb@cid, # in the order of the samples
                                            'NMF.consensus']) #get cluster membership
  annot.df = filter(driver.features.sigFeatOnly, id==ft) # filter annotations to relevant driver-feature
  
  if (!is.null(opt$gene_col)) { # if we have the gene column, use it for the boxplot title
    boxplot.title = glue('Expression of {annot.df[[opt$gene_col]]} in {annot.df$ome_type} ({annot.df$id_og})')
  } else { # otherwise just use the feature ID
    boxplot.title = glue('Expression of feature {annot.df$id_og} ({annot.df$ome_type})')
  }
  
  boxplot.colors = rep('black', opt$rank_top) #initialize color-scheme with all black
  boxplot.colors[[annot.df$cluster]] = "red" # color the expression red in the cluster the driver feature drives
  
  p <- ggplot(expr.df, aes(x = as.factor(cluster), y = ft.expression)) + # boxplot with ft-expression by ome
    geom_boxplot(color = boxplot.colors) + # boxplots of expression data
    geom_hline(yintercept = 0, color='grey', linetype='dashed') + # add dashed line at 0
    labs(x = "NMF Basis", y = "Feature Expression", # x axis and y axis labels
         title = boxplot.title, # add title
         subtitle = glue("Driver Feature in C{annot.df$cluster}\nFeature Score: {signif(annot.df$feature.scores, 3)}\nAdj-P-Val: {signif(annot.df$bh.pval, 3)}")) +
    theme_minimal() + theme(panel.grid.major.x = element_blank()) # adjust formatting
  plot(p) # plot as PDF page
}
dev.off()


###########################################################
##         TSNE Comparison of Clusters
###########################################################

#### try tSNE on All Features ####
expr_mat_core <- gct_expr_comb@mat[, rownames(filter(NMF.annots,  # expression matrix
                                                     NMF.core.member))] %>% # of JUST core features
  t() # samples as rows
tsne <- try(Rtsne(expr_mat_core,
                  perplexity = floor((nrow(expr_mat_core) - 1) / 3), # to avoid perplexity error
                  check_duplicates=F))
tsne.subtitle = "Clustered on all features" # add appropriate subtitle

#### if tSNE clustering failed on ALL features, try again on ONLY driver features ####
if(class(tsne) == 'try-error'){ # if tSNE failed on ALL features
  warning("TSNE clustering failed on all features! Trying again on driver-features only.") # try again on driver-features
  #### subset expression matrix to JUST driver features ####
  expr_mat_core <- gct_expr_comb@mat[driver.features.sigFeatOnly$id, # driver features only
                                     rownames(filter(NMF.annots,  # expression matrix
                                                     NMF.core.member))] %>% # of JUST core features
    t() # samples as rows
  #### try tSNE again ####
  tsne <- try(Rtsne(expr_mat_core,
                    perplexity = floor((nrow(expr_mat_core) - 1) / 3), # to avoid perplexity error
                    check_duplicates=F))
  
  tsne.subtitle = "Clustered on driver-features only." # add appropriate subtitle
} 

#### plot results ####
if (!class(tsne) == 'try-error') { # assuming tSNE ran properly
  #### create dataframe for plotting ####
  tsne.df <- data.frame(tsne_x = tsne$Y[,1],
                        tsne_y = tsne$Y[,2],
                        sample_id = rownames(expr_mat_core),
                        cluster = NMF.annots[rownames(expr_mat_core),
                                             "NMF.consensus"])
  #### plot and save ####
  p = ggplot(tsne.df, aes(tsne_x, tsne_y, color = cluster)) +
    geom_point() +
    labs(title = "tSNE Clustering on Core NMF Samples", # title
         subtitle = tsne.subtitle) + # subtitle, describing the dataset
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + # plot 0,0 axes in black
    theme_minimal()
  # ggsave("/opt/input/tmp.png", plot=p, device="png") # for manual docker testing
  ggsave(paste0(prefix, "tSNEclustering.pdf"),
         plot=p, device = "pdf") # saveas pdf
  ggsave(paste0(prefix, "tSNEclustering.png"),
         plot=p, device = "png") # as png
} else {
  warning("TSNE clustering failed.") # try again on driver-features
}


###########################################################
##         Tar all Files
###########################################################

output.files = list.files(pattern=prefix, full.names = T)
tar('NMF_results.tar.gz',
    files = output.files, compression = "gzip", tar="tar")


