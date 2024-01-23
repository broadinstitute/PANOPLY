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
  make_option( c("-n", "--nmf_results"), action='store', type='character',  dest='nmf_results', help='Tar file containing combined expression GCT, combined non-negative expression GCT, and Rdata object with NMF Clustering results (res.rank object).'),
  make_option( c("-r", "--rank_top"), action='store', type='numeric',  dest='rank_top', help='Best number of clustering / rank.'),
  # make_option( c("-e", "--expr_comb"), action='store', type='character',  dest='expr_comb', help='GCT file with combined expression data.'),
  # make_option( c("-f", "--expr_comb_nn"), action='store', type='character',  dest='expr_comb_nn', help='GCT file with combined (non-negative) expression data, used for NMF.'),
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest. If not provided, all annotations in the cdesc will be analyzed.'),
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  #### Post-Processing Parameters ####
  make_option( c("-p", "--pval_signif"), action='store', type='numeric', dest='pval_signif', help='P-value threshold for significant enrichement.'), 
  make_option( c("-q", "--feature_fdr"), action='store', type='numeric', dest='feature_fdr', help='Max FDR threshold for feature-selection 2-sample T-test.'),
  make_option( c("-l", "--max_annot_levels"), action='store', type='numeric', dest='max_annot_levels', help='Maximum number of levels an annotation can have and be considered discrete.'), # default='geneSymbol'),
  make_option( c("-t", "--top_n_features"), action='store', type='numeric', dest='top_n_features', help='Number of driver-features to create expression-boxplots of, per cluster.'), # default='geneSymbol'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters."),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # for testing arguments
                       # args = c('--nmf_results',"/opt/input/nmf_res.Rdata",
                       #          '--rank_top',"3",
                       #          '--expr_comb',"/opt/input/odg_test_combined_n61x26574.gct",
                       #          '--expr_comb_nn',"/opt/input/odg_test_combinedNonNegative_n61x53147.gct",
                       #          '-y',"/opt/input/master-parameters.yaml",
                       #          '-x',"odg_test")
)


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


#### Parse YAML Arguments ####
opt = opt_cmd # initialize options with command line options
if ( !is.null(opt$yaml_file) ) {
  #### read in yaml ####
  library(yaml)
  yaml_out <- read_yaml(opt$yaml_file)
  #### locate the NMF-postprocessing parameters section ####
  if ( !is.null(yaml_out$panoply_nmf_postprocess) ) { # if we have an NMF parameters section
    yaml_nmf =  yaml_out$panoply_nmf_postprocess # read in those parameters
  } else { # if the section is missing
    if ( !is.null(yaml_out$panoply_mo_nmf) ) { # check for the deprecated mo_nmf section
      cat(glue("\nWARNING: The parameter file '{opt$yaml_file}' is missing the 'panoply_nmf_postprocess' section, but includes the deprecated 'panoply_mo_nmf' parameter section. Parameters will be read in from 'panoply_mo_nmf', but please consider updating your parameters file!"))
      yaml_nmf =  yaml_out$panoply_mo_nmf # read in those parameters
      if (is.null(opt$top_n_features)) opt$top_n_features = 25 # manually provide default for top_n_features, since this parameter previously did not exist
    } else { # otherwise, stop
      stop(glue("The parameter file '{opt$yaml_file}' does not contain an NMF parameters section. Please check that the yaml file contains the appropraite 'panoply_nmf' section."))
    }
  }
  #### overwrite the command-line parameters ####
  # global parameters
  if (is.null(opt$gene_col)) opt$gene_col = yaml_out$global_parameters$gene_mapping$gene_id_col
  # postprocessing parameters
  if (is.null(opt$pval_signif)) opt$pval_signif = yaml_nmf$ora_pval
  if (is.null(opt$feature_fdr)) opt$feature_fdr = yaml_nmf$feature_fdr
  if (is.null(opt$max_annot_levels)) opt$max_annot_levels = yaml_nmf$ora_max_categories
  if (is.null(opt$top_n_features)) opt$top_n_features = yaml_nmf$top_n_features
} else { # if no YAML was provieded
  # check if any necessary parameters are missing
  if( any(sapply(list(opt$pval_signif, opt$feature_fdr,
                      opt$max_annot_levels, opt$top_n_features), 
                 is.null)) ) { # if we have at least one missing parameter
    stop("Master Parameter yaml-file is missing. Please either provide a master-parameters file, or manually provide all NMF parameters.") # error and stop
  }
}

# print parameters
cat("\n\n####################\nNMF POST-PROCESSING PARAMETERS\n\n")
print(opt)
save(file = "nmf_postprocess_opt.Rdata", opt)
cat("####################\n\n")



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
                             row_dend_side = 'right',
                             show_annotation_legend = TRUE,
                             show_heatmap_legend = FALSE,
                             row_split=NULL, 
                             column_split_annotCol="NMF.consensus",
                             show_annot = TRUE, ...){
  # create column-annotations separately
  if (show_annot) {
    ha <- HeatmapAnnotation(df=annot_df, col=annot_colors, na_col=na_color,
                            show_legend=show_annotation_legend)
  } else {ha=NULL}
  
  # Heatmap arguments
  args_list = list(matrix = matrix,
                   top_annotation = ha,
                   row_split = row_split, # class vector with groupings for rows
                   column_split = annot_df[[column_split_annotCol]], # class vector with grouping from columns; taken from provided annot_df
                   row_dend_side = row_dend_side,
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

cat("\n\n####################\nData Import\n\n")

# untar results file
untar(opt$nmf_results) # untar NMF results
# locate necessary inputs from tar, and error if ANY are missing
file_nmf_res = list.files(pattern = 'nmf_res.Rdata') # locate Rdata object containing res.rank
if (length(file_nmf_res) != 1) stop("Tarfile {opt$nmf_results} is missing nmf_res.Rdata object, containing res.rank output of nmf().")
file_expr_comb = list.files(pattern = "_combined_n.+.gct") # locate combined GCT file
if (length(file_expr_comb) != 1) stop("Tarfile {opt$nmf_results} is missing GCT file with combined expression data.")
file_expr_comb_nn = list.files(pattern = "_combinedNonNegative_n.+.gct") # locate combined nonnegative GCT file
if (length(file_expr_comb_nn) != 1) stop("Tarfile {opt$nmf_results} is missing GCT file with combined, non-negative expression data.")

  

prefix = paste0(opt$output_prefix,"_K", opt$rank_top, "_") # prefix for filenames
# cw=10; ch=cw # pheatmap cell size
sample_id_col = "Sample.ID" # id column used in groups file

#### NMF-Input Dataset Import ####
gct_expr_comb = parse_gctx(file_expr_comb)
gct_expr_comb_nn = parse_gctx(file_expr_comb_nn)

#### Annotations / Groups File Pre-processing ####
if (!is.null(opt$groups_file)) { # if we have a groups file
  #### import groups file ####
  groups = read_csv(opt$groups_file) %>% # import groups file
    column_to_rownames(sample_id_col) # set id column to rownames
  #### ensure that groups file contains all of our samples ####
  if ( length(setdiff(gct_expr_comb_nn@cid, rownames(groups))) > 0  ) { # if there are Sample.IDs that are not in our groups file
    stop(paste0("Invalid groups file; the following ",sample_id_col," are missing from groups file:\n\n", # throw an error
                paste(setdiff(gct_expr_comb_nn@cid, rownames(groups)),
                      collapse = ", "))) # and print the missing samples
  }
} else { # otherwise
  groups = gct_expr_comb@cdesc # otherwise use full cdesc
  rownames(groups) = gct_expr_comb@cid # set rownames to 'id' column / cid
}

#### Annotation Color Values ####

# automatically assign colors
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
# pull out discrete vs continuous annots
annots.is_discrete.full = names(idx.is_discrete[which(idx.is_discrete)])
annots.is_continuous = names(idx.is_discrete[which(!idx.is_discrete)])
# exclude "discrete" variables with too many categories
annots.excluded = sapply(annots.is_discrete.full, function(annot) {
  length(unique(groups[[annot]])) > opt$max_annot_levels # determine which discrete annotaitons have too many values
}) %>% which() %>% names() # convert named TRUE/FALSE vector to vector of names
if (length(annots.excluded) > 0) {
  cat(paste(glue("The following {length(annots.excluded)} annotations have more than {opt$max_annot_levels} discrete levels, and will be excluded from analysis & figures:"), "\n\n",
            paste(annots.excluded, collapse = ", "), "\n\n"))
}
annots.is_discrete = setdiff(annots.is_discrete.full, annots.excluded) # remove that annotation from annot.is_discrete


# overwrite automatic colors with YAML-file colors, if applicable
if (!is.null(opt$yaml_file)) {
  yaml_colors = yaml_out$groups.colors # get colors from YAML
  for (annot in annots.is_discrete) { # for each discrete-annotation we intend to plot
    if (!is.null( yaml_colors[[annot]] )) { # if the annotation has colors assigned in the YAML file
      colors[[annot]] <- unlist(yaml_colors[[annot]]) # overwrite the automatic values with the YAML values
    }
  }
}


#### Load NMF Results ####
load(file_nmf_res) # loads in res.rank object from Rdata file, stored in opt$nmf_results

#### Retrieve Basis & Coefficient Matrix ####

res = res.rank[[opt$rank_top]]
nmf_fit = res@fit

NMF.basis = predict(res)
NMF.consensus = predict(res, 'consensus') %>% # predict results
  names() %>% NMF.basis[.] # map back onto basis matrix
cat("\n\n####################\nCluster Membership-- Consensus Mapping\n\n")
consensusmap(res, filename = paste0(prefix, "consensusMap.pdf")) # idk what makes this meaningful compared to a different way of determining cluster membership but
consensusmap(res, filename = paste0(prefix, "consensusMap.png")) # idk what makes this meaningful compared to a different way of determining cluster membership but

basis.mat = nmf_fit@W
coef.mat = nmf_fit@H


###########################################################
##         H-Matrix Analysis / cluster membership
###########################################################

cat("\n\n####################\nCluster Membership-- H-Matrix Analysis\n\n")

#### determine cluster membership ####
coef_norm  <- apply(coef.mat, 2, function(x) x/sum(x)) # normalize each column to sum of column
coef_consensus <- apply(coef_norm, 2, which.max) # which cluster had the maximum membership
coef_membership <- apply(coef_norm, 2, max) # maximum membership fraction
# mindiff core membership
coef_core_member <- apply(coef_norm, 2, function(x) {
  if ( sum(x == max(x)) > 1 ) return(FALSE) # if more than one cluster-score is the maximum, return FALSE
  min_diff_vec = max(x) - setdiff(x, max(x)) # difference between our largest coefficient, and our other coefficients 
  min_diff = min(min_diff_vec) # take the smallest difference
  min_diff > 1/as.numeric(opt$rank_top) # is this difference less than 1/rank?
})

NMF.annots = data.frame(row.names = colnames(coef.mat),
                        NMF.consensus = factor(paste0('C',NMF.consensus), # use predict results, rather than max-membership
                                               levels = paste0("C",1:opt$rank_top)), # also set levels now, so that clusters are always ordered correctly
                        NMF.cluster.membership = round(coef_membership, 3),
                        NMF.core.member = coef_core_member) %>%
  arrange(NMF.consensus, desc(NMF.core.member), # sort by cluster (low->high) & core-membership (core->noncore)
          desc(NMF.cluster.membership)) # then by cluster-membership (high->low)
#### save NMF membership info to CSV ####
write.table(NMF.annots, sep = "\t", # sorted by cluster-membership
            paste0(prefix, "clusterMembership.tsv"))

#### add NMF.annots to groups-file ####
groups.full = merge(groups, NMF.annots, by='row.names') %>% # append NMF.annots to groups file / cdesc
  arrange(NMF.consensus) %>%  # order by cluster-membership
  mutate_at(annots.is_continuous, as.numeric) %>% # mutate all continuous annotations to numerics
  column_to_rownames('Row.names') # re-add rownames

#### assign colors for NMF annotations ####
colors.full.NMF = set_annot_colors(NMF.annots, continuous.return_function = TRUE)
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
gct.H@rid <- paste0('C',1:nrow(gct.H@mat)) # use C# format for rid
rownames(gct.H@mat) = gct.H@rid # overwrite rownames with rid
gct.H@cid <- colnames(gct.H@mat)
write_gct(gct.H, ofile=paste0(prefix, 'H')) # write H to GCT
## normalized H-matrix
gct.H.norm = gct.H # initialize
gct.H.norm@mat <- apply(gct.H@mat, 2, function(x) x/max(x)) # normalized to the max of each column
write_gct(gct.H.norm, ofile=paste0(prefix, 'H_normToMax')) # write H-norm to GCT

cat("\n\n####################\nCluster Membership-- Coeficient Matrix Heatmap\n\n")
#### plot H.norm heatmap, with and without clustering ####
for (cluster_cols in c(TRUE,FALSE)) {
  # get heatmap row/col dimensions
  dim_rows = dim(gct.H.norm@mat)[1] + # get the number of rows in our heatmap
    dim(select(gct.H.norm@cdesc, -all_of(annots.excluded)))[2]  # get number of annotation rows in our heatmap
  dim_cols = dim(gct.H.norm@mat)[2] # get the number of columns in our heatmap
  
  hm = MyComplexHeatmap(gct.H.norm@mat, # plot normalized heatmap
                        heatmap_col = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), # colorscale for heatmap
                        show_heatmap_legend = T, show_column_names = T, # show heatmap & sample IDs
                        select(gct.H.norm@cdesc, -all_of(annots.excluded)), # with all annotations
                        c(colors, colors.NMF), show_annot=T, # using appropriate annotation colors
                        show_annotation_legend = F, # hide annotations
                        show_row_names = T, row_names_side = "left", # show row name
                        cluster_rows = F, # don't cluster rows
                        name = "Normalized Coefficient Value",
                        column_title = "NMF Coefficient Matrix (Normalized Sample-wise)", column_title_rot = 0, # add column title-- and specify rotation argument, bc otherwise R thinks I'm being lazy and writing column_title_rot
                        column_title_gp = gpar(fontsize = 16, fontface = "bold"), # format column title as if its a header
                        cluster_columns = cluster_cols, #optionally cluster columns
                        heatmap_legend_param = list(direction = "horizontal"))
  
  
  # use hm dimensions to calculate figure dimension
  width = max(dim_cols/5/4 + 3, # set width equal to n_columns/5/4 (column width = row height/4) + padding for row titles
              8) # (enforce minimum width, to prevent title cutoff)
  height = dim_rows/5 + 3  # set height equal to n_rows/5 + padding for title / legend
  # open PDF to write heatmap to
  pdf(paste0(prefix, 'heatmap_coeficientMatrixNorm', # save to PDF
             ifelse(cluster_cols, '_clustered', ''),'.pdf'),
      width, height)
  draw(hm, heatmap_legend_side='bottom', # put heatmap legend at bottom so it doesn't overlap annots
       padding = unit(c(2, 2, 2, 2), "mm")) # add padding so annotations render properly
  dev.off()
  # open PNF to write heatmap to
  png(paste0(prefix, 'heatmap_coeficientMatrixNorm', # save to PNG
             ifelse(cluster_cols, '_clustered', ''),'.png'),
      width=width, height=height, units='in', res=300)
  draw(hm, heatmap_legend_side='bottom', # put heatmap legend at bottom so it doesn't overlap annots
       padding = unit(c(2, 2, 2, 2), "mm")) # add padding so annotations render properly
  dev.off()
  
  # obselete pheatmap code
  # pheatmap::pheatmap(gct.H.norm@mat, color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), # plot normalized coefficient matrix
  #                    annotation_col = select(gct.H.norm@cdesc, -c(all_of(annots.excluded), "NMF.core.member")), # add all annotations
  #                    # consider: binary NMF.core.member column seems to break pheatmap...
  #                    annotation_colors = c(colors, colors.NMF), # add all annotation colors
  #                    show_rownames = TRUE,
  #                    cluster_cols = cluster_cols, cluster_rows = FALSE,
  #                    cellwidth=cw, cellheight=ch,
  #                    gaps_col = cumsum(table(gct.H.norm@cdesc$NMF.consensus)),
  #                    filename = paste0(prefix, 'heatmap_coeficientMatrixNorm', # save to PDF
  #                                      ifelse(cluster_cols, '_clustered', ''),'.pdf')) # with or without clustering suffix
}

#### draw Annotation Legend (standalone) to file ####
cat("\n\n####################\nGenerate Standalone-Legend\n\n")

#### pretty version ####
legend = lapply(lapply(hm@top_annotation@anno_list, # for each annotation
                       function(annot) {annot@color_mapping}), # get the color mapping object
                color_mapping_legend) %>% # turn them into a list of legends
  packLegend(list = ., # and then pack them into a legend
             max_height = unit(12, 'in')) # prevent from being too tall
legend_width = convertX(legend@grob$vp$width, unitTo="in", valueOnly=TRUE)
legend_height = convertX(legend@grob$vp$height, unitTo="in", valueOnly=TRUE)
# pdf
pdf(paste0(prefix, "annotationLegend.pdf"),
    width = legend_width+0.5, 
    height = legend_height+0.5) 
draw(legend)
dev.off()
# png
png(paste0(prefix, "annotationLegend.png"),
    width=legend_width+0.5,
    height=legend_height+0.5,
    units = 'in', res=300)
draw(legend)
dev.off()


#### scrollable version ####
legend_scrollable = lapply(lapply(hm@top_annotation@anno_list, # for each annotation
                                  function(annot) {annot@color_mapping}), # get the color mapping object
                           color_mapping_legend) %>% # turn them into a list of legends
  packLegend(list = ., # and then pack them into a legend
             max_width = unit(2.5, 'in')) # only allow one column
legend_width = convertX(legend_scrollable@grob$vp$width, unitTo="in", valueOnly=TRUE)
legend_height = convertX(legend_scrollable@grob$vp$height, unitTo="in", valueOnly=TRUE)
# pdf
pdf(paste0(prefix, "annotationLegend_scrollable.pdf"),
    width = legend_width+0.5, 
    height = legend_height+0.5) # height is full-height, plus padding
draw(legend_scrollable)
dev.off()


###########################################################
##         enrichment analysis
###########################################################

cat("\n\n####################\nEnrichement Analysis\n\n")

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

cat("\n\n####################\nSilhouette Plots\n\n")

#consider: only valid for non-bayesian?
rank.sil <- silhouette(res)
# create PDF version
pdf(paste0(prefix, 'silhouettePlot_clusterMembership.pdf')) # open PDF
plot(rank.sil, main=paste('K=', opt$rank_top, sep=''),
     col=colors.full.NMF$NMF.consensus$colors)
dev.off() # close PDF
# create PNG version
png(paste0(prefix, 'silhouettePlot_clusterMembership.png')) # open PNG
plot(rank.sil, main=paste('K=', opt$rank_top, sep=''),
     col=colors.full.NMF$NMF.consensus$colors)
dev.off() # close PNG



###########################################################
##         Boxplots of Continuous Features
###########################################################

cat("\n\n####################\nBoxplots of Continuous Annotation\n\n")

if ( length(annots.is_continuous) > 0 ) {
  pdf(paste0(prefix, 'boxplots_continuousAnnots.pdf')) # open PDF
  for (annot in annots.is_continuous) { # for each continuous annotation
    p = ggplot(groups.core, aes(x=NMF.consensus, y=!!sym(annot), # compare the values across clusters
                                color = NMF.consensus)) + scale_color_manual(values = colors.full.NMF$NMF.consensus$colors) +
      geom_boxplot() + # make a boxplot
      ggtitle(paste0("Values of ",annot," by Cluster")) +
      theme_minimal() + # give consistent theme
      theme(plot.title = element_text(face="bold")) + # increase title fontsize / bold 
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

cat("\n\n####################\nDriver Features-- W-Matrix Analysis\n\n")

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
cat("\n\n####################\nDriver Features-- T-Test\n\n")
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
    mutate(cluster = paste0('C',clust)) # add cluster 
  expr_mat = gct_expr_comb@mat # original feature expression
  ## subset expression matrix to driver-features, and add id column
  ft.expr = expr_mat[ft.vec, , drop=FALSE] %>% # get corresponding SIGNED feature expression
    as.data.frame() %>% rownames_to_column("feature.id") # convert to dataframe and turn rownames into a column
  ## create class vector assigning sample-columns to ingroup vs outgroup
  class.vec = colnames(expr_mat) %in% rownames(filter(NMF.annots, NMF.consensus==paste0('C',clust))) %>% # get TRUE/FALSE for whether a sample is in this cluster
    ifelse(., paste0('C',clust), "allother") # convert TRUE to cluster-label and FALSE to "allother" (using expr_mat since ft.expr has a feature-id column)
  names(class.vec) <- colnames(expr_mat) # add sample_id names back in (using expr_mat since ft.expr has a feature-id column)
  # note: it's important that the class.vec and expr_mat have columns IN THE SAME ORDER; names() are not actually taken into account...

  #### perform moderated test ####
  res.tmp <-  try({
    label = glue("{clust}.over.rest")
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


#### simplify t-test results into dataframe ####
driver.features.df = lapply(driver.features.list, function(df) {
  if (is.na(df)) return(NULL)
  df %>% # subset each data-frame to
    dplyr::select(c(!ends_with(".over.rest"), # annotation columns
                    starts_with("P.Value."))) %>% # and p-value column
    dplyr::rename("pval" = starts_with("P.value.")) # rename pvalue column for Easy Stacking
}) %>% do.call(rbind, .) %>% # rbind list into dataframe
  dplyr::mutate(!!paste0("signif_at_",opt$feature_fdr) := pval < opt$feature_fdr) # boolean for significance

#### Benjamini-Hochberg correction ####
driver.features.bh = mutate(driver.features.df,
                            bh.pval = p.adjust(driver.features.df$pval, method = 'BH'), # apply BH correction
                            !!paste0("bh.signif_at_",opt$feature_fdr) := bh.pval < opt$feature_fdr) # boolean if significant




#### Write Driver Feature + P-value to CSV ####
cat("\n\n####################\nDriver Features-- CSV & Excel Export\n\n")
## all features
write.csv(driver.features.bh,
          paste0(prefix, "driverFeatures.csv"))
## significant only
driver.features.sigFeatOnly = filter(driver.features.bh, #filter just to significant features
                                     !!sym(paste0("bh.signif_at_",opt$feature_fdr))) 
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
  driver.features.list.bh[[clust]]$global.adj.P.Value <- dplyr::filter(driver.features.bh, cluster==paste0('C',clust))$bh.pval # append a column with the BH corrected pvalues
  driver.features.list.bh[[clust]]$global.adj.P.Value_signif <- dplyr::filter(driver.features.bh, cluster==paste0('C',clust))$bh.signif_at # append a column with whether the value was significant
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
#### create Significant-Driver-Feature expresion-matrix heatmap (for each cluster) ####
cat("\n\n####################\nDriver Features-- Generate Expression Heatmaps\n\n")
for (clust in unique(driver.features.sigFeatOnly$cluster)) { # for each cluster that had significant features
  #### subset expression matrix to drivers for this cluster ####
  gct_expr_driver = subset_gct(gct_expr_comb,# take combined expression data
                               rid = filter(driver.features.sigFeatOnly, cluster==clust)$id) # filter to significant features for this cluster
  gct_expr_driver@cdesc = groups.full[colnames(gct_expr_comb@mat), , drop = FALSE] # overwrite rdesc with groups.full
  #### write out GCT ####
  write_gct(gct_expr_driver, ofile=paste0(prefix,'driverFeature_expression_',
                                          clust)) # write to GCT
  ##### append signficant features to heatmap ####
  if (clust == unique(driver.features.sigFeatOnly$cluster)[1]) { # if we're on the first cluster-heatmap 
    # initialize heatmap dimensions
    dim_rows = dim(gct_expr_driver@mat)[1]/5 + # set row dimensions to the number of driver features (scaled down so rows take up less space)
      dim(select(gct_expr_driver@cdesc, -all_of(annots.excluded)))[2]  # add rows for the number of annotation rows
    dim_cols = dim(gct_expr_comb@mat)[2] # set column dimensions to the number of samples
    # initialize hm.concat
    hm.concat <- MyComplexHeatmap(gct_expr_driver@mat, NULL, # plot expression of drivers, using default color-palette
                                  select(gct_expr_driver@cdesc, -all_of(annots.excluded)), # with all annotations
                                  c(colors, colors.NMF), # with custom annot-colors
                                  row_split = data.frame(ome_type = gct_expr_driver@rdesc$ome_type), # split rows by ome-type
                                  # row_title = "%s_%s",
                                  # show_row_names = T, row_names_side = "left",
                                  column_title = "Expression of Significant Driver Features", column_title_rot = 0, # add column title-- and specify rotation argument, bc otherwise R thinks I'm being lazy and writing column_title_rot
                                  column_title_gp = gpar(fontsize = 16, fontface = "bold"), # format column title as if its a header
                                  name = clust,
                                  cluster_rows = T) # cluster rows/features
  } else { # otherwise
    # update heatmap dimensions
    dim_rows = dim_rows + dim(gct_expr_driver@mat)[1]/5 # add n_driver_feature rows (scaled down so rows take up less space)
    # append new heatmap to hm.concat
    hm.concat <- hm.concat %v% # append heatmaps vertically
      MyComplexHeatmap(gct_expr_driver@mat, NULL, # plot expression of drivers, using default color-palette
                       select(gct_expr_driver@cdesc, -all_of(annots.excluded)), # with all annotations
                       c(colors, colors.NMF), # with custom annot-colors
                       show_annot=FALSE, # don't re-plot annotations
                       row_split = data.frame(ome_type = gct_expr_driver@rdesc$ome_type), # split rows by ome-type
                       # row_title = "%s_%s",
                       # show_row_names = T, row_names_side = "left",
                       name = clust,
                       cluster_rows = T) # cluster rows/features
  }
}
save(file = paste0(prefix, "driverFeatures_hm.concat.Rdata"), hm.concat) # save hm.concat object for troubleshooting

#### draw Heatmap to Files ####
cat("\n\n####################\nDriver Features-- Compile and Export Expression Heatmap\n\n")
width = min(max(dim_cols/5/4 + 3, # set width equal to n_columns/5/4 (column width = row height/4) + padding for row titles
                (dim_rows/5+2)/2, 8), # enforce a minimum width of at least half the height of the figure, or 8"
            12) # limit  width to at most 12 inches
height = min(dim_rows/5 + 3,  # set height equal to # rows / 10, plus padding for title
             24) # limit  width to at most 32 inches
# open PDF to write heatmap to
pdf(paste0(prefix, "driverFeatures_complexHeatmap.pdf"),
    width,
    height)  
# pdf('/opt/input/outputs/tmp.pdf', 12, 12) # for manual docker testing
draw(hm.concat, #heatmap_legend_side='right',
     show_heatmap_legend = FALSE,
     ht_gap = unit(5, "mm")) #row_title = "Driver Features (grouped by Cluster)")
dev.off()

png(paste0(prefix, "driverFeatures_complexHeatmap.png"),
    width = width, 
    height = height,
    units="in", res=300)
draw(hm.concat, #heatmap_legend_side='right',
     show_heatmap_legend = FALSE)
dev.off()



###########################################################
#### create UpSet Plot of driver-feature-overlap b/t clusters ####
cat("\n\n####################\nDriver Features-- Upset Plots\n\n")
if ( length(unique(driver.features.bh$cluster)) > 1 ) { # if more than one cluster has driver features
  upset.mat <- driver.features.sigFeatOnly %>% # use only significant features
    dplyr::select(c(id, cluster)) %>% # select relevant columns
    dplyr::mutate(observed = 1) %>% # make placeholder "observed" column, to use as the value in the pivot_wider() function
    pivot_wider(id_col = 'id', # use feature-id as id
                names_from='cluster', # create column for each cluster with driver features
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
cat("\n\n####################\nDriver Features-- # per Cluster (per Ome-Type)\n\n")
ft_countByOme <- driver.features.sigFeatOnly %>% # use only significant features
  dplyr::select(c(ome_type, cluster)) %>% # select relevant columns
  group_by(cluster, ome_type) %>% summarise(n_feat = n()) %>% # get number of features per-ome per-cluster
  ungroup() %>% complete(ome_type, cluster, fill = list(n_feat = 0)) # add rows for ome_type/cluster combinations with n_feat = 0
#### write CSV of contributions per ome ####
write.csv(pivot_wider(ft_countByOme, id_cols=ome_type, # pivot ft_countByOme wider for easer viewing / interpretation
                      names_from=cluster, values_from=n_feat), # make a column for each cluster, with values being the number of features in that cluster/ome combo
          paste0(prefix, "driverFeatures_contributionsByOme.csv"))
#### plot Barplot of Significant Features per Ome ####
p <- ggplot(ft_countByOme, aes(x = cluster, y = n_feat, fill=ome_type)) +
  geom_bar(position="dodge", stat="identity") + # grouped barplot with n_feat per clut per ome
  facet_grid(~ cluster, space="free_x", scales="free_x", switch="x") + # facet by cluster to make the deliniation clearer
  scale_x_discrete(labels=NULL) + # remove x-axis labels so they aren't double-plotted
  geom_text(aes(label = n_feat), # add labels with number of features
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) + # position label above each column
  labs(x = "NMF Basis", y = "Number of Features") + # x axis and y axis labels
  ggtitle('Significant Driver-Features by Cluster') + # add title
  theme_minimal() +  # adjust formatting
  theme(panel.grid.major.x = element_blank(),  # remove x-axis gridlines
        plot.title = element_text(face="bold"), # increase title fontsize / bold
        strip.background = element_rect(fill=NA,colour="grey80")) # add border around cluster-label

# ggsave('/opt/input/tmp.pdf', # for manual docker testing
#        plot = p, device="pdf") # plot as PDF
ggsave(paste0(prefix, "driverFeatures_contributionsByOme_barplot.pdf"),
       plot = p, device="pdf") # plot as PDF
ggsave(paste0(prefix, "driverFeatures_contributionsByOme_barplot.png"),
       plot = p, device="png") # plot as PNG





###########################################################
#### Top Driver-Features Expression by Cluster (Boxplot & Heatmaps)
###########################################################
#### subset driver-features to top N features ####
driver.features.topNFeat <- driver.features.sigFeatOnly %>% # take significant driver-features
  group_by(cluster) %>% # group by cluster
  slice_min(order_by = bh.pval, # order by p-value
            n = opt$top_n_features) %>% # select n features with lowest p-values
  ungroup() # ungroup


#### make Boxplot PDF for each cluster. ####
cat(glue("\n\n####################\nDriver Features-- Top {opt$top_n_features} Expression Boxplots\n\n"))
for (cluster in unique(driver.features.topNFeat$cluster)) {
  pdf(paste0(prefix, "driverFeatures_expressionBoxplots_", cluster,".pdf"))
  # pdf(glue("/opt/input/tmp_{cluster}.pdf")) # for manual docker testing
  driverFeats = filter(driver.features.topNFeat, cluster==!!cluster)$id # get vector of driver features we wanna plot 
  for (ft in driverFeats) {
    #### collect the info we need for our boxplot ####
    expr.df = data.frame(ft.expression = gct_expr_comb@mat[ft,], # get expression data
                         # sample.id = gct_expr_comb@cid, # get corresponding sample name (unnecessary but convenient)
                         cluster = NMF.annots[gct_expr_comb@cid, # in the order of the samples
                                              'NMF.consensus']) #get cluster membership
    annot.df = filter(driver.features.sigFeatOnly, id==ft) # filter annotations to relevant driver-feature
    
    # boxplot title
    if (!is.null(opt$gene_col)) { # if we have the gene column, use it for the boxplot title
      boxplot.title = glue('Expression of {annot.df[[opt$gene_col]]} in {annot.df$ome_type} ({annot.df$id_og})')
    } else { # otherwise just use the feature ID
      boxplot.title = glue('Expression of feature {annot.df$id_og} ({annot.df$ome_type})')
    }
    
    # boxplot colors
    boxplot.colors = rep('black', opt$rank_top) #initialize color-scheme with all black
    names(boxplot.colors) = paste0("C", 1:opt$rank_top) # name vector w/ clusters
    boxplot.colors[[cluster]] = "red" # color the expression red in the cluster the driver feature drives
    
    #### create boxplot ####
    p <- ggplot(expr.df, aes(x = cluster, y = ft.expression)) + # boxplot with ft-expression by ome
      geom_boxplot(color = boxplot.colors, # boxplots of expression data
                   width = 0.5) + 
      geom_hline(yintercept = 0, color='grey', linetype='dashed') + # add dashed line at 0
      labs(x = "NMF Basis", y = "Feature Expression", # x axis and y axis labels
           title = boxplot.title, # add title
           subtitle = glue("Driver Feature in {annot.df$cluster}\nFeature Score: {signif(annot.df$feature.scores, 3)}\nAdj-P-Val: {signif(annot.df$bh.pval, 3)}")) +
      stat_summary(fun.data = function(y) {c(y = max(y), # add label slightly above box plot
                                             label = length(y))}, # label with the number of observations
                   geom = "text", vjust = -1,
                   position = position_dodge(width = 0.75)) +
      theme_minimal() + # adjust formatting
      theme(panel.grid.major.x = element_blank(), # remove lines
            plot.title = element_text(face="bold")) # increase title fontsize / bold
    plot(p) # plot as PDF page
  }
  dev.off()
}



#### make Heatmap PDF for each cluster ####
cat(glue("\n\n####################\nDriver Features-- Top {opt$top_n_features} Expression Heatmaps\n\n"))
for (cluster in unique(driver.features.topNFeat$cluster)) {
  driverFeats = filter(driver.features.topNFeat, cluster==!!cluster)$id # get vector of driver features we wanna plot 
  
  if (!is.null(opt$gene_col)) { # if we have the gene column
    goi = filter(driver.features.sigFeatOnly, id %in% driverFeats) %>%
      .[[opt$gene_col]] %>% unique() # and identify unique genes
    ft_list = lapply(goi, function(gene) {
      filter(driver.features.sigFeatOnly, !!as.name(opt$gene_col)==gene) %>% # filter for any significant drivers from the gene-of-interest
        .$id}) # get feature ID
  } else { # otherwise
    goi = NULL
    ft_list = as.list(driverFeats) # just use the driverFeats as is
  }
  
  # calculate heatmap dimensions
  dim_rows = 4 + # add space for 4 extra rows per dataset, as an estimate
    dim(select(groups.full[colnames(gct_expr_comb@mat), , drop = FALSE], -all_of(annots.excluded)))[2]  # add rows for the number of annotation rows
  dim_cols = dim(gct_expr_comb@mat)[2] # set column dimensions to the number of samples
  # open PDF to draw heatmap
  pdf(paste0(prefix, "driverFeatures_expressionHeatmaps_", cluster,".pdf"),
      max(dim_cols/5/4 + 1.5, # set width equal to n_columns/5/4 (column width = row height/4) + padding for row titles
          8), # (enforce minimum width, to prevent title cutoff)
      dim_rows/5 + 2) # set height equal to n_rows/10 x scaling + padding for title / legend
  # pdf('/opt/input/tmp.pdf', 12,8) # for manual docker testing
  for (i in 1:length(ft_list)) {
    fts = ft_list[[i]]
    
    if (!is.null(opt$gene_col)) { # if we have the gene column, use it for the boxplot title
      hm.title = glue('Expression of Driver Features for Gene {goi[i]}')
    } else { # otherwise just use the feature ID
      hm.title = glue('Expression of Driver Feature {fts}')
    }
    
    hm = MyComplexHeatmap(gct_expr_comb@mat[fts,,drop=F], NULL, # plot expression of drivers, using default color-palette
                          select(groups.full[colnames(gct_expr_comb@mat), , drop = FALSE],
                                 -all_of(annots.excluded)), # with all annotations
                          c(colors, colors.NMF), # with custom annot-colors
                          show_annot=TRUE, show_annotation_legend = F,
                          row_dend_side = 'left',
                          show_row_names = TRUE, row_names_side = "right",
                          column_title = hm.title, column_title_rot = 0, # add column title-- and specify rotation argument, bc otherwise R thinks I'm being lazy and writing column_title_rot
                          column_title_gp = gpar(fontsize = 16, fontface = "bold")) # format column title as if its a header
    draw(hm, padding = unit(c(2, 2, 2, 2), "mm"))
  }
  dev.off()
}


###########################################################
##         TSNE Comparison of Clusters
###########################################################

#### tSNE Clustering on NMF Core Samples (across all Features) ####
cat(glue("\n\n####################\nTSNE Core-Sample Comparisons\n\n"))
coreSample_expr <- gct_expr_comb@mat[, rownames(filter(NMF.annots,  # expression matrix
                                                       NMF.core.member))] %>% # of JUST core features
  t() # samples as rows
tsne_s <- tryCatch(Rtsne(coreSample_expr, check_duplicates=F), # try tSNE as is
                   # but if it fails
                   error = function(e) { #try again with perplexity set explicitly
                     message(conditionMessage(e))
                     message("Trying again, with perplexity set explitictly.")
                     try(Rtsne(coreSample_expr,
                               perplexity = floor((nrow(coreSample_expr) - 1) / 3), # to avoid perplexity error
                               check_duplicates=F))
                   })

#### plot results ####
if (!class(tsne_s) == 'try-error') { # assuming tSNE ran properly
  #### create dataframe for plotting ####
  tsne.df <- data.frame(tsne_x = tsne_s$Y[,1],
                        tsne_y = tsne_s$Y[,2],
                        sample_id = rownames(coreSample_expr),
                        cluster = NMF.annots[rownames(coreSample_expr),
                                             "NMF.consensus"])
  #### plot and save ####
  p = ggplot(tsne.df, aes(tsne_x, tsne_y, color = cluster)) +
    geom_point() +
    scale_color_manual(values = colors.NMF$NMF.consensus) + # color by NMF scheme
    labs(title = "tSNE Clustering on Core NMF Samples", # title
         subtitle = "Clustered across all features") + # subtitle, describing the dataset
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + # plot 0,0 axes in black
    theme_minimal() +
    theme(plot.title = element_text(face="bold")) # increase title fontsize / bold
  # ggsave("/opt/input/tmp.png", plot=p, device="png") # for manual docker testing
  ggsave(paste0(prefix, "tSNE_coreSamples.pdf"),
         plot=p, device = "pdf") # saveas pdf
  ggsave(paste0(prefix, "tSNE_coreSamples.png"),
         plot=p, device = "png") # as png
} else {
  cat("\nWARNING: TSNE clustering failed on NMF core samples.") # try again on driver-features
}


#### tSNE Clustering on significant Driver Features (across all Samples) ####
cat(glue("\n\n####################\nTSNE Driver-Feature Comparisons\n\n"))
driverFeat_expr <- gct_expr_comb@mat[driver.features.sigFeatOnly$id,] # driver features only, all samples
tsne_f <- tryCatch(Rtsne(driverFeat_expr, check_duplicates=F), # try tSNE as is
                   # but if it fails
                   error = function(e) { #try again with perplexity set explicitly
                     message(conditionMessage(e))
                     message("Trying again, with perplexity set explitictly.")
                     try(Rtsne(driverFeat_expr,
                               perplexity = floor((nrow(driverFeat_expr) - 1) / 3), # to avoid perplexity error
                               check_duplicates=F))
                   })
#### plot results ####
if (!class(tsne_f) == 'try-error') { # assuming tSNE ran properly
  #### create dataframe for plotting ####
  feat_id = rownames(driverFeat_expr) # get feature-ids from expression matrix
  cluster = driver.features.sigFeatOnly[match(feat_id, # match order of feature-ids
                                              driver.features.sigFeatOnly$id), # to the driver feature dataframe
                                        "cluster"] %>%
    factor(levels = paste0("C",1:opt$rank_top))
  tsne.df <- data.frame(tsne_x = tsne_f$Y[,1],
                        tsne_y = tsne_f$Y[,2],
                        feat_id = feat_id,
                        cluster = cluster) # pull out corresponding cluster
  #### plot and save ####
  p = ggplot(tsne.df, aes(tsne_x, tsne_y, color = cluster)) +
    geom_point() +
    scale_color_manual(values = colors.NMF$NMF.consensus) + # color by NMF scheme
    labs(title = "tSNE Clustering on significant Driver Features", # title
         subtitle = "Clustered across all samples") + # subtitle, describing the dataset
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + # plot 0,0 axes in black
    theme_minimal() +
    theme(plot.title = element_text(face="bold")) # increase title fontsize / bold
  # ggsave("/opt/input/tmp.png", plot=p, device="png") # for manual docker testing
  ggsave(paste0(prefix, "tSNE_driverFeatures.pdf"),
         plot=p, device = "pdf") # saveas pdf
  ggsave(paste0(prefix, "tSNE_driverFeatures.png"),
         plot=p, device = "png") # as png
} else {
  cat("\nWARNING: TSNE clustering failed on Driver Features.") # try again on driver-features
}


###########################################################
##         Tar all Files
###########################################################
cat(glue("\n\n####################\nTarring Files\n\n"))

output.files = c(list.files(pattern=prefix, full.names = T), # grab all figures
                 list.files("nmf_postprocess_opt.Rdata", full.names = T)) # also grab the opt parameters object
tar(paste0(opt$output_prefix,'_NMF_postprocess.tar.gz'),
    files = output.files, compression = "gzip", tar="tar")


