#!/usr/bin/env Rscript
#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  # make_option( c("-t", "--tar"), action='store', type='character',  dest='tar_file', help='tar file containing data tables in GCT v1.3 format.'),
  make_option( c("-d", "--ome_gcts"), action='store', type='character', dest='ome_gcts_string', help='String of GCT-files, comma-separated.') , ## default=0),
  make_option( c("-o", "--ome_labels"), action='store', type='character', dest='ome_labels_string', help='String of labels corresponding to GCT-files, comma-separated.') , ## default=0),
  #### Pre-process Parameters ####
  make_option( c("-f", "--sd_filt_min"), action='store', type='numeric', dest='sd_filt_min', help='Lowest standard deviation across columns percentile to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest sd will be removed. Will be applied before z-scoring (optional)') , ## default=0),
  make_option( c("-g", "--sd_filt_mode"), action='store', type='character', dest='sd_filt_mode', help='Determines how the low variable features (paramater "sd_filt_mode") will be removed from the dataset. "global": the filter will be applied to the entrire datasets. "separate": the will be applied separately to each data type (e.g. proteome, RNA, etc.). "equal": all data tables of different data types (e.g. proteome, RNA) will be filtered down to have the same number of features and ths is limited by the data type with the smallest number of features "N_feat". The data tables of other data types are filtered to retain the "N_feat" highest variable features.  Only relevant for multi-omics data.'), 
  make_option( c("-v", "--z_score_mode"), action='store', type='character', dest='z_score_mode', help='z-scoring mode: row, col, rowcol'), # default=TRUE),
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  #### NMF Parameters ####
  make_option( c("--kmin"), action='store', type='numeric',  dest='kmin', help='Minimal factorization rank.'),  # default = 2),
  make_option( c("--kmax"), action='store', type='numeric',  dest='kmax', help='Maximal factorization rank.'), # default = 4),
  make_option( c("-n", "--nrun"), action='store', type='numeric',  dest='nrun', help='Number of NMF runs with different starting seeds.'), # default = 5),
  make_option( c("-b", "--bayesian"), action='store', type='logical',  dest='bnmf', help='If TRUE Bayesian NMF to determine optimal rank will be used (package ccfindR).', default = FALSE),
  make_option( c("-m", "--nmf_method"), action='store', type='character',  dest='nmf_method', help='NMF Method.'),
  make_option( c("-s", "--seed"), action='store', type='character',  dest='seed', help='Seed method for NMF factorization.'), # default = 'random'),
  # make_option( c("-e", "--exclude_2"), action="store", dest='exclude_2', type="logical", help="If TRUE, a '2' will be excluded from calculation of the optimal rank."), #default=TRUE),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters.", default=NA),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.")
  #### ####
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # # for testing arguments
                       # args = c('-d',"/opt/input/LSCC_GCP_AllPlates_tumor_nat_n353x74.gct,/opt/input/LSCC_GCP_AllPlates_NATnorm_n155x74.gct",
                       #          '-o',"tumorOnly,natNorm",
                       #          '-a',"geneSymbol",
                       #          '-n','50', '-m','lee',
                       #          '--kmin','3', '--kmax','5',
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

library(parallel) # for detectCores()



###########################################################
##
##         Pre-Processing GCT Files
##
###########################################################

###################################################
##               Import Data
###################################################

import_gct_from_wdl <- function( ome_gcts_string, ome_labels_string ){
  
  # parse strings into vectors
  ome_gct_str = unlist(strsplit(ome_gcts_string, ","))
  ome_labels = unlist(strsplit(ome_labels_string, ","))
  # ome_labels = sapply(ome_gct_str, USE.NAMES = FALSE, function(filepath) { basename(filepath) %>% gsub(".gct","", .) }) # assume that file-name is ome-output_prefix
  
  if (length(ome_gct_str)!=length(ome_labels)) stop("Different number of GCT files and ome-labels detected. Please check your inputs.")
  
  ## import
  ome_gcts <- lapply(ome_gct_str, parse_gctx)
  names(ome_gcts) <- ome_labels
  
  out <- list(ome_gcts=ome_gcts, ome_labels=ome_labels)
  
  return(out)
}

gct_imp = import_gct_from_wdl( opt$ome_gcts_string, opt$ome_labels_string)
ome_gcts <- gct_imp$ome_gcts # list of GCT objects
ome_labels <- gct_imp$ome_labels # names associated with those GCTs

###################################################
##               Merge GCT Components
###################################################
for (ome in ome_labels) {
  if (ome==ome_labels[[1]]) { # if we're on our first iteration
    #### initialize each component of the GCT ####
    comb_cdesc = ome_gcts[[ome]]@cdesc # initialize cdesc
    comb_rdesc = data.frame(id = paste(ome, ome_gcts[[ome]]@rid, sep="_"), # initialize rdesc with id column
                            id_og = ome_gcts[[ome]]@rid, # initialize rdesc with id column
                            ome_type = ome) %>% # and ome
      mutate(!!opt$gene_col := ome_gcts[[ome]]@rdesc[[opt$gene_col]]) # add opt$gene_col column from rdesc
    comb_mat_raw = ome_gcts[[ome]]@mat # initialize matrix
  } else { # after the first instance
    #### merge each component of the GCT to the existing component ####
    # merge cdescs, prioritizing first entry
    comb_cdesc = merge(comb_cdesc, ome_gcts[[ome]]@cdesc, all=TRUE) %>% # merge cdescs across all shared columns, keeping all data
      distinct(id, .keep_all = TRUE) # filter out duplicated CIDs based on the id column
    # append new rdesc
    comb_rdesc = rbind(comb_rdesc, # append new entries to the end of the combined rdesc
                       data.frame(id = paste(ome, ome_gcts[[ome]]@rid, sep="_"), # initialize rdesc with id column
                                  id_og = ome_gcts[[ome]]@rid, # initialize rdesc with id column
                                  ome_type = ome) %>% # ome-output_prefix
                         mutate(!!opt$gene_col := ome_gcts[[ome]]@rdesc[[opt$gene_col]])) # add opt$gene_col column from rdesc # and ome
    # append matrices, based on shared samples / CIDs
    shared_samples <- intersect(colnames(comb_mat_raw), colnames(ome_gcts[[ome]]@mat))
    if (length(shared_samples)==0) stop(paste("The GCT for",ome,"does not share any samples in common with previously loaded GCTs"))
    comb_mat_raw = rbind(comb_mat_raw[,shared_samples],
                     ome_gcts[[ome]]@mat[,shared_samples]) # append new entries to the end of the combined matrix
    # overwrite matrix rownames with the fully-unique IDs
    rownames(comb_mat_raw) = comb_rdesc$id
  }
}

###################################################
##               Preprocess Matrix
###################################################

comb_mat = comb_mat_raw  # initialize from unprocessed-combined-matrix

#### remove features with missing values ####
keep_idx <- which(apply(comb_mat, 1, function(x) sum(is.na(x)))==0) # identify features that have no missing values
comb_mat <- comb_mat[keep_idx,] # subset matrix to fully quantified features

#### apply SD filter ####
mat_s = comb_mat # initialize
if(!is_null(opt$sd_filt_mode)){
  # todo: finish SD filter
  print("SD Filter has not been set up yet!")
}
comb_mat = mat_s # overwrite with (optionally) sd-filtered matrix

#### z-score multi-omic data matrix ####
mat_z = comb_mat # initialize
if(!is_null(opt$z_score_mode)){
  if(z_score_mode %in% c("row", "rowcol")) # if we are row-normalizing
    mat_z <- apply(mat_z, 1, function(x) (x-mean(x))/sd(x)) %>% # z-score rows
      t() # untransform matrix
  if(z_score_mode %in% c("col", "rowcol")) # if we are column normalizing
    mat_z <- apply(mat_z, 2, function(x) (x-mean(x))/sd(x)) # z-score cols (do samples AFTER feature-wise z-score)
}
comb_mat = mat_z # overwrite with (optionally) z-scored matrix

###################################################
##               Create Combined GCT Object
###################################################
gct_comb = GCT(mat = comb_mat, # fully quantified, optionally z-scored matrix
               rdesc = comb_rdesc[match(rownames(comb_mat), comb_rdesc$id), ],
               rid = rownames(comb_mat),
               cdesc = comb_cdesc[match(colnames(comb_mat), comb_cdesc$id), ], # pull the cdesc entries in the order they appear in the matrix
               cid = colnames(comb_mat))

write_gct(gct_comb, ofile=paste0(opt$output_prefix, "_combined"))

###################################################
##               Create Non-Negative GCT Object
###################################################
#### get all positive values ####
mat_pos = gct_comb@mat # initialize positive matrix
mat_pos[ mat_pos < 0 ] = 0 # replace all negative values with zero
rownames(mat_pos) = paste(rownames(mat_pos), 'pos', sep='_')
# create rdesc with sign-info
rdesc_pos = gct_comb@rdesc %>%
  mutate(sign = "pos", # save sign
         id_unsigned = id, # save unsigned ID to new column
         id = rownames(mat_pos)) # save unique (signed) ID to gct_comb

#### get all negative values ####
mat_neg = gct_comb@mat # initialize
mat_neg[ mat_neg > 0 ] = 0 # replace all positive values with zero
mat_neg = abs(mat_neg) # take abs() to get positives
rownames(mat_neg) = paste(rownames(mat_neg), 'neg', sep='_')
# create rdesc with sign-info
rdesc_neg = gct_comb@rdesc %>%
  mutate(sign = "neg", # save sign
         id_unsigned = id, # save unsigned ID to new column
         id = rownames(mat_neg)) # save unique (signed) ID to gct_comb

#### merge pos/neg into matrix ####
gct_NN = GCT(mat = rbind(mat_pos, mat_neg), # rbind matrixes
             rdesc = rbind(rdesc_pos, rdesc_neg), # rbind rdescs
             rid = rownames(rbind(mat_pos, mat_neg)), # rownames of rbinded matrices
             cdesc = gct_comb@cdesc, # same cdesc
             cid = gct_comb@cid) # same CID

#### remove features that only have zero (i.e. remove _neg entry for entirely positive feature, and visa versa) ####
keep.id <- names(which(apply(gct_NN@mat, 1, function(x) sum(x != 0) ) > 0)) # filter only for features with at least one non-zero value
gct_NN_fin = subset_gct(gct_NN, rid = keep.id)
#### write GCT file ####
write_gct(gct_NN_fin, ofile=paste0(opt$output_prefix, "_combinedNonNegative"))

###########################################################
##
##         Run NMF Analysis
##
###########################################################

#### Run NMF ####
res.rank = list() # initialize list for NMF results
if(opt$bnmf){
  print("This is a placeholder for Bayesian NMF. Someone sure should write this code!")
} else { # otherwise run normal NMF
  ranks = sort(opt$kmin:opt$kmax) # vector of ranks to try
  cores <- detectCores()
  opts <- paste('vp', cores,'t', sep='') # cores...?
  for(rank in ranks) # for each rank
    res.rank[[rank]] <-  nmf(gct_NN_fin@mat, rank=rank, method=opt$nmf_method, seed=opt$seed, nrun=opt$nrun, .options = opts)
}

#### Calculate Cluster Metrics ####
if (length(ranks) > 1) { # if we have more than one rank
  #### Calculate Cluster Metrics ####
  ## silhouette
  rank.sil <- lapply(res.rank[ranks], silhouette)
  rank.sil.avg <- lapply(rank.sil, function(x) tapply( x[,3], x[, 1], mean))
  
  ## cophenetic
  rank.coph <- sapply(res.rank[ranks], cophcor)
  
  ## dispersion of consensus matrix
  rank.disp <- sapply(res.rank[ranks], dispersion)
  
  ## combine
  rank.coph.disp <- rank.disp^(1-rank.coph)
  
  #### Determine Best Rank ####
  
  # it looks like we want the lowest number of clusters
  # so we say "if the next clustering is better, keep going"
  # else if it's the same or greater
  # # if it's within the relative increase threshold, keep going
  # # else declare current best our top rank
  
  # todo: figure out if I want to get rank a more.... elegant? way?
  top.rank = ranks[[which.min(rank.coph.disp)]] # for now literally just take the min

} else {top.rank = ranks}


#### Save NMF Object / Top Rank ####
save(file = "nmf_res.Rdata", res.rank)
write.table(top.rank, file = "nmf_best_rank.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
