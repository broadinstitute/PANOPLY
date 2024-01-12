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
  make_option( c("-u", "--z_score"), action='store', type='logical', dest='z_score', help='Should the data be z-scored?'), # default=TRUE),
  make_option( c("-v", "--z_score_mode"), action='store', type='character', dest='z_score_mode', help='z-scoring mode: row, col, rowcol'), # default=TRUE),
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  make_option( c("-i", "--org_id"), action='store', type='character', dest='organism', help='Organism type, used for gene-mapping if gene_col is provided. Support for \'human\' (Hs), \'mouse\' (Mm), or \'rat\' (Rn).'), # default='human'),
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
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters."),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.")
  #### ####
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # for testing arguments
                       # args = c('-d',"/opt/input/acetylome-filtered_table-output.gct,/opt/input/phosphoproteome-filtered_table-output.gct,/opt/input/proteome-filtered_table-output.gct,/opt/input/ubiquitylome-filtered_table-output.gct",
                       #          '-o',"acK,pSTY,prot,ubK",
                       #          '-y','/opt/input/master-parameters.yaml',
                       #          '-x',"odg_test")
                       )



#### Load Libraries ####
library(cmapR)
library(NMF)
library(tidyverse)
library(glue)

library(parallel) # for detectCores()



#### Parse YAML Arguments ####
opt = opt_cmd # initialize options with command line options
if ( !is.null(opt$yaml_file) ) {
  #### read in yaml ####
  library(yaml)
  yaml_out <- read_yaml(opt$yaml_file)
  #### locate the NMF parameters section ####
  if ( !is.null(yaml_out$panoply_nmf) ) { # if we have an NMF parameters section
    yaml_nmf =  yaml_out$panoply_nmf # read in those parameters
  } else { # if the section is missing
    if ( !is.null(yaml_out$panoply_mo_nmf) ) { # check for the deprecated mo_nmf section
      warning(glue("The parameter file '{opt$yaml_file}' is missing the 'panoply_nmf' section, but includes the deprecated 'panoply_mo_nmf' parameter section. Parameters will be read in from 'panoply_mo_nmf', but please consider updating your parameters file!"))
      yaml_nmf =  yaml_out$panoply_mo_nmf # read in those parameters
    } else { # otherwise, stop
      stop(glue("The parameter file '{opt$yaml_file}' does not contain an NMF parameters section. Please check that the yaml file contains the appropraite 'panoply_nmf' section."))
    }
  }
  #### overwrite the command-line parameters ####
  # global parameters
  if (is.null(opt$gene_col)) opt$gene_col = yaml_out$global_parameters$gene_mapping$gene_id_col
  if (is.null(opt$organism)) opt$organism = yaml_out$global_parameters$organism
  # nmf parameters
  if (is.null(opt$kmin)) opt$kmin = yaml_nmf$kmin
  if (is.null(opt$kmax)) opt$kmax = yaml_nmf$kmax
  if (is.null(opt$seed)) opt$seed = yaml_nmf$seed
  if (is.null(opt$nmf_method)) opt$nmf_method = yaml_nmf$method
  # if (is.null(opt$exclude_2)) opt$exclude_2 = yaml_nmf$exclude_2
  if (is.null(opt$nrun)) opt$nrun = yaml_nmf$nrun
  if (is.null(opt$sd_filt_min)) opt$sd_filt_min = yaml_nmf$sd_filt
  if (is.null(opt$sd_filt_mode)) opt$sd_filt_mode = yaml_nmf$filt_mode
  if (is.null(opt$z_score)) opt$z_score = yaml_nmf$z_score
  if (is.null(opt$z_score_mode)) opt$z_score_mode = yaml_nmf$z_score_mode
} 

# print parameters
cat("### NMF PARAMETERS ###\n")
print(opt)
cat("\n#####\n")

###########################################################
##
##         Pre-Processing GCT Files
##
###########################################################

###################################################
##               Import Data
###################################################

import.gct.from.wdl <- function( ome_gcts_string, ome_labels_string ){
  
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

gct_imp = import.gct.from.wdl( opt$ome_gcts_string, opt$ome_labels_string)
ome_gcts <- gct_imp$ome_gcts # list of GCT objects
ome_labels <- gct_imp$ome_labels # names associated with those GCTs

###################################################
##               Merge GCT Components
###################################################
#### Merge GCT Components ####
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
#### Add Gene Annotations to rdesc (requires gene column & organism-ID) ####
if( !is.null(opt$gene_col) & !is.null(opt$organism) ){ # if we have the required info to add gene-annotations
  # convert opt$organism to org_id
  org_id = recode(opt$organism,
                  "human" = "Hs",
                  "mouse" = "Mm",
                  "rat" = "Rn")
  # add annotations
  if ( is.null(comb_rdesc[[opt$gene_col]]) ) { # sanity check that gene-column is IN the rdesc
    warning(paste0("Gene Symbol Column '", opt$gene_col,
                   "' was not found in the rdesc. Gene-annotations cannot be added."))
    } else {  # if the gene_column is in the rdesc
      # load AnnotationDBI & organism package
      library(AnnotationDbi)
      library(paste0("org.",org_id,".eg.db"), character.only = TRUE) # load appropraite organism database
      #### add description and enzyme codes ####
      genes = unique(comb_rdesc[[opt$gene_col]]) # get unique gene-symbols to annotate
      gene_annot_df <- AnnotationDbi::select(eval(sym(paste0("org.",org_id,".eg.db"))),
                                             keys=genes, keytype='SYMBOL',
                                             column=c( 'GENENAME',  'ENZYME', 'ENTREZID'),
                                             multiVals='first')
      
      #### add cytoband annotations ####
      if (org_id=="Hs") { # if we're working humans
        ## add cytoband ID
        map <- unlist(as.list(org.Hs.egMAP[mappedkeys(org.Hs.egMAP)])) # get maping between entrez ID and cytoband ID; convert from S4 to vector by unlist(as.list()) 
        gene_annot_df$CYTOBAND <- map[ gene_annot_df$ENTREZID ] # map entrez ID onto cytoband map
        # consider: may need to have support for NULL entrez.id values
      }
      
      #### add annotations to rdesc ####
      comb_rdesc = cbind(comb_rdesc, # original rdesc
                         gene_annot_df[match(comb_rdesc[[opt$gene_col]], # sort gene_annot_df, to match comb_rdesc order
                                             gene_annot_df$SYMBOL), ])
    }
}

###################################################
##               Preprocess Matrix
###################################################

comb_mat = comb_mat_raw  # initialize from unprocessed-combined-matrix

#### remove features with missing values ####
cat("\n####################\nFilter Missing Values\n")
keep_idx <- which(apply(comb_mat, 1, function(x) sum(is.na(x)))==0) # identify features that have no missing values
# print results & overwrite comb_mat
cat(glue("\nRemoved {dim(comb_mat)[1]-length(keep_idx)} (of {dim(comb_mat)[1]}) features with missing values, leaving {length(keep_idx)} fully-quantified features.\n\n"))
comb_mat <- comb_mat[keep_idx,] # subset matrix to fully quantified features

#### apply SD filter ####
mat_s = comb_mat # initialize
if(!is_null(opt$sd_filt_mode)){
  cat("\n####################\nFilter Standard Deviation\n")
  #### sd-filter utility function ####
  sd_filter_featwise <- function(mat_s, sd_filt_min) { # features-wise filter to sd_filt_min percentile
    if (sd_filt_min==0) return(mat_s) # don't filer if sd_filt_min is 0
    
    feat_sd <- apply(mat_s, 1, sd) # calculate SD across all features in matrix
    sd_cutoff <- quantile(feat_sd, c(sd_filt_min)) # determine value for quantile-cutoff
    idx.keep <- which( feat_sd > sd_cutoff) # get ID of features to keep
    
    return(mat_s[idx.keep,]) # filter matrix to features that are above SD threshold
  }
  
  #### apply sd-filter ####
  if(opt$sd_filt_mode == 'global'){ # filter all features across all omes
    #### apply global filter ####
    mat_s = sd_filter_featwise(mat_s, opt$sd_filt_min)
  } else if (opt$sd_filt_mode == 'separate') { # filter each ome individually
    #### apply separate filters ####
    idx.keep = c() # initialize vector for features to keep
    for (ome in unique(comb_rdesc$ome_type)) {
      feat_byOme = filter(comb_rdesc, ome_type==ome)$id # get features-ids that correspond to this ome
      mat_tmp = mat_s[intersect(feat_byOme,rownames(mat_s)),] # subset expression-matrix to these features
      idx.keep = c(idx.keep, rownames(sd_filter_featwise(mat_tmp, opt$sd_filt_min))) # append features to idx.keep
    }
    mat_s = mat_s[idx.keep, ] # subset expression matrix to idx.keep-features
  } else if (opt$sd_filt_mode == 'equal') { # filter each ome individually
    #### apply global filter ####
    mat_s = sd_filter_featwise(mat_s, opt$sd_filt_min)
    
    #### determine size of smallest ome ####
    size_smallest_ome = filter(comb_rdesc, id %in% rownames(mat_s))$ome %>% # filter rdesc to mat_s
      table() %>% # sum remaining features
      min(na.rm=TRUE) # get size of smallest ome (post filtering)
    
    #### apply separate filters, to ensure equal-sized omes ####
    idx.keep = c() # initialize vector for features to keep
    for (ome in unique(comb_rdesc$ome_type)) {
      feat_byOme = filter(comb_rdesc, ome_type==ome)$id # get features-ids that correspond to this ome
      mat_tmp = mat_s[intersect(feat_byOme,rownames(mat_s)),] # subset expression-matrix to these features
      # get percentile-equivalent for size_smallest_ome
      sd_filt_equal = 1-size_smallest_ome/dim(mat_tmp)[1] # percentile to filter (1 - number of features / current size)
      # re-filter
      idx.keep = c(idx.keep, rownames(sd_filter_featwise(mat_tmp, sd_filt_equal))) # append features to idx.keep
    }
    mat_s = mat_s[idx.keep, ] # subset expression matrix to idx.keep-features
    
  } else {
    #### print error if incorrect filter selected ####
    stop(paste0("SD-filter mode '", opt$sd_filt_mode,"' not found. Please select either 'global', 'separate', or 'equal'"))
  }
  
  #### print-out post-filtering results ####
  # print results & overwrite comb_mat
  cat(glue("\nRemoved {dim(comb_mat)[1]-dim(mat_s)[1]} (of {dim(comb_mat)[1]}) features with standard-deviation below threshold ({opt$sd_filt_mode} filtering), leaving {dim(mat_s)[1]} fully-quantified features.\n\n"))
  ## filter data
  filter.df = data.frame(Status = factor( c(rep("Unfiltered", dim(comb_mat)[1]), # unfiltered vs filtered label
                                            rep("Filtered", dim(mat_s)[1])),
                                          levels = c("Unfiltered", "Filtered")), # factor so its ordered
                         id = c( rownames(comb_mat), # rownames of unfilt / filt features
                                 rownames(mat_s)),
                         sd = c( apply(comb_mat, 1, sd),# sd() of unfilt / filt features
                                 apply(mat_s, 1, sd))) %>%
    mutate( ome_type = comb_rdesc[match(id, comb_rdesc$id), "ome_type"]) # pull ome-type from rdesc
  ## print to boxplot
  pdf(paste0(opt$output_prefix,'_boxplots_SDfilter_',opt$sd_filt_mode,'.pdf')) # open PDF
  p = ggplot(filter.df, aes(x=ome_type, y=sd, fill=Status)) +
    geom_boxplot() + # boxplot 
    stat_summary(fun.data = function(x){ c(y = ceiling(max(filter.df$sd)), label = length(x)) }, # n_features labels
                 geom = "label", alpha=0.5, hjust = 0.5, position = position_dodge(1)) +
    xlab("-ome type") + ylab("Standard Deviation") +
    scale_fill_manual(values=c("#4477AA", "#66CCEE")) + theme_minimal() # color formatting
  plot(p)
  dev.off() # close PDF
}
comb_mat = mat_s # overwrite with (optionally) sd-filtered matrix

#### z-score multi-omic data matrix ####
mat_z = comb_mat # initialize
if(!is_null(opt$z_score_mode)){
  cat(glue("\n####################\nApplying Z-Scoring {opt$z_score_mode}\n"))
  if(opt$z_score_mode %in% c("row", "rowcol")) # if we are row-normalizing
    mat_z <- apply(mat_z, 1, function(x) (x-mean(x))/sd(x)) %>% # z-score rows
      t() # untransform matrix
  if(opt$z_score_mode %in% c("col", "rowcol")) # if we are column normalizing
    mat_z <- apply(mat_z, 2, function(x) (x-mean(x))/sd(x)) # z-score cols (do samples AFTER feature-wise z-score)
  
  # if an invalid method of zscoring was provided, complain
  if(! opt$z_score_mode %in% c("row", "col", "rowcol")) {
    #### print error if incorrect filter selected ####
    stop(paste0("Z-score mode '", opt$z_score_mode,"' not found. Please select either 'row', 'col', or 'rowcol'"))
  }
  
  # PDF with outliers
  pdf(paste0(opt$output_prefix,'_boxplots_zscore_',opt$z_score_mode,'.pdf'))
  boxplot(comb_mat, cex=0.5, main='before z-scoring') # create boxplot object
  boxplot(mat_z, cex=0.5, main=glue('after z-scoring ({opt$z_score_mode})'))#, outline=F)
  dev.off()
}
comb_mat = mat_z # overwrite with (optionally) z-scored matrix

###################################################
##               Create Combined GCT Object
###################################################
gct_comb = GCT(mat = comb_mat, # fully quantified, optionally z-scored matrix
               rdesc = comb_rdesc[match(rownames(comb_mat), comb_rdesc$id), ], # subset rdesc to match final comb_mat
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
  opts <- paste('vp', detectCores(),'t', sep='') # pass n.cores to NMF
  for(rank in ranks) # for each rank
    res.rank[[rank]] <-  nmf(gct_NN_fin@mat, rank=rank, method=opt$nmf_method, seed=opt$seed, nrun=opt$nrun, .options = opts)
}

#### Calculate Cluster Metrics ####
if (length(ranks) > 1) { # if we have more than one rank
  #### Calculate Cluster Metrics ####
  ## silhouette
  rank.sil <- lapply(res.rank[ranks], silhouette)
  rank.sil.avg <- lapply(rank.sil, function(x) tapply( x[,3], x[, 1], mean))
  
  ## cophenetic correlation
  # measure of clustering stability -- range [0,1] perfect  = 1
  rank.coph <- sapply(res.rank[ranks], cophcor)
  
  ## dispersion of consensus matrix
  # measure of clustering reproducibility -- range [0,1] perfect = 1
  rank.disp <- sapply(res.rank[ranks], dispersion)
  
  ## combine coph & disp
  # range [0,1] perfect = 1
  rank.coph.disp <- rank.disp^(1-rank.coph)
  
  
  #### Determine Best Rank ####
  # helper function to determine index of best score
  get.best.score.index <- function(scores, # vector of scores, where higher score = better score
                                   inc_thresh = 1e-6 # percent-increase threshold before a new score is considered "better" than the previous top score
  ) {
    top.i = 1 # initialize best score-index at 1
    for (i in 2:length(scores)) { # for every other score
      rel_inc = (scores[i] - scores[top.i]) / abs(scores[top.i]) # calculate relative increase (or possibly decrease)
      if ( rel_inc > inc_thresh ) # if the new score surpasses old score by a certain % increase
        top.i = i # record this index as the new top-score
    }
    return(top.i) # return the best score-index
  }
  # calculate top.rank (i.e. best cluster assignment)
  top.rank = get.best.score.index(scores = rank.coph.disp) %>% # determine index of best cluster-score
    ranks[.] # pull out cluster-rank corresponding to the best cluster-score
  

} else {top.rank = ranks}


#### Save NMF Object / Top Rank ####
save(file = "nmf_res.Rdata", res.rank)
write.table(top.rank, file = "nmf_best_rank.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)



###########################################################
##         Tar PDF Files
###########################################################

output.files = list.files(pattern=paste0(opt$output_prefix,".+.pdf"), full.names = T)
tar('NMF_preprocessing_figures.tar.gz',
    files = output.files, compression = "gzip", tar="tar")


