#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))

# specify command line arguments
option_list <- list(
  #### Input Parameters ####
  make_option( c("-p", "--phosphoproteome_gct"), action='store', type='character',  dest='gct_pSTY', help='Path to input phosphoproteome GCT file.'),
  make_option( c("-a", "--acetylome_gct"), action='store', type='character',  dest='gct_acK', help='Path to input acetylome GCT file.'),
  make_option( c("-u", "--ubiquitylome_gct"), action='store', type='character',  dest='gct_ubK', help='Path to input ubiquitylome GCT file.'),
  #### Analysis Parameters ####
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest. If not provided, all annotations in the cdesc will be analyzed.'),
  make_option( c("-s", "--sample_id_col"), action='store', type='character',  dest='sample_id_col', help='Sample ID column in groups-file.'),
  make_option( c("-f", "--fdr_cutoff"), action='store', type='numeric', dest='fdr_cutoff', help='FDR cutoff for significant genes.'),
  make_option( c("-m", "--min_samples"), action='store', type='numeric', dest='min_samples', help='Minimum number of samples an annotation-subgroup can have before it is dropped from analysis (5 recommended, 3 minimum).'),
  make_option( c("-l", "--max_annot_levels"), action='store', type='numeric', dest='max_annot_levels', help='Maximum number of levels an annotation can have and be considered discrete.'), # default='10'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action='store', type='character',  dest='libdir', help='Folder to source from.'),
  make_option( c("-y", "--yaml_file"), action='store', type='character',  dest='yaml_file', help='yaml parameter file.', default = 'NA')
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # # for testing arguments
                       # args = c('--phosphoproteome_gct',"opt/input/ODG-v3-phosphoproteome-SpectrumMill-ratio-QCfilter-NArm.gct",
                       #          '--acetylome_gct',"opt/input/ODG-v3-acetylome-SpectrumMill-ratio-QCfilter-NArm.gct",
                       #          '--ubiquitylome_gct',"opt/input/ODG-v3-ubiquitylome-SpectrumMill-ratio-QCfilter-NArm.gct",
                       #          '-g',"opt/input/groups-subset.csv",
                       #          '-f',"0.05",
                       #          '-y',"opt/input/master-parameters.yaml",
                       #          '-x',"ODG_v3")
)


#### Parse YAML Arguments ####
opt = opt_cmd # initialize options with command line options
if ( !is.null(opt$yaml_file) ) {
  #### read in yaml ####
  library(yaml)
  yaml_out <- read_yaml(opt$yaml_file)
  #### overwrite the command-line parameters ####
  # # global parameters
  if (is.null(opt$sample_id_col)) opt$sample_id_col = yaml_out$DEV_sample_annotation$sample_id_col_name
  # diff-exp parameters
  if (is.null(opt$fdr_cutoff)) opt$fdr_cutoff = yaml_out$panoply_clumps_ptm$diff_exp$fdr_cutoff # toDo: consider changing to a clumps-ptm parameter
  if (is.null(opt$min_samples)) opt$min_samples = yaml_out$panoply_clumps_ptm$diff_exp$min_samples # toDo: consider changing to a clumps-ptm parameter
  if (is.null(opt$max_annot_levels)) opt$max_annot_levels = yaml_out$panoply_clumps_ptm$diff_exp$max_annot_levels # toDo: consider changing to a clumps-ptm parameter
} else { # if no YAML was provieded
  # check if any necessary parameters are missing
  if( any(sapply(list(), 
                 is.null)) ) { # if we have at least one missing parameter
    stop("Master Parameter yaml-file is missing. Please either provide a master-parameters file, or manually provide all parameters.") # error and stop
  }
}

# print parameters
cat("\n####################\nPARAMETERS:\n\n")
print(opt)
cat("####################\n")



library(tidyverse)
library(glue)
library(cmapR)
library(limma)
library(SimDesign) # for quiet() function


###################################################
##      Data Import / Preprocessing
###################################################

#### Source Files ####
source('https://raw.githubusercontent.com/broadinstitute/protigy/master/src/modT.R') # for modT.test.2class()

#### Parameter Wrangling ####

# import GCT files
if ( is.null(opt$gct_pSTY) && is.null(opt$gct_ubK) && is.null(opt$gct_acK) ) stop("No GCT files provided. Please provide at least one PTM GCT file.")
ome_list = list()
if (!is.null(opt$gct_pSTY)) ome_list[['pSTY']] = list(type = 'phosphoproteome', gct = parse_gctx(opt$gct_pSTY))
if (!is.null(opt$gct_ubK)) ome_list[['acK']] = list(type = 'ubiquitylome', gct = parse_gctx(opt$gct_ubK))
if (!is.null(opt$gct_acK)) ome_list[['ubK']] = list(type = 'acetylome', gct = parse_gctx(opt$gct_acK))

#### Read in Annotations & ID Mapping ####
if (!is.null(opt$groups_file)) { annots = read.csv(opt$groups_file) } else { # read in groups file or
  annots = ome_list[[1]]$gct@cdesc; cat(glue("\nWARNING: No groups file provided; module will attempt to use {ome_list[[1]]$type} gct@cdesc.\n")) } # print warning and use cdescs



###################################################
##      Differential Expression
###################################################


# helper function to determine if annot is blank / should be skipped
skip.annot = function(value_of_interest) {
  return( is.na(value_of_interest) ||
            value_of_interest=="" || # match blank
            grepl("^[nN]\\.*[aA]\\.*$",value_of_interest) || # match NA / N.A. of any case
            grepl("^[nN][aA][nN]$",value_of_interest)) # match NaN
}

# initialize log object to keep track of which annotations have what files
log_file = data.frame(annot.name = character(0),
                      subvalue.name = character(0), # valid unique subvalue name (for file lookup)
                      subvalue.value = character(0), # original subvaluename
                      valid.subvalue = logical(0),
                      valid.results = logical(0))


for (annot_of_interest in names(annots)) {
  # skip annotation if it has too many values, or too few
  if (length(unique(annots[[annot_of_interest]])) > opt$max_annot_levels) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; too many annotation-values ({length(unique(annots[[annot_of_interest]]))}) to be considered discrete (>{opt$max_annot_levels}).\n\n")); next }
  if (length(unique(annots[[annot_of_interest]]))  == 1 ) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; only one unique annotation-value.\n\n")); next }
  
  # filter out samples with missing / empty / NA values
  cid.keep = annots[['Sample.ID']][!sapply(annots[[annot_of_interest]], skip.annot)]
  if (length(cid.keep) < opt$min_samples*2) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation; only {length(cid.keep)} non-missing annotation values, which is not enough for any valid comaprisons.\n\n")); next }
  
  #### Set Up Log-File / Directory ####
  # initialize log_file with all FALSE
  log_file.tmp = data.frame(annot.name = annot_of_interest,
                            subvalue.name = make.names(sort(unique(annots[[annot_of_interest]]))), # valid unique subvalue name (for file lookup)
                            subvalue.value = sort(unique(annots[[annot_of_interest]])), # original subvaluename
                            # use subvalue.name as row.names for quick lookup / edits
                            row.names = make.names(sort(unique(annots[[annot_of_interest]]))), # valid unique subvalue name (for file lookup)
                            valid.subvalue = FALSE,
                            valid.results = FALSE)
  # initialize dataframe for 
  df_full = data.frame(index = character(0),
                       id = character(0),
                       id.name = character(0),
                       feature = character(0),
                       logFC = numeric(0),
                       P.Value = numeric(0),
                       adj.P.Val = numeric(0))
  
  # otherwise run analysis for every unique annotation subvalue
  cat(glue("\n\n####################\nAnalyzing '{annot_of_interest}' Annotation \n####################\n\n"))
  for (ome in ome_list) {
    # wrangle data & annotations
    gct.tmp = subset_gct(ome$gct, cid=cid.keep) # filter gct to match annots file / drop blank samples
    if (length(gct.tmp@cid) < opt$min_samples*2) { cat(glue("\n\n####################\nSkipping '{annot_of_interest}' annotation for {ome$type}; only {length(gct.tmp@cid)} non-missing annotation values, which is not enough for any valid comaprisons.\n\n")); next }
    annots.tmp = annots[match(gct.tmp@cid, annots[[opt$sample_id_col]]), ] # filter annots to match gct file & sort in same order
    # sanity check
    if (sum(annots.tmp[[opt$sample_id_col]]!=gct.tmp@cid)!=0) stop(glue("GCT samples do not align with annotation file, for the  '{annot_of_interest}' annotation in the {ome$type}. Something has gone terribly wrong!"))
    
    # analyze each annotation subvalue
    for (value_of_interest.name in rownames(log_file.tmp)) {
      value_of_interest = log_file.tmp[value_of_interest.name,"subvalue.value"] # get original annot value for use in analysis
      if (skip.annot(value_of_interest)) { cat(glue("\nSkipping '{value_of_interest}' annotation.\n\n")); next }
      log_file.tmp[value_of_interest.name,"valid.subvalue"] = TRUE # mark subvalue as a valid subvalue
      
      # create binary CLS
      cls = annots.tmp[[annot_of_interest]] %>%
        { ifelse(!is.na(.) & .==value_of_interest, ., glue("{value_of_interest}_not")) } %>% # NA -> no_<value_of_interest>
        factor(levels = c(glue("{value_of_interest}_not"), value_of_interest)) # factor cls with _not first, to force order of comparison
      if ( sum(cls==value_of_interest) < opt$min_samples ) { cat(glue("\nSkipping '{value_of_interest}' annotation for {ome$type}; too few samples ({sum(cls==value_of_interest)}) to perform a valid comparison.\n\n")); next }
      
      #### T-Test ####
      d = rownames_to_column(as.data.frame(gct.tmp@mat), "feature_id")
      out = quiet(modT.test.2class(d, 'tmp', groups=cls, id.col = "feature_id")) # wrapped in quiet() to suppress repeated printouts
      out_df = out$output %>%
        dplyr::arrange(P.Value) %>%
        # { if (!is.null(n_feat)) {top_n(., n=n_feat, wt = -P.Value)} else {.} }
        { if (!is.null(opt$fdr_cutoff)) { dplyr::filter(., adj.P.Val < opt$fdr_cutoff) } else {.} } # filter to significant logFC
      df = data.frame(index = rownames(out_df['logFC']),
                      logFC = round(out_df[['logFC']], 5),
                      P.Value = out_df[['P.Value']],
                      adj.P.Val = out_df[['adj.P.Val']]) %>%
        mutate(id = value_of_interest,
               id.name = value_of_interest.name, # make.names(value_of_interest)
               feature = ome$type)
      
      if ( dim(df)[1]>0 ) {
        df_full = rbind(df_full, df[names(df_full)]) # append to full dataframe
        log_file.tmp[value_of_interest.name,"valid.results"] = TRUE # mark subvalue as having valid results
      } else {
        cat(glue("\nAnnotation subvalue '{value_of_interest}' had no significant features in {ome$type}.\n\n"))
      }
    }
  }
  
  
  # write log file
  if (dim(df_full)[1]>0) {
    write.csv(df_full, sep='\t',
              file = glue("{opt$output_prefix}_{annot_of_interest}_diff_exp.tsv"))
  } else {
    cat(glue("\nNo significant features found for '{annot_of_interest}' across any ome.\n\n"))
  }

  
  
  #### Add to Log File (even if there was no enrichment) ####
  rownames(log_file.tmp) = NULL # drop rownames
  log_file = rbind(log_file, log_file.tmp) # append tmp logfile to logfile
}

# write log file
write.csv(log_file,
          file = glue("{opt$output_prefix}_log_file.csv"))



