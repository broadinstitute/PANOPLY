#!/usr/bin/env Rscript
#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### Input Parameters ####
  make_option( c("-m", "--metabolome_gct"), action='store', type='character',  dest='metabolome_gct', help='GCT file containing expression data for the metabolome.'),
  make_option( c("-n", "--metabolome_norm_gct"), action='store', type='character',  dest='metabolome_norm_gct', help='GCT file containing expression data for the NORMALIZED metabolome.'),
  make_option( c("-c", "--contrast_gct"), action='store', type='character',  dest='contrast_gct', help='GCT file containing contrast values for the metabolome, across some annotation.'),
  make_option( c("-o", "--proteome_gct"), action='store', type='character',  dest='proteome_gct', help='GCT file containing expression data for the proteome'),
  make_option( c("-r", "--rna_gct"), action='store', type='character',  dest='rna_gct', help='GCT file containing expression data for the transcriptome'),
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest. If not provided, all annotations in the cdesc will be analyzed.'),
  #### Analysis ####
  make_option( c("-a", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  make_option( c("-p", "--pval_signif"), action='store', type='numeric', dest='pval_signif', help='P-value threshold for significant enrichement.'), 
  make_option( c("-q", "--feature_fdr"), action='store', type='numeric', dest='feature_fdr', help='Max FDR threshold for feature-selection 2-sample T-test.'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-y", "--yaml"), action="store", dest='yaml_file', type="character", help="Path to .yaml file with parameters."),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # for testing arguments
                       args = c('--metabolome_gct',"opt/input/ODG-v2_2-metabolomics-with-NMF-HMDB_UNIQUE.gct",
                                '--metabolome_norm_gct',"opt/input/ODG-v2_2-metabolomics_log_norm-HMDB_UNIQUE.gct",
                                '--contrast_gct',"opt/input/contrasts/metabolomics_logNorm-NMF.consensus-contrast.gct",
                                '--proteome_gct',"opt/input/ODG-v2_2-proteome-SpectrumMill-ratio-QCfilter-NArm-with-NMF.gct",
                                '--rna_gct',"opt/input/ODG-v2_2-rnaseq-expression-TPM-protein-coding-log2-median-norm-NArm-with-NMF.gct",
                                # '-g',"opt/input/groups-subset.csv",
                                '-y',"opt/input/master-parameters.yaml",
                                '-x',"ODG_v2_2")
)
opt = opt_cmd # ToDo: Add YAML parameters (temporarily setting opt straight from command-line params)



library(MetaboAnalystR) # analysis tools for Metabolomic Data
# library(RefMet) # mapping between RefMet and other ID types
library(glue)
library(cmapR)
library(tidyverse)

gct_meta = parse_gctx(opt$metabolome_gct)
gct_meta_norm = parse_gctx(opt$metabolome_norm_gct)
gct_prot = parse_gctx(opt$proteome_gct)
gct_rna = parse_gctx(opt$rna_gct)

gct_contr = parse_gctx(opt$contrast_gct)


mSetPrep <- function(gct, anal.type,
                     metabolite_id_col, id_type, sample_id_col,
                     cov_of_interest, value_of_interest=NULL, 
                     rowNorm="NULL", # Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, "CompNorm" for Normalization by a reference feature, "SumNorm" for Normalization to constant sum, "MedianNorm" for Normalization to sample median, and "SpecNorm" for Normalization by a sample-specific factor.
                     transNorm="NULL", # Select option to transform the data, "LogNorm" for Log Normalization, and "CrNorm" for Cubic Root Transformation.
                     scaleNorm="NULL", # Select option for scaling the data, "MeanCenter" for Mean Centering, "AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.)
                     match_no_hits = F, label = NULL) {
  if (is.null(label)) {label = anal.type}
  # write .txt file for Read.TextData() to parse
  df = t(gct@mat) %>% # transform matrix to samples as rows
    as.tibble(df) %>% # convert to tibble to prevent special-character replacement
    # set_names(gct@rid) %>% # set colnames as original Metabolite IDs
    set_names(gct@rdesc[[metabolite_id_col]]) %>% # set colnames as original Metabolite IDs
    mutate(Sample.ID = gct@cdesc[[sample_id_col]]) %>% # add Sample.ID as column
    mutate(cls = gct@cdesc[[cov_of_interest]]) %>%
    {if (!is.null(value_of_interest)) {mutate(., cls = ifelse(cls==value_of_interest, cls,
                                                              glue("not_{value_of_interest}")))} else {.} } %>%
    select(Sample.ID, cls, everything())
  input_fn = tempfile(pattern = anal.type, fileext = ".csv")
  write_csv(df, file = input_fn)
  
  # Initialize mSetObj & add data from Contrasts GCT
  mSet=InitDataObjects("conc", anal.type, FALSE)
  mSet=Read.TextData(mSet, input_fn, "rowu", "disc")
  
  # Perform cross-referencing of compound names
  mSet=CrossReferencing(mSet, id_type)
  
  # match_no_hits=FALSE # toggle for metabolites without hits should be soft-searched against MetaboAnalyst's databases
  no_hits = mSet$name.map$query.vec[!mSet$name.map$match.state] # get list of metabolites with no hits
  if (length(no_hits)>0 && match_no_hits) { for (id in no_hits) {
    cat(glue("\n\n#####\nMetabolite '{id}' has the following potential matches:\n\n"))
    print(PerformDetailMatch(mSet, id)$dataSet$candidates)
  }} else {
    cat(glue("\n\n#####\n{length(no_hits)} of {length(mSet$name.map$query.vec)} metabolites did not map to known compound names\n\n"))
  }
  
  if (anal.type=="msetqea") { # if we're doing mSetQEA
    tmp = anal.type; assign('anal.type', 'tmp', envir = .GlobalEnv)  # temporarily change anal.type to get around hardcoded mSet global-variable reference to in CreateMappingResultTable()
    anal.type_resetToggle=T # and set a toggle to ensure it gets fixed after 
  } else {anal.type_resetToggle=F}
  mSet<-CreateMappingResultTable(mSet)
  if (anal.type_resetToggle) { assign('anal.type', tmp, envir = .GlobalEnv) } # reset to original analysis type
  
  
  #### Preprocessing Steps ####
  # Mandatory check of data 
  mSet1=SanityCheckData(mSet)
  if (is.numeric(mSet1)) { stop("Error code in SanityCheckData().") } else { mSet=mSet1 }
  # Replace missing values with minimum concentration levels (required, even when no missing exist)
  # ToDo: NOTE that this replaces things with half the smallest positive value-- it ASSUMES we have all positive, which doesn't work for z-scores 
  mSet=ReplaceMin(mSet)
  # Perform no normalization
  mSet=PreparePrenormData(mSet)
  mSet=Normalization(mSet, rowNorm=rowNorm, transNorm=transNorm, scaleNorm=scaleNorm, ratio=FALSE, ratioNum=20)
  
  
  #### Normalization Plots ####
  # Plot normalization
  mSet=PlotNormSummary(mSet, glue("{label}_norm_0_"), "png", 72, width=NA)
  # Plot sample-wise normalization
  mSet=PlotSampleNormSummary(mSet, glue("{label}_snorm_0_"), "png", 72, width=NA)
  
  return(mSet)
}

mSetCreateReport <- function(mSetObj) {
  if (exists('mSet', envir = .GlobalEnv)) {
    message('\nWARNING: Global mSet variable already exist. Object will be temporarily saved and reinstated on.exit().\nExisting MetaboAnalystR analyses will likely be disrupted.\n')
    tmp_global_mSet = mSet # store mSet object
    rm(envir = .GlobalEnv, list=c('mSet')) # remove global mSet object to avoid conflicts with new analysis
    on.exit(expr = assign('mSet', tmp_global_mSet, envir = .GlobalEnv), add = T) # and reassign it on-exit
  }
  
  on.exit(expr = {
    cat("Removing intermediary global variables set during MetaboAnalyst Report-Generation.\n########################\n")
    rm(list = c('fig.count', 'rnwFile', 'table.count'), envir = .GlobalEnv)
  }, add=T)
  
  assign('mSet', mSetObj, envir = .GlobalEnv) # briefly set mSet as a global variable
  PreparePDFReport(mSet, 'PANOPLY')
  rm(envir = .GlobalEnv, list=c('mSet')) # remove that mSet object to avoid problems later
}


enrichementAnalysis <- function(gct, metabolite_id_col, sample_id_col, id_type, 
                                cov_of_interest, value_of_interest=NULL, 
                                pathway = "smpdb_pathway", label = "results",
                                rowNorm="NULL", # Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, "CompNorm" for Normalization by a reference feature, "SumNorm" for Normalization to constant sum, "MedianNorm" for Normalization to sample median, and "SpecNorm" for Normalization by a sample-specific factor.
                                transNorm="NULL", # Select option to transform the data, "LogNorm" for Log Normalization, and "CrNorm" for Cubic Root Transformation.
                                scaleNorm="NULL", # Select option for scaling the data, "MeanCenter" for Mean Centering, "AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.
                                match_no_hits = FALSE, create_report = TRUE) {
  #### On-Exit Behavior (MetaboAnalyst Cleanup) ####
  on.exit(expr = {
    # rm(list = c('mSetQEA')) # remove mSetQEA object; these functions seem to only re
    cat("\n########################\nRemoving all intermediary .qs files.\n")
    file.remove(list.files(pattern="\\.qs")) # remove all .qs files
    # cat("\n\n#####\nRemoving intermediary mSetQEA object.\n\n")
    # rm(list = c('mSetQEA'), envir = .GlobalEnv)
    cat("Removing intermediary global variables set by MetaboAnalyst.\n########################\n")
    rm(list = c('anal.type', 'api.base', 'current.msetlib', 'mdata.all',
                'mdata.siggenes', 'meta.selected', 'metaboanalyst_env',
                'module.count', 'plink.path', 'primary.user', 
                'smpdbpw.count', 'url.pre'), envir = .GlobalEnv)
  }, add=T)
  
  mSetQEA = mSetPrep(gct, anal.type = 'msetqea', label=label,
                     metabolite_id_col=metabolite_id_col, id_type=id_type,
                     sample_id_col=sample_id_col,
                     cov_of_interest=cov_of_interest, value_of_interest=value_of_interest,
                     rowNorm=rowNorm,transNorm=transNorm,scaleNorm=scaleNorm, match_no_hits=match_no_hits)
  
  #### Finalize Enrichement Setup & Visualize ####
  # Set the metabolome filter
  mSetQEA=SetMetabolomeFilter(mSetQEA, F) # turn metabolite filter off
  # Set the metabolite set library to pathway
  mSetQEA=SetCurrentMsetLib(mSetQEA, pathway, 0) # might be able to change pathways here
  # Calculate the global test score
  mSetQEA=CalculateGlobalTestScore(mSetQEA)
  
  # Plot the QEA Barchart
  mSetQEA=PlotQEA.Overview(mSetQEA, glue("{label}_qea_{pathway}_barPlot_"), "bar", "pdf", width=NA)
  mSetQEA=PlotEnrichDotPlot(mSetQEA, imgName = glue("{label}_qea_{pathway}_dotPlot_"), enrichType = "qea", "pdf", width=NA)
  # mSetQEA=with(environment(MetaboAnalystR::PlotEnrichNet.Overview), PlotEnrichPieChart(mSetQEA, imgName = glue("{label}_qea_{pathway}_pieChart_"), enrichType = "qea", "png", 72, width=NA)) # i dunno exactly what this shows in this context
  
  # Plot the QEA Network Chart
  folds = mSetQEA$analSet$qea.mat[, 3]/mSetQEA$analSet$qea.mat[, 4]
  GetShortNames = with(environment(MetaboAnalystR::PlotEnrichNet.Overview), 
                       GetShortNames) # steal the GetShortNames function from the MetaboAnalystR environment
  names(folds) = GetShortNames(rownames(mSetQEA$analSet$qea.mat)) 
  pvals = mSetQEA$analSet$qea.mat[, "Raw p"]
  mSetQEA$analSet$enrich.net = PlotEnrichNet.Overview(folds, pvals)
  file.copy("msea_network.json", glue("{label}_qea_{pathway}_network.json"), overwrite = T)
  
  
  #### Pull Out Enrichement Analysis Data ####
  qea.mat = mSetQEA$analSet$qea.mat # Q and P values
  qea.pvals = mSetQEA$analSet$qea.pvals # p-values of individual metabolites
  
  write.csv(qea.mat, glue("{label}_qea_{pathway}_results.csv"))
  
  if(create_report) {
    mSetCreateReport(mSetQEA)
    file.copy('Analysis_Report.pdf', glue("{label}_qea_{pathway}_report.pdf"), overwrite = T)
  }
  
  return(mSetQEA)
}


networkAnalysis <- function() {
  #### On-Exit Behavior (MetaboAnalyst Cleanup) ####
  on.exit(expr = {
    # rm(list = c('mSetQEA')) # remove mSetQEA object; these functions seem to only re
    cat("\n########################\nRemoving all intermediary .qs files.\n")
    file.remove(list.files(pattern="\\.qs")) # remove all .qs files
    # cat("\n\n#####\nRemoving intermediary mSetQEA object.\n\n")
    # rm(list = c('mSetQEA'), envir = .GlobalEnv)
    # cat("Removing intermediary global variables set by MetaboAnalyst.\n########################\n")
    # rm(list = c('anal.type', 'api.base', 'current.msetlib', 'mdata.all',
    #             'mdata.siggenes', 'meta.selected', 'metaboanalyst_env',
    #             'module.count', 'plink.path', 'primary.user', 
    #             'smpdbpw.count', 'url.pre'), envir = .GlobalEnv)
  }, add=T)
  
  if (exists('mSet', envir = .GlobalEnv)) {
    message('\nWARNING: Global mSet variable already exist. Object will be temporarily saved and reinstated on.exit().\nExisting MetaboAnalystR analyses will likely be disrupted.\n')
    tmp_global_mSet = mSet # store mSet object
    rm(envir = .GlobalEnv, list=c('mSet')) # remove global mSet object to avoid conflicts with new analysis
    on.exit(expr = assign('mSet', tmp_global_mSet, envir = .GlobalEnv), add = T) # and reassign it on-exit
  }
  
  
  
  mSetNet<-PrepareNetworkData(mSetNet);
  
  
  PerformDSPC<-with(environment(MetaboAnalystR::SetKEGG.PathLib),
                    PerformDSPC)
  mSetNet<-PerformDSPC(mSetNet)
  
  mSetNet<-CreateGraph(mSetNet)
  
  if(create_report) {
    assign('mSet', mSetNet, envir = .GlobalEnv) # briefly set mSet as a global variable
    PreparePDFReport(mSet, 'PANOPLY')
    rm(envir = .GlobalEnv, list=c('mSet')) # remove that mSet object to avoid problems later
    file.copy('Analysis_Report.pdf', glue("{label}_network_report.pdf"))
    # }
  }
}

pathwayAnalysis <-  function(gct, metabolite_id_col, sample_id_col, id_type, 
                             cov_of_interest, value_of_interest=NULL, 
                             label = "results",
                             rowNorm="NULL", # Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, "CompNorm" for Normalization by a reference feature, "SumNorm" for Normalization to constant sum, "MedianNorm" for Normalization to sample median, and "SpecNorm" for Normalization by a sample-specific factor.
                             transNorm="NULL", # Select option to transform the data, "LogNorm" for Log Normalization, and "CrNorm" for Cubic Root Transformation.
                             scaleNorm="NULL", # Select option for scaling the data, "MeanCenter" for Mean Centering, "AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.
                             match_no_hits = FALSE, create_report = TRUE) {
  #### On-Exit Behavior (MetaboAnalyst Cleanup) ####
  on.exit(expr = {
    # rm(list = c('mSetQEA')) # remove mSetQEA object; these functions seem to only re
    cat("\n########################\nRemoving all intermediary .qs files.\n")
    file.remove(list.files(pattern="\\.qs")) # remove all .qs files
    # cat("\n\n#####\nRemoving intermediary mSetQEA object.\n\n")
    # rm(list = c('mSetQEA'), envir = .GlobalEnv)
    cat("Removing intermediary global variables set by MetaboAnalyst.\n########################\n")
    rm(list = c('anal.type', 'api.base', #'current.msetlib',
                'mdata.all', 'mdata.siggenes', 'meta.selected', 'metaboanalyst_env',
                'mem.univ', 'msg.vec',
                'module.count', 'plink.path', 'primary.user', 
                'smpdbpw.count', 'url.pre'), envir = .GlobalEnv)
  }, add=T)
  
  if (exists('mSet', envir = .GlobalEnv)) {
    message('\nWARNING: Global mSet variable already exist. Object will be temporarily saved and reinstated on.exit().\nExisting MetaboAnalystR analyses will likely be disrupted.\n')
    tmp_global_mSet = mSet # store mSet object
    rm(envir = .GlobalEnv, list=c('mSet')) # remove global mSet object to avoid conflicts with new analysis
    on.exit(expr = assign('mSet', tmp_global_mSet, envir = .GlobalEnv), add = T) # and reassign it on-exit
  }
  
  mSetPTW = mSetPrep(gct, anal.type = 'pathqea', label=label,
                     metabolite_id_col=metabolite_id_col, id_type=id_type,
                     sample_id_col=sample_id_col,
                     cov_of_interest=cov_of_interest, value_of_interest=value_of_interest,
                     rowNorm=rowNorm,transNorm=transNorm,scaleNorm=scaleNorm, match_no_hits=match_no_hits)
  
  mSetPTW<-SetKEGG.PathLib(mSetPTW, "hsa", "current")
  # mSetPTW<-SetSMPDB.PathLib(mSetPTW, "hsa")
  mSetPTW<-SetMetabolomeFilter(mSetPTW, F);
  
  mSetPTW1<-CalculateQeaScore(mSetPTW, "rbc", "gt") # API issues here... kinda at a loss.....
  if (is.numeric(mSetPTW1)) { stop("Error code in CalculateQeaScore().") } else { mSetPTW=mSetPTW1 }
  mSetPTW<-PlotPathSummary(mSetPTW, F, "path_view_0_", "png", 72, width=NA, NA, NA )
  
  
  ComputePathHeatmap<-with(environment(MetaboAnalystR::SetKEGG.PathLib),
                           ComputePathHeatmap)
  ComputePathHeatmap(mSetPTW, "hsa",glue("{label}_pathwayqea_heatmap.json"), "pathqea")
  
  
  mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA, NA, NA )
  
}


# variables for enrichement analysis
gct = gct_meta
cov_of_interest = "Type"
id_type = "hmdb"
metabolite_id_col = "HMDB_ID"
sample_id_col = "Sample.ID"
pathway = "smpdb_pathway"
# pathway = "kegg_pathway"

rowNorm="NULL" # Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, "CompNorm" for Normalization by a reference feature, "SumNorm" for Normalization to constant sum, "MedianNorm" for Normalization to sample median, and "SpecNorm" for Normalization by a sample-specific factor.
transNorm="logNorm" # Select option to transform the data, "LogNorm" for Log Normalization, and "CrNorm" for Cubic Root Transformation.
scaleNorm="NULL" # Select option for scaling the data, "MeanCenter" for Mean Centering, "AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.


# #### ID Mapping ####
# gct_full@rdesc$Standardized.name = refmet_map_df(gct_full@rdesc[[metabolite_id_col]])$Standardized.name # get standardized names with RefMet::refmet_map_df(), and add to rdesc
# refmet_map = read.csv(file.path(opt$lib_dir, 'refmet.csv'))
# gct_full@rdesc = left_join(gct_full@rdesc, refmet_map, by = c("Standardized.name" = "refmet_name"))
# # gct_full@rdesc$Metabolite[which(gct_full@rdesc$Standardized.name=="-")] # get IDs that didn't map
# # prune no-match or duplicates IDs
# # duplicated_pubchem_id = unique(gct_full@rdesc$pubchem_cid[which(duplicated(gct_full@rdesc$pubchem_cid))])
# # gct = subset_gct(gct_full, rid = which(!(gct_full@rdesc$pubchem_cid %in% duplicated_pubchem_id)))
# # metabolite_id_col = "pubchem_cid"
# # id_type = "pubchem"

#### Enrichement Analysis ####
for (value_of_interest in unique(gct@cdesc[[cov_of_interest]])) {
  out = enrichementAnalysis(gct, cov_of_interest=cov_of_interest, value_of_interest = value_of_interest,
                            metabolite_id_col = metabolite_id_col, sample_id_col = sample_id_col, id_type=id_type, 
                            rowNorm=rowNorm, transNorm = transNorm, scaleNorm=scaleNorm,
                            label = glue("{opt$output_prefix}_{cov_of_interest}_{value_of_interest}"),
                            pathway=pathway)
}

# file-copy for docker testing
file.copy(list.files(pattern=opt$output_prefix), 'opt/input/', overwrite=T)
file.copy(list.files(path='/tmp/', full.names=T, recursive = 1), '/opt/input/', overwrite=T)



#### Pathway Analysis ####
for (value_of_interest in unique(gct@cdesc[[cov_of_interest]])) {
  out = pathwayAnalysis(gct, cov_of_interest=cov_of_interest, value_of_interest = value_of_interest,
                        metabolite_id_col = metabolite_id_col, sample_id_col = sample_id_col, id_type=id_type, 
                        rowNorm=rowNorm, transNorm = transNorm, scaleNorm=scaleNorm,
                        label = glue("{opt$output_prefix}_{cov_of_interest}_{value_of_interest}"),
                        pathway=pathway)
}


#### LogFC ####

source('https://raw.githubusercontent.com/broadinstitute/protigy/master/src/modT.R')
# install.packages('statmod')
library(limma)
for (value_of_interest in unique(gct@cdesc[[cov_of_interest]])) {
  df_list = list(list(ome = "metabolome",
                      ft_id_col = metabolite_id_col,
                      gct = gct_meta_norm),
                 list(ome = "proteome",
                      ft_id_col = "geneSymbol",
                      gct = gct_prot),
                 list(ome = "RNA",
                      ft_id_col = "geneSymbol",
                      gct = gct_rna))
  
  for (data in df_list) {
    d = as.data.frame(data$gct@mat) %>% # transform matrix into tibble
      mutate(feature_id = data$gct@rdesc[[data$ft_id_col]]) %>% # add feature_id column
      select(feature_id, everything()) 
    d_subset = filter(d, !(feature_id %in% d[which(duplicated(d$feature_id)), 'feature_id']))
    cat(glue("{dim(d)[1]-dim(d_subset)[1]} features dropped (from original {dim(d)[1]} features) for duplicated '{data$ft_id_col}' IDs"))
    groups = data$gct@cdesc[[cov_of_interest]] %>%
      { ifelse(.==value_of_interest, ., glue("not_{value_of_interest}")) }
    
    out = modT.test.2class(d_subset, 'tmp', groups=groups, id.col = "feature_id")
    # n_feat = NULL # number of features
    out_df = out$output %>%
      arrange(P.Value) %>%
      # { if (!is.null(n_feat)) {top_n(., n=n_feat, wt = -P.Value)} else {.} }
      filter(adj.P.Val < 0.05) # filter to significant logFC
    df = data.frame(id = rownames(out_df['logFC']),
                    logFC = round(out_df[['logFC']], 5))
    write.table(df, file=glue("{opt$output_prefix}_logFC_signFeat_{data$ome}_{cov_of_interest}-{value_of_interest}.vs.{value_of_interest}_not.tsv"),
                sep = '\t', quote=FALSE, row.names = FALSE)
  }
}
file.copy(list.files(pattern=glue("{opt$output_prefix}_logFC_signFeat")), 'opt/input/', overwrite=T)

# note: can run internal functions with the following:
# with(environment(MetaboAnalystR::PlotEnrichDotPlot), GetMyHeatCols(92))

# # minimal code for creates current.msetlib object
# mSet=InitDataObjects("conc", "mSetQEA", FALSE)
# mSet = SetCurrentMsetLib(mSet, pathway, 0)
# current.msetlib[,3]

# # or, create it manually
# libname = pathway
# destfile <- paste(libname, ".qs", sep = "")
# my.qs <- paste("https://www.metaboanalyst.ca/resources/libs/msets/",
#                destfile, sep = "")
# if (!file.exists(destfile)) {
#   download.file(my.qs, destfile, method = "curl")
# }
# current.msetlib <- qs::qread(destfile)
# 
# ms.list <- iconv(current.msetlib[, 3], from = "utf8",
#                  to = "utf8")
# ms.list <- lapply(ms.list, function(x) unique(unlist(strsplit(x,
#                                                               "; ", fixed = TRUE))))
# names(ms.list) <- current.msetlib[, 2]
