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
                       args = c('--metabolome_gct',"opt/input/ODG-v2_2-metabolomics-with-NMF.gct",
                                #'--metabolome_gct',"opt/input/ODG-v2_2-metabolomics_log_norm.gct",
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
gct_prot = parse_gctx(opt$proteome_gct)
gct_rna = parse_gctx(opt$rna_gct)

gct_contr = parse_gctx(opt$contrast_gct)

enrichementAnalysis <- function(gct, cov_of_interest, value_of_interest, label,
                                id_type = "name", pathway = "smpdb", match_no_hits = FALSE,
                                rowNorm="NULL", # Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, "CompNorm" for Normalization by a reference feature, "SumNorm" for Normalization to constant sum, "MedianNorm" for Normalization to sample median, and "SpecNorm" for Normalization by a sample-specific factor.
                                transNorm="NULL", # Select option to transform the data, "LogNorm" for Log Normalization, and "CrNorm" for Cubic Root Transformation.
                                scaleNorm="NULL", # Select option for scaling the data, "MeanCenter" for Mean Centering, "AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.
                                metabolite_id_col="Metabolite", sample_id_col="Sample.ID") {
  # write .txt file for Read.TextData() to parse
  df = t(gct@mat) %>% # transform matrix to samples as rows
    as.tibble(df) %>% # convert to tibble to prevent special-character replacement
    set_names(gct@rid) %>% # set colnames as original Metabolite IDs
    # set_names(gct@rdesc[[metabolite_id_col]]) %>% # set colnames as original Metabolite IDs
    mutate(Sample.ID = gct@cdesc[[sample_id_col]]) %>% # add Sample.ID as column
    mutate(cls = gct@cdesc[[cov_of_interest]]) %>%
    mutate(cls = ifelse(cls==value_of_interest, cls, glue("not_{value_of_interest}"))) %>% # add Sample.ID as column
    select(Sample.ID, cls, everything())
  write_csv(df, file = "tmp.csv")
  
  rm(list = c('mSet')) # clear previous mSet object
  # Initialize mSetObj & add data from Contrasts GCT
  mSet <- InitDataObjects("conc", "msetqea", FALSE) 
  mSet <- Read.TextData(mSet, "tmp.csv", "rowu", "disc")
  
  # Perform cross-referencing of compound names
  mSet <- CrossReferencing(mSet, id_type, lipid=F);
  
  # match_no_hits=FALSE # toggle for metabolites without hits should be soft-searched against MetaboAnalyst's databases
  no_hits = mSet$name.map$query.vec[!mSet$name.map$match.state] # get list of metabolites with no hits
  if (length(no_hits)>0 && match_no_hits) { for (id in no_hits) {
    cat(glue("\n\n#####\nMetabolite '{id}' has the following potential matches:\n\n"))
    print(PerformDetailMatch(mSet, id)$dataSet$candidates)
  }} else {
    cat(glue("\n\n#####\n{length(no_hits)} of {length(mSet$name.map$query.vec)} metabolites did not map to known compound names\n\n"))
  }
  
  mSet<-CreateMappingResultTable(mSet) # creates the mSet$dataSet$map.table entry (map between ID types)
  
  
  #### Preprocessing Steps ####
  # Mandatory check of data 
  mSet1<-SanityCheckData(mSet);
  if (is.numeric(mSet1)) { stop("Error code in SanityCheckData().") } else { mSet <- mSet1 }
  # Replace missing values with minimum concentration levels (required, even when no missing exist)
  # ToDo: NOTE that this replaces things with half the smallest positive value-- it ASSUMES we have all positive, which doesn't work for z-scores 
  mSet<-ReplaceMin(mSet);
  # Perform no normalization
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, rowNorm=rowNorm, transNorm=transNorm, scaleNorm=scaleNorm, ratio=FALSE, ratioNum=20)
  
  
  #### Normalization Plots ####
  # Plot normalization
  mSet<-PlotNormSummary(mSet, glue("{label}_norm_0_"), "png", 72, width=NA)
  # Plot sample-wise normalization
  mSet<-PlotSampleNormSummary(mSet, glue("{label}_snorm_0_"), "png", 72, width=NA)
  
  
  #### Finalize Enrichement Setup & Visualize ####
  # Set the metabolome filter
  mSet<-SetMetabolomeFilter(mSet, F); # turn metabolite filter off
  # Set the metabolite set library to pathway
  mSet<-SetCurrentMsetLib(mSet, glue("{pathway}_pathway"), 0); # might be able to change pathways here
  # Calculate the global test score
  mSet<-CalculateGlobalTestScore(mSet)
  # Plot the QEA
  mSet<-PlotQEA.Overview(mSet, glue("{label}_qea_0_"), "bar", "png", 72, width=NA)
  
  
  #### Pull Out Enrichement Analysis Data ####
  qea.mat = mSet$analSet$qea.mat # Q and P values
  qea.pvals = mSet$analSet$qea.pvals # p-values of individual metabolites
  
  write.csv(qea.mat, glue("{label}_qea_results.csv"))
}


# variables for enrichement analysis
gct = gct_meta
cov_of_interest = "Type"
id_type = "name"
metabolite_id_col = "Metabolite"
sample_id_col = "Sample.ID"
pathway = "smpdb"
# pathway = "kegg"

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
  enrichementAnalysis(gct, cov_of_interest, value_of_interest, id_type = id_type,
                      metabolite_id_col = metabolite_id_col, sample_id_col = sample_id_col,
                      transNorm = "LogNorm",
                      label = glue("{opt$output_prefix}_{value_of_interest}"), pathway=pathway)
}

# # file-copy for docker testing
# file.copy(list.files(pattern=opt$output_prefix), 'opt/input/', overwrite=T)

# # Set organism to human, at the moment only human data can be accomodated
# mSet<-SetOrganism(mSet, "hsa")



