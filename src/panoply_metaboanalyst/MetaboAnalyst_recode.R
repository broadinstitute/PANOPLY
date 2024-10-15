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
  make_option( c("-o", "--proteome_gct"), action='store', type='character',  dest='proteome_gct', help='GCT file containing expression data for the proteome'),
  make_option( c("-r", "--rna_gct"), action='store', type='character',  dest='rna_gct', help='GCT file containing expression data for the transcriptome'),
  make_option( c("-w", "--pathway_gmt"), action='store', type='character',  dest='pathway_gmt', help='GMT file containing pathways of interest.'),
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
                       args = c(#'--metabolome_gct',"opt/input/ODG-v2_2-metabolomics-with-NMF-HMDB_UNIQUE.gct",
                                '--metabolome_gct',"opt/input/HMDB_ID_GCTs/ODG-v2_2-metabolomics_log_norm-HMDB_UNIQUE.gct",
                                '--proteome_gct',"opt/input/ODG-v2_2-proteome-SpectrumMill-ratio-QCfilter-NArm-with-NMF.gct",
                                '--rna_gct',"opt/input/ODG-v2_2-rnaseq-expression-TPM-protein-coding-log2-median-norm-NArm-with-NMF.gct",
                                '--pathway_gmt',"opt/input/smpdb_pathway.gmt",
                                '-g',"opt/input/groups-subset.csv",
                                '-y',"opt/input/master-parameters.yaml",
                                '-x',"ODG_v2_2")
)
opt = opt_cmd # ToDo: Add YAML parameters (temporarily setting opt straight from command-line params)



library(glue)
library(cmapR)
library(tidyverse)
# BiocManager::install('globaltest')
library(globaltest)
# BiocManager::install('ActivePathways')
library(ActivePathways)


gct_meta = parse_gctx(opt$metabolome_gct)
# gct_meta_norm = parse_gctx(opt$metabolome_norm_gct)
gct_prot = parse_gctx(opt$proteome_gct)
gct_rna = parse_gctx(opt$rna_gct)

gct_contr = parse_gctx(opt$contrast_gct)

gmt = read.GMT(opt$pathway_gmt)


gct = gct_meta
# id_filter = read.csv('/opt/input/name_map.csv') %>%
#   dplyr::filter(as.logical(Comment)) %>%
#   .[['Query']]
# gct = subset_gct(gct_meta, rid = id_filter)

msea.data = as.matrix(t(gct@mat))
cov_of_interest = "Type"
value_of_interest="ODG"
cls_type = gct@cdesc[[cov_of_interest]] %>%
  { ifelse(!is.na(.) & .==value_of_interest, ., glue("{value_of_interest}_not")) } # NA -> no_<value_of_interest>
cls = gct@cdesc[[cov_of_interest]] %>%
  { ifelse(!is.na(.) & .==value_of_interest, 1, 0) } # NA -> no_<value_of_interest>

# data = data.frame(msea.data) %>%
#   mutate(cls = cls_type) %>%
#   rownames_to_column('Sample.ID') %>%
#   relocate(Sample.ID, cls)
# write.csv(data, file = "/opt/input/qea_input_ODG.csv", row.names = F)


#### QEA // Global Test ####

pathway_length = lapply(gmt, function(mset) {
  length(mset$genes)
})
hits = lapply(gmt, function(mset) {
  x = mset$genes
  x[x %in% colnames(msea.data)]
})

# if (sum(lapply(hits, length)==1) > min.overlap) {
gt.obj <- globaltest::gt(cls, msea.data, subsets = hits)
gt.res <- globaltest::result(gt.obj)


# rename matrix
raw.p <- gt.res[, 1]
bonf.p <- p.adjust(raw.p, "holm")
fdr.p <- p.adjust(raw.p, "fdr")
res.mat <- cbind(pathway_length, gt.res[, 5], gt.res[, 2],
                 gt.res[, 3], raw.p, bonf.p, fdr.p)
rownames(res.mat) <- rownames(gt.res)
colnames(res.mat) <- c("Total Cmpd", "Hits", "Statistic Q",
                       "Expected Q", "Raw p", "Holm p", "FDR")
# reorder based on p-value
hit.inx <- res.mat[, 2] > 0
res.mat <- res.mat[hit.inx, , drop = FALSE]
ord.inx <- order( unlist(res.mat[, 5]) )
res.mat <- res.mat[ord.inx, , drop = FALSE]

# dg=8;round(unlist(res.mat[,5]),dg)==round(unlist(res.mat.filter[,5]),dg)

####
