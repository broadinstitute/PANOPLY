#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for NMF Clustering module ###

args = commandArgs(TRUE)


rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### NMF Outputs ####
  make_option( c("-n", "--nmf_results"), action='store', type='character',  dest='nmf_results', help='Tar file containing combined expression GCT, combined non-negative expression GCT, and Rdata object with NMF Clustering results (res.rank object).'),
  make_option( c("-r", "--rank_top"), action='store', type='numeric',  dest='rank_top', help='Best number of clustering / rank.'),
  #### Postprocessing Outputs ####
  make_option( c("-t", "--postprocess_tar"), action='store', type='character',  dest='postprocess_tar', help='Tar file containing figures and analyses for NMF results.'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='label', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list),
                   # # for testing arguments
                   # args = c('--nmf_results',"/opt/input/odg_all-mo_nmf_NMF_results.tar.gz",
                   #          '--rank_top',"3",
                   #          '-t',"/opt/input/odg_all-mo_nmf_NMF_postprocess.tar.gz",
                   #          '-x',"odg_test")
                   )


pwd=getwd()
rmarkdown::render(file.path(opt$lib_dir,"nmf_rmd.rmd"),
                  params = list(title = paste0("NMF Report - ", opt$label),
                                nmf_results = opt$nmf_results,
                                label = opt$label,
                                rank_top = opt$rank_top,
                                postprocess_tar = opt$postprocess_tar),
                  output_file = file.path(pwd,paste0(opt$label,"_nmf_report.html")))


