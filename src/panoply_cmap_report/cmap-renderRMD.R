#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for CMAP module ###

args = commandArgs(TRUE)


rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### CMAP Analysis Outputs ####
  make_option( c("-c", "--cmap_results"), action='store', type='character',  dest='cmap_results', help='Tar file containing results from the CMAP annotate module.'),
  make_option( c("-s", "--cmap_ssgsea_results"), action='store', type='character',  dest='cmap_ssgsea_results', help='Tar file containing ssGSEA analysis run on the CMAP results.'),
  #### General Parameters ####
  make_option( c("-l", "--output_prefix"), action='store', type='character',  dest='label', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list),
                   # # for testing arguments
                   # args = c('-c',"/opt/input/panoply_cmap-annotate-output.tar",
                   #          '-s',"/opt/input/panoply_cmap_annotate-ssgsea.tar",
                   #          '-x',"odg_test")
                   )

source(file.path(opt$lib_dir, 'rmd-ssgsea-functions.R')) # source pw_hm() function

pwd=getwd()
rmarkdown::render(file.path(opt$lib_dir,"cmap_rmd.rmd"),
                  params = list(title = paste0("CMAP Report - ", opt$label),
                                cmap_results = opt$cmap_results,
                                cmap_ssgsea_results = opt$cmap_ssgsea_results,
                                label = opt$label),
                  output_file = file.path(pwd,paste0(opt$label,"_CMAP_report.html")))


