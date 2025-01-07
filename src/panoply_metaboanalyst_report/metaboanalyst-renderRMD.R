#
# Copyright (c) 2024 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for MetaboAnalyst Clustering module ###

args = commandArgs(TRUE)


rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### MetaboAnalyst Outputs ####
  make_option( c("-n", "--metaboanalyst_results"), action='store', type='character',  dest='metaboanalyst_results', help='Tar file containing figures and analyses for MetaboAnalyst results.'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='label', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list),
                   #' # for testing arguments
                   #' args = c('--metaboanalyst_results',"/opt/input/ODG_v3_prelim_MetaboAnalyst 2.tar.gz",
                   #'          '-x',"odg_test")
                   #'          #'-z',"/opt/input/")
                   )


# file.copy("/opt/input/metaboanalyst_rmd.rmd", file.path(opt$lib_dir,"metaboanalyst_rmd.rmd"), overwrite=T)

pwd=getwd()
fn = paste0(opt$label,"_metaboanalyst_report.html")
rmarkdown::render(file.path(opt$lib_dir,"metaboanalyst_rmd.rmd"),
# rmarkdown::render("/opt/input/metaboanalyst_rmd.rmd",
                  params = list(title = paste0("MetaboAnalyst Report - ", opt$label),
                                metaboanalyst_results = opt$metaboanalyst_results,
                                label = opt$label),
                  output_file = file.path(pwd,fn)
)

# file.copy(fn, '/opt/input', overwrite = T)
