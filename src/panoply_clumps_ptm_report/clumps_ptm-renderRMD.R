#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for clumps_ptm Clustering module ###

args = commandArgs(TRUE)


rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))


#### Command Line Arguments ####
option_list <- list(
  #### clumps_ptm Outputs ####
  make_option( c("-i", "--clumps_ptm_results"), action='store', type='character',  dest='clumps_ptm_results', help='Tar file(s) containing figures and analyses for ClumpsPTM results, comma-separated.'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='label', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action="store", dest='lib_dir', type="character", help="the src directory.", default='/prot/proteomics/Projects/PGDAC/src')
  #### ####
)

opt <- parse_args( OptionParser(option_list=option_list),
                   #' # for testing arguments
                   #' args = c('--clumps_ptm_results',"/opt/input/ODG_v3_NMF.consensus.core.k3_clumps_ptm_full_results.tar,/opt/input/ODG_v3_NMF.consensus.core.k6_clumps_ptm_full_results.tar",
                   #'          '-x',"ODG_v3")
                   #'          #'-z',"/opt/input/")
                   )


# file.copy("/opt/input/clumps_ptm_rmd.rmd", file.path(opt$lib_dir,"clumps_ptm_rmd.rmd"), overwrite=T)

pwd=getwd()
fn = paste0(opt$label,"_clumps_ptm_report.html")
rmarkdown::render(file.path(opt$lib_dir,"clumps_ptm_rmd.rmd"),
# rmarkdown::render("/opt/input/clumps_ptm_rmd.rmd",
                  params = list(title = paste0("Clumps-PTM Report - ", opt$label),
                                clumps_ptm_results = opt$clumps_ptm_results,
                                label = opt$label),
                  output_file = file.path(pwd,fn)
)

# file.copy(fn, '/opt/input', overwrite = T)
