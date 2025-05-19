#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("optparse"))

options( warn = -1, stringsAsFactors=F )

# specify command line arguments
option_list <- list(
  make_option( c("-t", "--tar_file"), action='store', type='character',  dest='tar_file', help='Path to panoply_ssgsea result file.'),
  make_option( c("-l", "--label"), action='store', type='character',  dest='label', help="label/file prefix used in ssgsea GCT files.", default='NA'),
  make_option( c("-f", "--fdr"), action='store', type='character',  dest='fdr', help="max. FDR", default='NA'),
  make_option( c("-n", "--top_n"), action='store', type='character',  dest='top_n', help="Max. number of significant hits to plot/label", default='NA'),
  make_option( c("-y", "--yaml_file"), action='store', type='character',  dest='yaml_file', help='yaml parameter file.', default = NA),
  make_option( c("-p", "--ptmsigdb"), action='store', type='logical',  dest='ptmsigdb', help='PTMsigDB?', default = FALSE),
  make_option( c("-z", "--libdir"), action='store', type='character',  dest='libdir', help='Folder to source from.', default = 'NA')
)

################################################
## funtion to parse parse and update parameters
## - cmd line
## - yaml file
## parameters in yaml file will be updated with 
## parameters specified on cmd
parse_param_ssgsea_report <- function(cmd_option_list, yaml_section='panoply_ssgsea_report'){
  
  ## #########################################################
  # parse command line parameters
  opt_cmd <- parse_args( OptionParser(option_list=option_list) ,
                         # ## optional testing arguments
                         # args = c("--tar_file", "/opt/input/ODG_v3-NMF.k3.core-proteome-KEGGpathways.tar.gz",
                         #          "--label", "KEGG_MetaboAnalyst",
                         #          "--yaml_file", "/opt/input/master-parameters.yaml",
                         #          "--libdir", "/script_source")
                         )
  
  ############################################################
  ## parse yaml file
  if(!is.na(opt_cmd$yaml_file)) {
    
    if(file.exists(opt_cmd$yaml_file)){
      
      p_load(yaml)
      
      ## import yaml
      opt_yaml <- read_yaml(opt_cmd$yaml_file)
      
      ## extract relevant section
      opt_yaml <- opt_yaml[[yaml_section]]
      
      ## parse cmd params
      cat('\n\nparsing command line parameters:\n', paste0(rep('-', 60), collapse = ''), '\n', sep='')
      for(x in names(opt_cmd))
        cat('---', x, opt_cmd[[x]], '; prefer yaml?', opt_cmd[[x]] == 'NA','\n')
      
      ## update yaml with parameters specified on cmd line 
      cat('\n\nUpdating parameter file with command line parameters:\n', paste0(rep('-', 60), collapse = ''), '\n', sep='')
      ## cmd parameters
      cmd_not_null <- which( !sapply(opt_cmd, function(x) x == 'NA' ) )
      cmd_to_update <- intersect( names(opt_cmd)[ cmd_not_null], names(opt_yaml) )
      cmd_to_add <- setdiff( names(opt_cmd), names(opt_yaml) )
      
      
      sapply(cmd_to_update, function(x) cat(x, ':', opt_yaml[[x]], '->', opt_cmd[[x]], '\n'))
      
      ## update yaml by cmd 
      opt_yaml[cmd_to_update] <- opt_cmd[cmd_to_update]
      
      cat(paste0(rep('-', 60), collapse = ''), '\n\n')
      
      ## add parameters only specified on cmd
      if(length(cmd_to_add) > 0){
        opt_cmd_to_add <- opt_cmd[cmd_to_add]
        opt_yaml <- append(opt_yaml, opt_cmd_to_add)
      }
      
      ## updated params
      opt <- opt_yaml
      
    }    
  } else {
    ## no yaml file
    opt <- opt_cmd
  }
  
  #########################################
  ## force correct mode
  opt$tar_file <- as.character(opt$tar_file)
  opt$label <- as.character(opt$label)
  opt$fdr <- as.numeric(opt$fdr)
  opt$top_n <- as.numeric(opt$top_n)
  opt$ptmsigdb <- as.logical(opt$ptmsigdb)
  opt$libdir <- as.character(opt$libdir)
  
  return(opt)
} 


# parse command line parameters
opt <- parse_param_ssgsea_report(option_list)
source(file.path(opt$libdir, 'rmd-ssgsea-functions.R'))


#require(pacman)
p_load(rmarkdown)
p_load(cmapR)
p_load(dplyr)
p_load(glue)
p_load(knitr)
p_load(kableExtra)



#### Create Figures & Rdata object for RMD repoort ####

tar.file=opt$tar_file
label=opt$label
fdr=opt$fdr
top.n=opt$top_n

#tmp.dir <- tempdir()
tmp.dir <- 'tmp'
wd <- getwd()

# ## prepare log file
# logfile=paste0(label, '_rmd-', label.rmd)
# start.time <- Sys.time()
# cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_', label.rmd,'\'--\n\n', file=logfile )
# cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)


## #################################
## extract tar ball
if(!dir.exists(tmp.dir))
  dir.create(tmp.dir)
# cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
untar(tar.file, exdir=tmp.dir)

#####################################
## identify -combinded.gct
gct.comb <- dir(tmp.dir, pattern = '-combined.gct$', full.names = T)

## parameter file
param <- dir(tmp.dir, pattern = 'parameters.txt$', full.names = T)

######################################
## import and save as .Rdata
if(file.exists(gct.comb))
  gct.comb <- parse.gctx(gct.comb)
if(file.exists(param))
  param <- readLines(param)

## create the figure
pw_hm(gct.comb, fdr.max=fdr, n.max=top.n, ptmsigdb=opt$ptmsigdb)

## copy to tmp.dir
fn.png <- dir('.', pattern='.png$')
file.copy(fn.png, tmp.dir)



#### Generate RMD Report ####

pwd=getwd()
rmarkdown::render(file.path(opt$libdir, "ssgsea_rmd.rmd"),
                  #"/script_source/ssgsea_rmd.rmd",
                  params = list(title = paste0("ssGSEA Report - ", label),
                                label = label,
                                fdr = fdr,
                                top.n = top.n,
                                param = param,
                                fn.png = file.path(pwd, fn.png)),
                  output_file = file.path(pwd,paste0(label,"_ssGSEA_rmd.html")))

# file.copy(file.path(pwd,paste0(label,"_ssGSEA_rmd.html")),
#           '/opt/input/', overwrite=T)
