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
  make_option( c("-f", "--fdr"), action='store', type='double',  dest='fdr', help="max. FDR"),
  make_option( c("-n", "--top_n"), action='store', type='integer',  dest='top_n', help="Max. number of significant hits to plot/label"),
  make_option( c("-c", "--cluster_rows"), action='store', type='logical',  dest='cluster_rows', help="If TRUE, rows will be clustered using distance matrix, with ser.meth as seriation method."),
  make_option( c("-m", "--ser_meth"), action='store', type='character',  dest='ser_meth', help="Seriation method used for clustering rows. The default in pw_hm() is 'ARSA'"),
  make_option( c("-g", "--geneset_groups_file"), action='store', type='character',  dest='geneset_groups_file', help="CSV file containing a mapping between all genesets, and some category they should be grouped into in the heatmap."),
  # make_option( c("-p", "--ptmsigdb"), action='store', type='logical',  dest='ptmsigdb', help='PTMsigDB?', default = FALSE),
  make_option( c("-s", "--split_by_prefix"), action='store', type='logical',  dest='split_by_prefix', help="If TRUE, separate heatmaps will be created for pathways with unique prefixes '<prefix>-<pathway_name>'. On by default for PTM-SEA pathways; can be manually overridden with this parameter."),
  make_option( c("-y", "--yaml_file"), action='store', type='character',  dest='yaml_file', help='yaml parameter file.', default = NA),
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
                         # args = c("--tar_file", "/opt/input/ODG_v3_proteome.tar.gz",
                         #          "--label", "ODG_v3",
                         #          "--yaml_file", "/opt/input/master-parameters.yaml",
                         #          # "--geneset_groups_file", "/opt/input/mitocarta_geneset_groups.csv",
                         #          "--libdir", "/home/pgdac/src/")
                         #          # "--libdir", "/Users/wcorinne/Git/panoply-sandbox/src/panoply_ssgsea_report")
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
      
    } else {
      warning("WARNING: YAML file does not exist. Using command-line parameters only.")
      ## no yaml file
      opt <- opt_cmd
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
  opt$cluster_rows <- as.logical(opt$cluster_rows)
  opt$ser_meth <- as.character(opt$ser_meth)
  if(!is.null(opt$geneset_groups_file)) opt$geneset_groups_file <- as.character(opt$geneset_groups_file) # either retain NULL, or coerce character
  # opt$ptmsigdb <- as.logical(opt$ptmsigdb)
  if(!is.null(opt$split_by_prefix)) opt$split_by_prefix <- as.logical(opt$split_by_prefix) # either retain NULL, or coerce character
  
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

tar_file=opt$tar_file
label=opt$label
fdr.max=opt$fdr
n.max=opt$top_n
cluster.rows <- opt$cluster_rows; if (length(cluster.rows)==0) cluster.rows=F # if we didn't set cluster_rows, set to FALSE (legacy behavior)
ser.meth <- opt$ser_meth
geneset_groups_file=opt$geneset_groups_file
split.by.prefix=opt$split_by_prefix

tmp.dir <- tempdir()
# tmp.dir <- 'tmp'
wd <- getwd()

# ## prepare log file
# logfile=paste0(label, '_rmd-', label.rmd)
# start.time <- Sys.time()
# cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_', label.rmd,'\'--\n\n', file=logfile )
# cat('## parameters\ntar file:', tar_file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)


## #################################
## extract tar ball
if(!dir.exists(tmp.dir))
  dir.create(tmp.dir)
# cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
untar(tar_file, exdir=tmp.dir)

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
pw_hm(output.prefix=gct.comb, fdr.max=fdr.max, n.max=n.max,
      geneset_groups_file = geneset_groups_file,
      split.by.prefix = split.by.prefix,#, ptmsigdb=opt$ptmsigdb)
      cluster.rows = cluster.rows, ser.meth = ser.meth)

## copy to tmp.dir
fn.png <- dir('.', pattern='.png$')
if ( length(fn.png)==0 ) stop("No heatmap figures were created. Something likely went wrong in pw_hm().")
file.copy(fn.png, tmp.dir)



#### Generate RMD Report ####

pwd=getwd()
rmarkdown::render(file.path(opt$libdir, "ssgsea_rmd.rmd"),
                  #"/script_source/ssgsea_rmd.rmd",
                  params = list(title = paste0("ssGSEA Report - ", label),
                                label = label,
                                fdr = fdr.max,
                                top_n = n.max,
                                param = param,
                                fn.png = file.path(pwd, fn.png)),
                  output_file = file.path(pwd,paste0(label,"_ssGSEA_rmd.html")))

# file.copy(file.path(pwd,paste0(label,"_ssGSEA_rmd.html")),
#           '/opt/input/', overwrite=T)
