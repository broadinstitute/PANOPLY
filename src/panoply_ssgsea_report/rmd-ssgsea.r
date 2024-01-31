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
  opt_cmd <- parse_args( OptionParser(option_list=cmd_option_list) )
  
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

## ###################################################################
##      create a Rmarkdown report for sGSEA results
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output
rmd_ssgsea <- function(opt){    
  
  tar.file=opt$tar_file
  label=opt$label
  fdr=opt$fdr
  top.n=opt$top_n
  label.rmd <- 'report'
  
  #tmp.dir <- tempdir()
  tmp.dir <- 'tmp'
  wd <- getwd()
  
  ## prepare log file
  logfile=paste0(label, '_rmd-', label.rmd)
  start.time <- Sys.time()
  cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_', label.rmd,'\'--\n\n', file=logfile ) 
  cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
  
  
  ## #################################
  ## extract tar ball
  if(!dir.exists(tmp.dir))
    dir.create(tmp.dir)
  cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
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
  
  save(label, label.rmd, fdr, top.n, opt,
       gct.comb, param,
       fn.png,
       file=file.path( tmp.dir, 'data.RData'))                               
  
  ## ####################################################
  ## Rmarkdown
  cat('## Generating Rmarkdown file\n', file=logfile, append=T)
  

  rmd <- paste0("---
title: ssGSEA report - ", label,"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
---

## Overview

This document provides a high-level summary of the ```panoply_ssgsea``` module which performs single sample Gene Set Enrichment analysis (ssGSEA) on each column of the input data matrix. For more information about the module please visit the [PANOPLY Wiki page](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea). 

```{r echo=F, warning=F, message=F}\n
options(knitr.table.format = 'pipe', stringsAsFactors=F)
library(pacman)
p_load(plotly)
p_load(dplyr)
p_load(kableExtra)

source(file.path(opt$libdir, 'rmd-ssgsea-functions.R'))

## ##############################
## prepare data set
load('data.RData') ## import data 

# counter variable for figures and tables
fig_count <- 1
tab_count <- 1
```

\n\n
***
                
\n")
  
  
  rmd <- paste0( rmd, "\n
  
## Pathway heatmap

The heatmap in **Figure `r fig_count`** depicts normalized enrichment scores (NES) of gene sets (y-axis) across columns of the input data matrix (x-axis). NES were calculated according to the type of analysis used to derive the input data matrix:

* <u>NMF-clustering:</u> Columns represent NMF cluster. NES were derived from ssGSEA applied to feature weights determined by NMF. 

* <u>Consenus clustering:</u> Columns represent consensus cluster. NES were derived from ssGSEA applied to signed, -log transformed p-values of a two sample comparison  (one vs. other).

* <u>ssGSEA projection:</u> Columns represent individual samples. NES were derived from ssGSEA applied the expression data.

\n
\n

***

```{r, include=TRUE, fig.align='left', fig.cap=paste0('**Figure ', fig_count,'**: Normalized enrichment scores (NES) of gene sets across columns of the input data matrix. Asterisks indicate gene sets with FDR < ', fdr,'. Shown are up to ', top.n,' most significant gene sets per column.'), echo=FALSE}
knitr::include_graphics(fn.png)
```


```{r inc_fig_1, echo=F}
## increment
fig_count <- fig_count + 1
```\n

\n\n

***

", sep="")
  
  
  ######################################
  ## parameters
  rmd <- paste(rmd, "\n\n
\n## Parameters
\n```{r param, warning=F, message=F, results='as.is', echo=F}
kbl(param, col.names=' ')
\n```
\n", sep='')  
  
  
  
  ## ####################################
  ## Footer
  rmd <- paste(rmd, "\n\n
<br>
\n
\n***
\n
\n**_Created on ", Sys.time(),"_**
\n
\n***
\n<br>
\n", sep='')
  
  
  cat('## Rendering Rmarkdown file\n', file=logfile, append=T)
  writeLines(rmd, con=paste(tmp.dir, paste0( label.rmd, '_', label,'.rmd'), sep='/'))
  rmarkdown::render(paste(tmp.dir, paste(label.rmd, '_', label,'.rmd', sep=''), sep='/'))
  file.copy(paste(tmp.dir, paste(label.rmd, '_', label,'.html', sep=''), sep='/'), paste(wd, paste(label.rmd, '_', label,'.html', sep=''), sep='/'))
  
  cat('## Output written to:', paste(wd, paste(label.rmd, '_', label,'.html', sep=''), sep='/'), file=logfile, append=T)
  cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )
  
  
}


## ################################################
## run
rmd_ssgsea(opt)
