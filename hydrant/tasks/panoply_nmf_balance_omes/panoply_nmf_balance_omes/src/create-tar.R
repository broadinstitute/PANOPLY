#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
label <- args[1]

library(glue)

## function to create tar ball as input 
## for panoply_mo_nmf
create_tar <- function( 
  label='test',               ## character, used in filename of tar file
  str_list=list(
    prot='prote_ome.txt',     ## text file containing the absolute file path to PROTEOME GCT file
    pSTY='phospho_ome.txt',   ## text file containing the absolute file path to PHOSPHO GCT file
    acK='acetyl_ome.txt',     ## text file containing the absolute file path to ACETYL GCT file
    ubK='ubiquityl_ome.txt',  ## text file containing the absolute file path to UBIQUITYL GCT file
    glyco="glyco_ome.txt",    ## text file containing the absolute file path to GLYCO GCT file
    RNA='rna_ome.txt',        ## text file containing the absolute file path to RNA GCT file
    CNV='cna_ome.txt'         ## text file containing the absolute file path to CNV GCT file
    )
){
  
  ## non-empty files
  keep_idx <- which(sapply(str_list, function(x) ifelse(nchar(readLines(x)) > 0, TRUE, FALSE)))
  str_list <- str_list[ keep_idx ]

  ## file paths
  file_list <- sapply(str_list, readLines)
  
  ## copy locally
  file.copy(file_list, '.')
  
  ## update paths
  file_list <- sapply(file_list, function(x) sub('.*/', '', x))
  
  ## create mapping file
  file.create('nmf.conf')
  writeLines(paste(names(file_list), file_list, sep='\t'), con='nmf.conf')
  
  ## create tar archive
  #fn <- tolower( paste(names(file_list), collapse='-') )
  fn <- glue("{label}.tar")
  system( glue("tar -cvf {fn} nmf.conf {paste(file_list, collapse=' ')}" ))

  return(0)
}

## run
create_tar(label)

