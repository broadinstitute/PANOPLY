#!/usr/bin/env Rscript
options( warn = -1 )
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("pacman"))

## 20180309
## preprocess GCT file to:
## 1. make sure ids are unique
## 2. create PTM site-centric GCT reports

# parse the directory this file is located
#this.file.dir <- commandArgs()[4]
#this.file.dir <- sub('^(.*(/|\\\\)).*', '\\1', sub('.*?\\=','', this.file.dir))
#cat('DIRECTORY: ',this.file.dir, '\n')

# specify command line arguments
option_list <- list(
  make_option( c("-i", "--input"), action='store', type='character',  dest='gct.str', help='Path to input GCT file.'),
  make_option( c("-l", "--level"), action='store', type='character',  dest='level', help='Mode of report, "site" or "gene".', default='site'),
  make_option( c("-r", "--residue"), action='store', type='character',  dest='mod.res', help='Modified residues, e.g. "S|T|Y" or "K".', default='S|T|Y'),
  make_option( c("-m", "--modification"), action='store', type='character',  dest='mod.type', help='Type of modification, e.g "-p" or "-ac".', default = '-p'),
  make_option( c("-a", "--accession"), action='store', type='character',  dest='acc.type', help='Type of protein accession number.', default = 'uniprot'),
  make_option( c("-o", "--organism"), action='store', type='character',  dest='org', help='Organism.', default = 'human'),
  make_option( c("-d", "--dimension"), action='store', type='logical',  dest='appenddim', help='Logical, should matrix dimensions be included in file name?', default = TRUE),
  make_option( c("-z", "--libdir"), action='store', type='character',  dest='libdir', help='folder to source from', default = 'c:/Users/karsten/Dropbox/Devel/PGDAC/src/pgdac_ssgsea/')
  
  )
# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

p_load(doParallel)
p_load(foreach)
p_load(glue)

source(glue("{opt$libdir}/pgdac_ptmsea_functions.R"))


###################################################
##       run the function
preprocessGCT(gct.str = opt$gct.str,
              level=opt$level,
              mode=opt$mode,
              mod.res = opt$mod.res,
              mod.type=opt$mod.type,
              acc.type=opt$acc.type,
              org=opt$org,
              appenddim=opt$appenddim)
