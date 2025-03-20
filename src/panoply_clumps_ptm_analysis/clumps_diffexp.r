#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#
rm(list=ls())
options( warn = -1, stringsAsFactors = F )
suppressPackageStartupMessages(library("optparse"))

# specify command line arguments
option_list <- list(
  make_option( c("-p", "--phosphoproteome_gct"), action='store', type='character',  dest='gct_pSTY', help='Path to input phosphoproteome GCT file.'),
  make_option( c("-a", "--acetylome_gct"), action='store', type='character',  dest='gct_acK', help='Path to input acetylome GCT file.'),
  make_option( c("-u", "--ubiquitylome_gct"), action='store', type='character',  dest='gct_ubK', help='Path to input ubiquitylome GCT file.'),
  make_option( c("-g", "--groups_file"), action='store', type='character',  dest='groups_file', help='Groups-file, i.e. an annotations file subsetted to annotations of interest. If not provided, all annotations in the cdesc will be analyzed.'),
  make_option( c("-c", "--gene_column"), action='store', type='character', dest='gene_col', help='Column name in rdesc in the GCT that contains gene names.'), # default='geneSymbol'),
  make_option( c("-f", "--fdr_assoc"), action='store', type='character', dest='fdr_assoc', help='FDR value for marker selection.'), # default='geneSymbol'),
  #### General Parameters ####
  make_option( c("-x", "--output_prefix"), action='store', type='character',  dest='output_prefix', help='Label associated with this run.'),  # default = 2),
  make_option( c("-z", "--libdir"), action='store', type='character',  dest='libdir', help='Folder to source from.'),
  make_option( c("-y", "--yaml_file"), action='store', type='character',  dest='yaml_file', help='yaml parameter file.', default = 'NA')
)


#### Parse Command-Line Arguments ####
opt_cmd <- parse_args( OptionParser(option_list=option_list),
                       # for testing arguments
                       args = c('--phosphoproteome_gct',"opt/input/ODG-v3-phosphoproteome-SpectrumMill-ratio-QCfilter-NArm.gct",
                                '--acetylome_gct',"opt/input/ODG-v3-acetylome-SpectrumMill-ratio-QCfilter-NArm.gct",
                                '--ubiquitylome_gct',"opt/input/ODG-v3-ubiquitylome-SpectrumMill-ratio-QCfilter-NArm.gct",
                                '-g',"opt/input/groups-subset.csv",
                                '-y',"opt/input/master-parameters.yaml",
                                '-x',"ODG_v3")
)


#### Parse YAML Arguments ####
opt = opt_cmd # initialize options with command line options
if ( !is.null(opt$yaml_file) ) {
  #### read in yaml ####
  library(yaml)
  yaml_out <- read_yaml(opt$yaml_file)
  #### overwrite the command-line parameters ####
  # global parameters
  if (is.null(opt$gene_col)) opt$gene_col = yaml_out$global_parameters$gene_mapping$gene_id_col
  if (is.null(opt$sample_na_max)) opt$sample_na_max = yaml_out$global_parameters$missing_values_and_filtering$sample_na_max
  # config.r parameters
  if (is.null(opt$id_col)) opt$id_col = yaml_out$DEV_sample_annotation$gct_file_ids$phosphoproteome$id_col
  if (is.null(opt$desc_col)) opt$desc_col = yaml_out$DEV_sample_annotation$gct_file_ids$phosphoproteome$desc_col
  # diff-exp parameters
  if (is.null(opt$fdr_assoc)) opt$fdr_assoc = yaml_out$panoply_association$fdr_assoc # toDo: consider changing to a clumps-ptm parameter
} else { # if no YAML was provieded
  # check if any necessary parameters are missing
  if( any(sapply(list(opt$gene_col), 
                 is.null)) ) { # if we have at least one missing parameter
    stop("Master Parameter yaml-file is missing. Please either provide a master-parameters file, or manually provide all NMF parameters.") # error and stop
  }
}


library(tidyverse)
library(glue)
library(cmapR)


#### Taken from config.R ####
## Path for R-utilites (for I/O and other misc functions)
Rutil.path <- switch (Sys.info()[['sysname']],
                      Windows = {'//flynn-cifs/prot_proteomics/Projects/R-utilities'},
                      Darwin = {'/Volumes/prot_proteomics/Projects/R-utilities'},
                      Linux = {'/prot/proteomics/Projects/R-utilities'})

Source <- function (f) {
  # adaptive source file location -- from current directory or from R-utiliites
  if (file.exists (f)) {source (f)}
  else {
    f.os <- file.path (Rutil.path, f)
    if (file.exists (f.os)) source (f.os)
    else stop (paste ("Can't find file", f.os, '-- R-utilities missing?'))
  }
}

Source ('gct-io.r')
Source ('io.r')

#### Source Files ####

Source ('markersel-classify.r')
Source ('stats.r')


#### Parameter Wrangling ####

if ( is.null(opt$gct_pSTY) && is.null(opt$gct_ubK) && is.null(opt$gct_acK) ) stop("No GCT files provided. Please provide at least one PTM GCT file.")

# organize parameters for use
assoc.subgroups=opt$groups_file
assoc.fdr = opt$fdr_assoc
# config.R params
gene.id.col=opt$gene_col
sample.na.max = opt$sample_na_max
id.col=opt$id_col
desc.col=opt$desc_col
ome_list = list(list(type = 'phosphoproteome', gct = opt$gct_pSTY),
                list(type = 'acetylome', gct = opt$gct_acK),
                list(type = 'ubiquitylome', gct = opt$gct_ubK))

#### Differential Expression Analysis ####

# NOTE: this code is largely based on assoc-analysis.r
# must be used with the association docker from PANOPLY v1_6 or later
# v1_5 and earlier had a bug where logFC values were not calculated in the correct direction
for (ome in ome_list) { # use variable t instead of type, since config.r con
  # MODIFIED: removed original method of retrieving GCT
  type = ome$type
  gct.file = ome$gct
  
  
  run.marker.selection <- function (input.gct.file, input.cls.file, prefix, run.1vAll=FALSE) {
    # runs marker selection on input data for given class vector in input.cls
    tryCatch (marker.selection.and.classification (input.gct.file, input.cls.file, paste (prefix, '-analysis', sep=''),
                                                   gene.id.col=gene.id.col, id.col=id.col, desc.col=desc.col, gsea=FALSE,
                                                   id.to.gene.map=NULL,   # GeneSymbol already present in GCT v1.3 input
                                                   duplicate.gene.policy=duplicate.gene.policy,
                                                   impute.colmax=sample.na.max,
                                                   official.genenames=file.path ('..', 'data', 'gene-symbol-map.csv'),
                                                   fdr=assoc.fdr,
                                                   models=c("pls","rf","glmnet")),
              error = function(cond) {
                message(paste("Failed to complete marker selection for ", prefix))
                message(cond)
              })
    
    if (run.1vAll) {
      # if class vector has > 2 classes, and sufficient numbers per class,
      # run 1 vs. all marker selection for each class
      cls <- read.cls (input.cls.file)
      if (nlevels (factor (cls)) > 2 && min (summary (factor (cls))) >= 5) { ### MODIFIED: reduced min samples per subvalues
        cls.1vAll <- classes.1vAll (cls, sort_other_first=TRUE)
        for (i in 1:ncol(cls.1vAll)) {
          prefix.1vA <- paste (prefix, colnames(cls.1vAll)[i], sep='-')
          new.clsf <- paste (prefix.1vA, '.cls', sep='')
          write.cls (cls.1vAll[,i], new.clsf, respect_factor_order = TRUE)
          marker.selection.and.classification (input.gct.file, new.clsf, paste (prefix.1vA, '-analysis', sep=''),
                                               gene.id.col=gene.id.col, id.col=id.col, desc.col=desc.col, gsea=FALSE,
                                               id.to.gene.map=NULL,   # GeneSymbol already present in GCT v1.3 input
                                               duplicate.gene.policy=duplicate.gene.policy,
                                               impute.colmax=sample.na.max,
                                               official.genenames=file.path ('..', 'data', 'gene-symbol-map.csv'),
                                               fdr=assoc.fdr,
                                               models=c("pls","rf","glmnet"))
        }
      } else {
        print(glue::glue("\n\nWARNING: only {nlevels (factor (cls))} levels and min {min (summary (factor (cls)))} samples per subvalue.\n\n"))
      }
    }
  }
  
  # this module ASSUMES that a groups file (assoc.subgroups) is specified
  if (! exists ("assoc.subgroups")) {
    stop("A groups file must be provided. This file should include Sample.ID,
       additional row-description columns to be analyzed for enrichent.")
  } else {
    # (groups file format is similar to expt-design-file, with Sample.ID and additional columns)
    # association analysis will be run for each additional column, excluding samples marked 'ignore'
    # (different columns cannot have the same subgroup name)
    subgroup.table <- read.csv (assoc.subgroups)
    rownames (subgroup.table) <- subgroup.table [, 'Sample.ID']
    cls.list <- setdiff (colnames (subgroup.table), 'Sample.ID')
    
    ds <- parse.gctx (gct.file)
    sample.order <- ds@cid
    
    if (length (cls.list) > 0) {
      for (g in cls.list) {
        #replace "" with NA, otherwise it gets kept
        group <- subgroup.table[sample.order, g]
        group[group==""]=NA
        group <- make.names (group)  # use make.names to convert text to proper class labels
        #REMOVE NAs!
        subsamp <- !group%in%c('ignore',"NA.")
        # write subsamples dataset and class labels
        f <- paste (type, '-', g, sep='')
        ds.g <- col.subset.gct (ds, subsamp)
        gct.g <-  sprintf ("%s.gct", f)
        cls.g <- sprintf ("%s.cls", f)
        write.gct (ds.g, gct.g, appenddim=FALSE)
        write.cls (group [subsamp], cls.g)
        
        if (min (summary (factor (read.cls (cls.g)))) < 3) {
          # too few members in some class(es)
          warning ( paste (g, "has classes with < 3 members ... skipping") )
          next
        }
        run.marker.selection (gct.g, cls.g, f, run.1vAll=TRUE) ### MODIFIED: run.1vAll=TRUE
        
        ### MODIFIED: compile and save results into a single CSV
        {
          files = list.files(pattern=glue('^{type}-{g}-class.+-analysis-markers-all.csv'))
          fc.full = data.frame(id=character(0))
          p.full = data.frame(id=character(0))
          adj.p.full = data.frame(id=character(0))
          t.score.full = data.frame(id=character(0))
          df.full = data.frame(id=character(0))
          for (fn in files) {
            class = gsub(".+-class.(.+)-analysis-markers-all.csv", "\\1", fn)
            df.tmp = read.csv(fn) %>%
              dplyr::rename(id = Gene.ID)
            # combine into one full table
            df.full = full_join(df.full,
                                dplyr:: select(df.tmp, c(id, Fold.Change, P.Value, adj.P.Val, limma.score)) %>%
                                  dplyr::rename( t.score = limma.score ) %>% # rename limma.score to t.score
                                  dplyr::rename_at(2:5, ~ paste0(., glue(".{class}")) ), 
                                by = 'id')
          }
          
          write.csv(df.full, glue('{type}-{g}-1vAll-Fold_Change_and_P_Val.csv'))
        }
      }
      
    } else {
      warning ("No cls files to run association analysis on ... ignoring")
    }
  }
}