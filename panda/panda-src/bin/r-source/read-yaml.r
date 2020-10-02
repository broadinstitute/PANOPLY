#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( yaml )
p_load( optparse )
p_load( glue )

opt <- list()

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-i", "--input"  ),
      dest   = "input",
      action = "store",
      type   = 'character' ) )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

main <- function(){
  set_arguments()
  yaml  <- read_yaml( opt$input )
  shell <- file( "config.sh" )
  lines <- "#!/bin/bash"

  lines <- c( lines, "\n## Location of input files" )
  lines <- c( lines, glue(
    "freeze_path=\"{yaml$freeze.path}\"\n" ) )

  lines <- c( lines, "\n## Columns to select for groups.csv file" )
  lines <- c( lines, glue(
    "groups_cols=\"{paste0( yaml$groups.cols, collapse = ';' )}\"\n" ) )


  yaml$typemap.gct <- yaml$typemap.gct[lapply( yaml$typemap.gct, length ) > 0]
  yaml$typemap.csv <- yaml$typemap.csv[lapply( yaml$typemap.csv, length ) > 0]
  gct.files <- paste0( unlist( yaml$typemap.gct ), collapse = ';' )
  csv.files <- paste0( unlist( yaml$typemap.csv ), collapse = ';' )
  gct.types <- paste0( names ( yaml$typemap.gct ), collapse = ';' )
  csv.types <- paste0( names ( yaml$typemap.csv ), collapse = ';' )


  lines <- c( lines, "\n## Input files and extensions" )
  lines <- c( lines, glue(
    "gct_files=\"{gct.files}\"\n" ) )
  lines <- c( lines, glue(
    "csv_files=\"{csv.files}\"\n" ) )
  lines <- c( lines, glue(
    "gct_types=\"{gct.types}\"\n" ) )
  lines <- c( lines, glue(
    "csv_types=\"{csv.types}\"\n" ) )

  lines <- c( lines, "\n## Terra parameters" )
  wkspace <- gsub( pattern = ".", replacement = "-", x = yaml$wkspace )
  lines <- c( lines, glue(
    "wkspace=\"{yaml$wkspace}\"\n" ) )
  lines <- c( lines, glue(
    "project=\"{yaml$globals$project}\"\n" ) )
  lines <- c( lines, glue(
    "group=\"{yaml$globals$group}\"\n" ) )
  lines <- c( lines, glue(
    "meth_space=\"{yaml$globals$meth_space}\"\n" ) )

  lines <- c( lines, "\n## Sample set parameters" )
  sets <- names( yaml$sample.sets )
  setline <- ""
  for ( x in sets ) {
    setline <- glue( "{setline}{x}" )
    setline <- glue( "{setline},{yaml$sample.sets[[x]]$fil_col}" )
    setline <- glue( "{setline},{yaml$sample.sets[[x]]$fil_val}" )
    ## for the next line to work the config.yaml file must have
    ## for each sample set, two entries, fil_col and fil_val before
    ## adding any other attributes
    other.attr <- names( yaml$sample.sets[[x]] )[-c( 1, 2 )]
    for ( attr in other.attr ) {
      setline <- glue( "{setline},{attr}@" )
      setline <- glue( "{setline}{yaml$sample.sets[[x]][[attr]]}" )
    }
    setline <- glue( "{setline}%")
  }

  setline <- paste0( unlist( strsplit( setline, split = '%' ) ),
                     collapse = '%' )
  setline <- glue( "sets=\"{setline}\"" )
  lines <- c( lines, setline )

  writeLines( lines, shell )
  close( shell )
}

if ( !interactive() ){
  main()
}
