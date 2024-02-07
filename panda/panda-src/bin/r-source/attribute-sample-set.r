#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( glue )

sample.set.member.file <- "sample_set_membership.tsv"
opt <- list()
acc.patterns <- '\\.csv|\\.gct'

set_arguments <- function() {
  optionList <- list(
    make_option(
      c( "-a", "--bucket" ),
      dest = "bucket",
      action = "store",
      type = 'character'),

    make_option(
      c( "-b", "--wk.space" ),
      dest = "wk.space",
      action = "store",
      type = 'character'),

    make_option(
      c( "-c", "--project" ),
      dest = "project",
      action = "store",
      type = 'character'),

    make_option(
      c( "-d", "--csv.types" ),
      dest = "csv.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-e", "--gct.types" ),
      dest = "gct.types",
      action = "store",
      type = 'character') )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
  opt$agg.suffix <<- "ss"
}

read_sets <- function( sample.set.member.file ){
  
  meta.file <- read.csv(glue( "{sample.set.member.file}" ), sep="\t") # sample-set membership file
  set.names = unique(meta.file[[1]]) # get list of sets from the first column
  
  #  convert the meta.file into a list
  sets_list = sapply( set.names, function(set) { # for each set
    meta.file[meta.file[[1]]==set,2] # return a vector with the samples in that set
  }, simplify = FALSE, USE.NAMES = TRUE) # return as list, and label elements with set.names 
  
  return( sets_list )
}


process_other_attributes <- function( set ){
  ## Read other attributes for the sample sets such that
  ## other.files is a data frame with col.names and file.names
  ss.other.file.name <- glue( "other-attributes/{set}-other.txt" )

  if( !( file.exists( ss.other.file.name ) ) )
    return()

  other.files <- read.delim( ss.other.file.name,
                             header = F,
                             sep = '=',
                             stringsAsFactors = F )
  colnames( other.files ) <- c( 'col.name', 'file.name' )

  ## iterate over all the files/attributes per sample-set
  ## copy the associated file to the directory structure in the
  ## sample-set folder
  for ( rIdx in 1:nrow( other.files ) )
  {
    attr.col.name  <- other.files$col.name[rIdx]
    buck.file.name <- other.files$file.name[rIdx]
    locl.file.name <- tail( unlist( strsplit( buck.file.name,
                                              split = '/' ) ), 1 )
    file.copy( from = glue( "pipeline-input/{buck.file.name}" ),
               to = glue( "subsets/{set}/{locl.file.name}" ) )
  }

  ## use gsutil to copy all the subsets and sample-set
  ## related data files to the google bucket
  bucket.cp <- glue( 'gsutil -m cp subsets/{set}/* ',
                     'gs://{opt$bucket}/sample_sets/{set}/' )
  system( bucket.cp )

  ## iterate over "other" attributes once again and set attribute column
  ## on Terra's sample-set table to the corresponding google bucket file
  for ( rIdx in 1:nrow( other.files ) )
  {
    attr.col.name <- other.files$col.name[rIdx]
    buck.file.name <- other.files$file.name[rIdx]
    locl.file.name <- tail( unlist( strsplit( buck.file.name,
                                              split = '/' ) ), 1 )
    attr.set <- glue(
      'fissfc -V -y attr_set -w {opt$wk.space} -p {opt$project}',
      ' -a {attr.col.name}',
      ' -v gs://{opt$bucket}/sample_sets/{set}/{locl.file.name}',
      ' -t sample_set -e {set}' )
    system( attr.set )
  }
}

process_additional_parameters <- function( set ){
  ## use gsutil to copy all the additional parameter files to the
  ## google bucket
  if( length( list.files( glue( "additional-params/{set}/" ) ) ) > 0 ){
    bucket.cp  <- glue( 'gsutil -m cp additional-params/{set}/* ',
                        'gs://{opt$bucket}/sample_sets/{set}/' )
    system( bucket.cp )
  }

  types <- c( opt$csv.types, opt$gct.types )
  for ( type in types ){
    if ( type %in% c( 'cna', 'rna', 'annotation', 'groups' ) ) next

    ## set links between additional parameters for -ome
    get.addp <- glue(
      'gsutil ls gs://{opt$bucket}/sample_sets/{set}/{type}.r' )
    addp.path <- system( get.addp, intern = T )
    if ( length( addp.path ) == 0 ) next
    attr.set <- glue(
      'fissfc -V -y attr_set -w {opt$wk.space} ',
      '-p {opt$project} ',
      '-a {type}_add_params_ss ',
      '-v gs://{opt$bucket}/sample_sets/{set}/{type}.r ',
      '-t sample_set -e {set}' )
    system( attr.set )
  }
}

process_types <- function( set ){
  ## use gsutil to copy files to the google bucket
  bucket.cp <- glue( 'gsutil -m cp subsets/{set}/* ',
                     'gs://{opt$bucket}/sample_sets/{set}/' )
  system( bucket.cp )
  
  ## iterate over all types and set links between type files in the google
  ## bucket and in the Terra table
  types <- c( opt$csv.types, opt$gct.types )
  for ( type in types )
  {
    attr.name <- glue( "{type}_{opt$agg.suffix}" )
    if ( type %in% opt$gct.types ) ext <- 'gct' else ext <- 'csv'

    attr.set <- glue(
      'fissfc -V -y attr_set -w {opt$wk.space} ',
      '-p {opt$project} ',
      '-a {attr.name} ',
      '-v gs://{opt$bucket}/sample_sets/{set}/{type}-subset.{ext} ',
      '-t sample_set -e {set}' )
    system( attr.set )
  }
}

process_gmts <- function( set ){
  ## set links for GMT files in the table
  get.gmt    <- glue( 'gsutil ls gs://{opt$bucket}/sample_sets/{set}/*.gmt' )
  gmt.path   <- system( get.gmt, intern = T )
  if ( length( gmt.path ) == 0 ) return()
  attr.set   <- glue( 'fissfc -V -y attr_set -w {opt$wk.space} ',
                      '-p {opt$project} ',
                      '-a geneset_db -v {gmt.path} ',
                      '-t sample_set -e {set}' )
  system( attr.set )
}

sample_set_attr <- function( sample.set.member.file )
{
  sets <- read_sets( sample.set.member.file )
  for ( set in names( sets ) )
  {
    process_other_attributes( set )
    process_additional_parameters( set )
    process_types( set )
    # process_gmts( set )

  }
}

main <- function()
{
  set_arguments()
  opt$csv.types <<- unlist( strsplit( opt$csv.types, split = ';' ) )
  opt$gct.types <<- unlist( strsplit( opt$gct.types, split = ';' ) )
  sample_set_attr( sample.set.member.file )
}

if ( !interactive() )
{
  main()
}
