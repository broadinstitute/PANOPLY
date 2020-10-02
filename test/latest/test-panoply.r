#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
## Script to check your test input tarball against a gold-standard tarball
## Inputs:
## gold.tar  : gold standard tarball
## test.tar  : test tarball
## tolerance : threhsold for comparing numerical matrices,
##             by default set to 1.490116e-08 or .Machine$double.eps ^ 0.5
## file.list : a YAML file listing patterns
## size.diff : a float number indicating a percent size change treshold for comparing image files 
##
## More on the YAML file:
##  1. keys should be sub-directories inside the tarball
##  2. each sub-directory should have a non-nested list of patterns
##  3. patterns should be regular expressions (for file matching w/ list.files())
##  4. patterns can also be of the form "!<pattern>" (ex. "!sample-info.csv or
##     "!*.tsv" )
##  5. By default this script will compare all CSV files assuming that
##     first row will be its header
##  6. By default thos script will compare all GCT file "mat" attributes

if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( cmapR )
p_load( glue )
p_load( optparse )
p_load( yaml )

opt <- list()

set_arguments <- function() {
  optionList <- list(
    make_option(
      c( "-g", "--gold.tar" ),
      dest = "gold.tar",
      action = "store",
      type = 'character'),

    make_option(
      c( "-t", "--test.tar" ),
      dest = "test.tar",
      action = "store",
      type = 'character'),

    make_option(
      c( "-f", "--file.list" ),
      dest = "file.list",
      action = "store",
      type = 'character'),

    make_option(
      c( "-d", "--tolerance" ),
      dest = "tolerance",
      action = "store",
      type = 'double'),
    
    make_option(
      c( "-s", "--size.diff" ),
      dest = "size.diff",
      action = "store",
      type = 'double') )
  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

## write out all reports to a test-report file
report <- function(){
  log.lines <- c( pass.lines, miss.lines, fail.lines )
  log.lines <- paste0( log.lines, collapse = '\n\n' )
  log.conn <- file( "test-report.txt")
  writeLines( log.lines, log.conn )
  close( log.conn )
}

## collect pass / fail checks for gct files
compare_gct <- function( gct.file ){
  gold.gct <- parse.gctx( glue( "{gold.dir}/{gct.file}" ) )
  test.gct <- parse.gctx( glue( "{test.dir}/{gct.file}" ) )
  check <- all.equal( gold.gct@mat, test.gct@mat, tolerance = tolerance )
  if( check == T ) pass.lines <<- c( pass.lines, glue( "++ pass. {gct.file}\n" ) )
  else fail.lines <<- c( fail.lines, glue( "XX fail. {gct.file}:\n\t{check}\n" ) )
}

## collect pass / fail checks for csv files
compare_csv <- function( csv.file ){
  gold.csv <- read.csv( glue( "{gold.dir}/{csv.file}" ),
                        header = T, stringsAsFactors = F )
  test.csv <- read.csv( glue( "{test.dir}/{csv.file}" ),
                        header = T, stringsAsFactors = F )
  check <- all.equal( gold.csv, test.csv, tolerance = tolerance )
  if( check == T ) pass.lines <<- c( pass.lines, glue( "++ pass. {csv.file}\n" ) )
  else fail.lines <<- c( fail.lines, glue( "XX fail. {csv.file}:\n\t{check}\n" ) )
}

## collect missing reports for data files
compare_files <- function( files ){
   for ( file in files ){
    print( glue( "Testing {file}..." ) )
    if ( !( file.exists( glue( "{test.dir}/{file}" ) ) ) ){
      miss.lines <<- c( miss.lines, glue( "-- missing. {file}\n" ) )
      next
    }
    ext <- tail( unlist( strsplit( file, split = '\\.' ) ), n = 1 )
    if ( ext == "gct" )
      compare_gct( file )
    else compare_csv( file )
   }
}

## parse patterns in your YAML file and get a list of files that need to be
## compared / tested
get_compare_file_list <- function( patterns, sub.dir ){
  files <- c()
  ignore <- c()
  for ( pattern in patterns ){
    ## if it is a regular non-ignore pattern
    if ( !( substring( pattern, 1, 1 ) == '!' ) ){
      files <- c( files, list.files( glue( "{gold.dir}/{sub.dir}" ),
                                     pattern = pattern, full.names = F ) )
    } else ignore <- c( ignore, substring( pattern, 2, nchar( pattern ) ) )
  }
  for ( ig in ignore )
    files <- grep( x = files, pattern = ig, invert = T, value = T )
  files <- unlist( lapply( files, function( x ) glue( "{sub.dir}/{x}" ) ) )
  return( files )
}


## compare data matrices
check_delta <- function(){
  for ( sub.dir in gold.dir.list ){
    if ( make.names( sub.dir ) %in% names( test.rules ) )
      patterns <- test.rules[[make.names( sub.dir )]] else
        patterns <- c( "*.gct", "*.csv", "*.txt" )
    files <- get_compare_file_list( patterns, sub.dir )
    compare_files( files )
  }
}

## collect missing reports for any images / plots in all sub directories
## collect unacceptable file size differences for any images / plots in all sub directories
## collect missing reports for a sub diectory itself
check_missing <- function(){
  for ( sub.dir in gold.dir.list ) {
    ## check sub directory
    if ( !( dir.exists( glue( "{test.dir}/{sub.dir}" ) ) ) ){
      miss.lines <<- c( miss.lines, glue( "-- missing. {sub.dir}\n" ) )
      next
    }
    
    ## check for plots and images
    gold.plot.list <- list.files( glue( "{gold.dir}/{sub.dir}/" ),
                                  pattern = "*\\.pdf$|\\.png$|\\.jpeg$",
                                  full.names = F, recursive = T )
    check.plot.list <- unlist( lapply( gold.plot.list, function( x )
      glue( "{sub.dir}/{x}" ) ) )
    
    for ( plot.file in check.plot.list ) {
      test.size <- file.size( glue( "{test.dir}/{plot.file}" ) )
      gold.size <- file.size( glue( "{gold.dir}/{plot.file}" ) )
      
      if ( !( file.exists( glue( "{test.dir}/{plot.file}" ) ) ) )
        miss.lines <<- c( miss.lines, glue( "-- missing. {plot.file}\n" ) )
      else if ( test.size/gold.size < ( 1 - diff.allow ) )
        miss.lines <<- c( miss.lines, glue( "-- Warning: file size too small. {plot.file}\n" ) )
    }
  }
}

untar <- function( tarball, dir ){
  if( dir.exists( dir ) )
    unlink( dir, recursive = T )
  dir.create( dir )
  untar.comm <- glue( "tar -xf {tarball} -C {dir} --strip-components 1" )
  system( untar.comm )
}

## untar input files and set up global variables
test_setup <- function(){
  if ( !( "tolerance" %in% names( opt ) ) )
    tolerance <<- .Machine$double.eps ^ 0.5 else
      tolerance <<- opt$tolerance
  gold.dir <<- "gold-dir"
  test.dir <<- "test-dir"
  untar( opt$gold.tar, gold.dir )
  untar( opt$test.tar, test.dir )
  if ( !( "file.list" %in% names( opt ) ) )
    test.rules <<- list() else
      test.rules <<- read_yaml( opt$file.list )
  gold.dir.list <<- list.dirs( glue( "{gold.dir}/" ),
                                full.names = F, recursive = F )
  if ( !( "size.diff" %in% names( opt ) ) )
    diff.allow <<- 0.5 else
      diff.allow <<- opt$size.diff
    
  pass.lines <<- c()
  miss.lines <<- c()
  fail.lines <<- c()
}

main <- function(){
  set_arguments()
  test_setup()
  check_missing()
  check_delta()
  report()
}

if ( !interactive() ){
  main()
}
