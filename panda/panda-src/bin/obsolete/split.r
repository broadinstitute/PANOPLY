#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( cmapR )
p_load( optparse )
p_load( glue )

opt <- list()
input.gct.file <- ""

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-d", "--data.type" ),
      dest   = "data.type",
      action = "store",
      type   = 'character' ),

    make_option(
      c( "-e", "--ext"  ),
      dest   = "ext",
      action = "store",
      type   = 'character' ),

    make_option(
      c( "-i", "--input"  ),
      dest   = "input",
      action = "store",
      type   = 'character' ),

    make_option(
      c( "-o", "--output.location"  ),
      dest   = "output.location",
      action = "store",
      type   = 'character' ) )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

split_gct <- function()
{
  data <- parse.gctx( opt$input )
  invisible( sapply( 1:length( data@cid ), function( x ){
    cid <- c( data@cid[x] )
    mat <- as.matrix( data@mat[, x] )
    colnames( mat ) <- cid
    file.name <- glue( "{opt$output.location}/{cid}-{opt$data.type}.gct" )
    if ( nrow( data@cdesc ) == 0 )
      cdesc <- data.frame() else  cdesc <- data@cdesc[x, ]
    if ( nrow( data@rdesc ) == 0 )
      rdesc <- data.frame() else  rdesc <- data@rdesc
    rdesc[] <- lapply( rdesc, as.character )
    cdesc[] <- lapply( cdesc, as.character )
    gctclass <- new(
      "GCT",
      mat = mat,
      cdesc = cdesc,
      rdesc = rdesc,
      rid = data@rid,
      cid = cid,
      src = file.name )
    write.gct( gctclass, file.name, ver = 3, appenddim = FALSE )
  } ) )
}

split_csv <- function()
{
  data <- read.csv( opt$input, header = T, stringsAsFactors = F )
  rownames( data ) <- data$Sample.ID
  invisible( sapply( 1:nrow( data ), function( x ){
    scsv <- as.data.frame( data[ x, ], stringsAsFactors = F )
    file.name <- glue(
      "{opt$output.location}/{data$Sample.ID[x]}-{opt$data.type}.csv" )
    write.csv( scsv, file = file.name, row.names = F, quote = T )
  } ) )
}

main <- function(){
  set_arguments()
  if ( opt$ext == "csv" ) split_csv()
  else if ( opt$ext == "gct" ) split_gct()
}

if ( !interactive() ){
  main()
}
