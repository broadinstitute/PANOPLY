#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( cmapR )
p_load( glue )

opt <- list()

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-i", "--file.id" ),
      dest = "file.id",
      action = "store",
      type = 'character'),

    make_option(
      c( "-c", "--cols" ),
      dest = "cols",
      action = "store",
      type = 'character') )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

write_groups <- function(){
  pome <- parse.gctx( glue( "{opt$file.id}-proteome.gct" ) )
  keep <- c( 'Sample.ID', unlist( strsplit( opt$cols, split = ';' ) ) )
  keep <- unlist( lapply( keep,
    function( x ) which( colnames( pome@cdesc ) == x ) ) )
  groups <- pome@cdesc[, keep]
  write.csv( groups, glue( "{opt$file.id}-groups.csv" ), row.names = F,
             quote = T )
}

set_arguments()
write_groups()
