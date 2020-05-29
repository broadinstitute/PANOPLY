if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( cmapR )
p_load( glue )

opt <- list()

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-o", "--location"  ),
      dest   = "location",
      action = "store",
      type   = 'character'),

    make_option(
      c( "-w", "--wkspace" ),
      dest = "wkspace",
      action = "store",
      type = 'character' ),

    make_option(
      c( "-c", "--cols" ),
      dest = "cols",
      action = "store",
      type = 'character') )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

write_groups <- function(){
  ## annotation file
  annot.file <- glue( "{opt$location}/pipeline-input" )
  annot.file <- glue( "{annot.file}/{opt$wkspace}-annotation.csv" )
  annot.data <- read.csv( annot.file, header = T, stringsAsFactors = F )
  keep <- c( 'Sample.ID', unlist( strsplit( opt$cols, split = ';' ) ) )
  keep <- unlist( lapply( keep,
    function( x ) which( colnames( annot.data ) == x ) ) )
  groups <- annot.data[, keep]
  write.csv( groups,
             glue( "{opt$location}/pipeline-input/{opt$wkspace}-groups.csv" ),
             row.names = F,
             quote = T )
}

set_arguments()
write_groups()
