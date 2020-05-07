if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( cmapR )
p_load( glue )

opt <- list()

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-l", "--gct.file.list"  ),
      dest = "gct.file.list",
      action = "store",
      type = 'character'),

    make_option(
      c( "-g", "--gct.types"  ),
      dest = "gct.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-c", "--csv.types"  ),
      action = "store",
      dest = "csv.types",
      type = 'character'),

    make_option(
      c( "-o", "--output.location"  ),
      dest   = "output.location",
      action = "store",
      type   = 'character') )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
  descriptors <<- c( unlist( strsplit( opt$csv.types, split = ';' ) ),
                     unlist( strsplit( opt$gct.types, split = ';' ) ) )
  opt$gct.file <<- unlist( strsplit( opt$gct.file.list, split = ';' ) )[1]
}

get_ids <- function(){
  data <- parse.gctx( opt$gct.file )
  id_table <- data@cdesc[, c( 'Sample.ID', 'Participant' )]
  colnames( id_table ) <- c( 'Sample.ID', 'Participant.ID' )
  id_table$Sample.ID <- gsub( '\\.', '-', id_table$Sample.ID )
  id_table$Participant.ID <- gsub( '\\.', '-', id_table$Participant.ID )
  return( id_table )
}

get_normal_ids <- function(){
  return( which( data@cdesc[,'Type'] == 'Normal' ) )
}

populate_participant <- function(){
  id_table <- as.data.frame(
    unique( get_ids()$Participant.ID ), stringsAsFactors = F )
  colnames( id_table ) <- c( 'entity:participant_id' )
  write.table( id_table, glue( "{opt$output.location}/participant.tsv" ),
               quote = F, sep = '\t', na = 'NA', row.names = F, col.names = T  )
}

populate_sample <- function()
{
  # Assuming opt$gct.file does not change the order of its sample IDs
  sample_meta <- cbind( get_ids(), parse.gctx( opt$gct.file )@cdesc[, 'Type'] )
  fill_nas <- rep( NA, nrow( sample_meta ) )
  colnames( sample_meta ) <- c( 'entity:sample_id', 'participant', 'type' )
  for ( x in 1:length( descriptors ) )
    sample_meta <- cbind( sample_meta, fill_nas )
  colnames( sample_meta )[-c( 1:3 )] <- descriptors
  write.table( sample_meta, glue( "{opt$output.location}/samples.tsv" ),
               quote = F, sep = '\t', na = 'NA', row.names = F, col.names = T  )
}

main <- function(){
  set_arguments()
  populate_participant()
  populate_sample()
}

if ( !interactive() ){
  main()
}
