if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( cmapR )
p_load( glue )
opt <- list()

set_arguments <- function() {
  optionList <- list(
    make_option(
      c( "-g", "--gct.types" ),
      dest = "gct.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-c", "--csv.types" ),
      dest = "csv.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-s", "--set.dir" ),
      dest = "set.dir",
      action = "store",
      type = 'character') )
  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

set_descriptors <- function(){
  opt$suffix <<- "ss"
  opt$gct.types <- unlist( strsplit( opt$gct.types, split = ';' ) )
  opt$csv.types <- unlist( strsplit( opt$csv.types, split = ';' ) )
  descriptors <<- c( opt$csv.types, opt$gct.types )
  output_links <- c( "complete_runtime_results",
                     "summary_results",
                     "rna_corr_report",
                     "cna_corr_report",
                     "sample_qc_report" )
  temp_desc <- descriptors
  for ( type in opt$gct.types ){
    if ( type %in% c( 'cna', 'rna' ) ) next
    temp_desc <- c(
      temp_desc,
      unlist( lapply( output_links, function( x ) glue( "{type}_{x}" ) ) ),
      glue( "{type}_add_params" ) )
  }
  descriptors   <<- temp_desc
  rm( temp_desc )
}

populate_sample_sets <- function(){
  membership <- data.frame()
  set.entity <- data.frame()

  set.files <- list.files( opt$set.dir, pattern = "*.txt" )

  for ( set.file.name in set.files ){
    set_name <- unlist( strsplit( set.file.name, split = '.txt' ) )[1]
    set.members <- readLines( glue( "{opt$set.dir}/{set.file.name}" ), n = -1 )
    set.members <- unlist(
      lapply( set.members, function( x ) gsub( '\\.', '-', x ) ) )
    if( identical( membership, data.frame() ) ){
      membership <- cbind(
        rep( set_name, length( set.members ) ),
        set.members )
      set.entity <- cbind.data.frame(
        split(
          c( set_name, rep( 'NA', length( descriptors ) ) ),
          rep( 1:( length( descriptors ) + 1 ) ) ),
        stringsAsFactors = F )
    } else {
      membership <- rbind(
        membership,
        cbind( rep( set_name, length( set.members ) ), set.members ) )
      set.entity <- rbind(
        set.entity,
        cbind.data.frame(
          split(
            c( set_name, rep( 'NA', length( descriptors ) ) ),
            rep( 1:(length(descriptors)+1) ) ),
          stringsAsFactors = F ) )
    }
  }

  colnames( membership ) <- c( 'membership:sample_set_id', 'sample_id' )
  colnames( set.entity ) <- c( 'entity:sample_set_id',
                               paste0( descriptors, '_', opt$suffix ) )

  write.table( membership, 'sample_set_membership.tsv', quote = F, sep = '\t',
               na = 'NA', row.names = F, col.names = T, append = F  )
  write.table( set.entity, 'sample_set_entity.tsv', quote = F, sep = '\t',
               na = 'NA', row.names = F, col.names = T, append = F  )
}

main <- function(){
  set_arguments()
  set_descriptors()
  populate_sample_sets()
}

if ( !interactive() ){
  main()
}
