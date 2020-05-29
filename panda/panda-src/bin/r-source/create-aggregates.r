if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( cmapR )
p_load( dplyr )
p_load( glue )

opt <- list()

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-a", "--csv.types"  ),
      dest = "csv.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-b", "--gct.types"  ),
      dest = "gct.types",
      action = "store",
      type = 'character' ),

    make_option(
      c( "-o", "--location"  ),
      dest   = "location",
      action = "store",
      type   = 'character'),

    make_option(
      c( "-w", "--wkspace" ),
      dest = "wkspace",
      action = "store",
      type = 'character' )

    )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
  opt$meta.file <<- glue( "{opt$location}/sample_set_membership.tsv" )
}

read_sets <- function(){
  meta <- read.delim( opt$meta.file, header = T, sep = '\t',
                      stringsAsFactors = F )
  sets <- list()
  for ( rIdx in 1:nrow( meta ) ){
    set_id <- meta$membership.sample_set_id[rIdx]
    if ( set_id %in% names( sets ) )
      sets[[set_id]] <- make.names( c( sets[[set_id]], meta$sample[rIdx] ) )
    else sets[[set_id]] <- make.names( c( meta$sample[rIdx] ) )
  }
  return( sets )
}

combine_gcts <- function( file.list, type, set ){
  mat <- data.frame()
  cdesc <- data.frame()
  for ( file.name in file.list ){
    if ( !file.exists( file.name ) ){
      print( glue( "{file.name} does not exist" ) )
      next
    }
    gct.data <- parse.gctx( file.name )
    if ( identical ( mat, data.frame() ) ){
      mat   <- gct.data@mat
      cdesc <- gct.data@cdesc
    } else {
      mat   <- cbind( mat, gct.data@mat )
      cdesc <- rbind( cdesc, gct.data@cdesc )
    }
  }
  if ( exists( "gct.data" ) ){
    rdesc <- gct.data@rdesc
    aggrt <- new( "GCT", mat = as.matrix( mat ), cdesc = cdesc, rdesc = rdesc,
                  rid = rownames( mat ), cid = colnames( mat ), src = out.file )
    set.path <- glue( "{opt$location}/aggregates/{set}" )
    dir.create( set.path, showWarnings = T )
    out.file <- glue( "{set.path}/{type}-aggregate.gct" )
    write.gct( aggrt, out.file, ver = 3, appenddim = FALSE )

    if ( type == 'rna' ){
      out.file <- glue( "{set.path}/{type}-v2-aggregate.gct" )
      write.gct( aggrt, out.file, ver = 2, appenddim = FALSE )
    }
  }
}

combine_csv <- function( file.list.all, type, set ){
  edat <- data.frame()
  for ( file.name in file.list.all ){
    if ( identical( edat, data.frame() ) ){
      edat <- read.csv( file.name, header = T, stringsAsFactors = F )
    } else edat <- rbind( edat, read.csv( file.name, header = T,
                                          stringsAsFactors = F ) )
  }
  set.path <- glue( "{opt$location}/aggregates/{set}" )
  out.file <- glue( "{set.path}/{type}-aggregate.csv" )
  write.csv( edat, out.file, row.names = F, quote = T, na = "" )
}

download_files <- function(){
  addr.file <- glue( "{opt$location}/filled_samples.tsv" )
  addr <- read.csv( addr.file, header = T, stringsAsFactors = F, sep = '\t' )
  types <- c( opt$gct.types, opt$csv.types, "groups" )
  invisible( lapply( types, function( t ){
    dir.create( glue( "{opt$location}/split-data/{t}" ) )
    links <- addr[[t]]
    invisible( lapply( links, function( link ){
      system( glue( "gsutil cp {link} {opt$location}/split-data/{t}/." ) )
    } ) )
  } ) )
}

aggregate_all <- function(){
  sets <- read_sets()
  download_files()
  for ( set in names( sets ) ){

    file.list.all <- sapply(
      opt$gct.types, function( x ){
        unlist( lapply( sets[[set]], function( y )
          glue( "{opt$location}/split-data/{x}/{y}-{x}.gct" ) ) ) } )

    invisible( lapply( opt$gct.types, function( x )
      combine_gcts( file.list.all[, x], x, set ) ) )

    file.list.all <- sapply(
      opt$csv.types, function( x ){
        unlist( lapply( sets[[set]], function( y )
          glue( "{opt$location}/split-data/{x}/{y}-{x}.csv" ) ) ) } )

    invisible( lapply( opt$csv.types, function( x )
      combine_csv( file.list.all[, x], x, set ) ) )
  }
}

main <- function(){
  set_arguments()
  opt$csv.types <<- unlist( strsplit( opt$csv.types, split = ';' ) )
  opt$gct.types <<- unlist( strsplit( opt$gct.types, split = ';' ) )
  aggregate_all()
}

if ( !interactive() ){
  main()
}
