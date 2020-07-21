if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( glue )

sample.entity <- "samples.tsv"
opt <- list()
acc.patterns <- '\\.csv|\\.gct'

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-b", "--bucket" ),
      dest = "bucket",
      action = "store",
      type = 'character' ),

    make_option(
      c( "-w", "--wk.space" ),
      dest = "wk.space",
      action = "store",
      type = 'character'),

    make_option(
      c( "-p", "--project" ),
      dest = "project",
      action = "store",
      type = 'character'),

    make_option(
      c( "-c", "--csv.types" ),
      dest = "csv.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-g", "--gct.types" ),
      dest = "gct.types",
      action = "store",
      type = 'character'),

    make_option(
      c( "-o", "--location"  ),
      dest   = "location",
      action = "store",
      type   = 'character') )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

samples_attr <- function()
{
  sample.entity <- glue( "{opt$location}/{sample.entity}" )
  samples <- read.delim(
    sample.entity, header = T, sep = '\t', stringsAsFactors = F  )
  samples <- samples[['entity.sample_id']]
  types <- c( opt$csv.types, opt$gct.types )
  print( glue(
    '{length( samples )} samples present ',
    'in firecloud data model.' ) )
  for ( type in types ){
    ls.commd <- glue( 'gsutil ls gs://{opt$bucket}/{type}/' )
    ls.samples <- system( ls.commd, intern = T )
    ls.samples <- unlist( lapply(
      ls.samples,
      function( x ) grep(
        tail( unlist( strsplit( x, split = '/' ) ), 1 ),
        pattern = acc.patterns, value = T ) ) )
    ls.samples <- unlist( lapply(
      ls.samples,
      function( x ) gsub(
        '\\.', '-',
        head( unlist( strsplit( head( unlist(
          strsplit( x, split = acc.patterns ) ), 1 ),
          split = glue( '-{type}' ) ) ), 1 ) ) ) )
    ed.samples <- intersect( samples, ls.samples )
    print( glue( '{length( ed.samples )} files available ',
                 'for {type} in google bucket.' ) )
    for ( sample in ed.samples ){
      if ( type %in% opt$gct.types )  ext <- 'gct' else ext <- 'csv'
      print(glue("{sample}, {make.names(sample)}"))
      attr.set <- glue(
        'fissfc -V -y attr_set -w {opt$wk.space} ',
        '-p {opt$project} ',
        '-a {type} ',
        '-v gs://{opt$bucket}/{type}/{make.names( sample )}-{type}.{ext} ',
        '-t sample -e {sample}' )
      system( attr.set )
    }
  }
}

main <- function(){
  set_arguments()
  opt$csv.types <<- unlist( strsplit( opt$csv.types, split = ';' ) )
  opt$gct.types <<- unlist( strsplit( opt$gct.types, split = ';' ) )
  samples_attr()
}

if ( !interactive() ) {
  main()
}
