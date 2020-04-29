if( !suppressMessages( require( "pacman" ) ) )
  install.packages( "pacman" )
p_load( optparse )
p_load( glue )

opt <- list()

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-p", "--panoply" ),
      dest = "panoply",
      action = "store",
      type = 'character' ),
    make_option(
      c( "-t", "--task" ),
      dest = "task",
      action = "store",
      type = 'character' ) )
  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

create_map <- function(){
  task.path <- glue( "{opt$panoply}/hydrant/tasks" )
  task.list <- grep( 'panoply*',
                     list.dirs( task.path, full.names = F, recursive = F ),
                     value = TRUE)
  map <- list()
  for ( task in task.list ){
    dockerfile <- glue( "{task.path}/{task}/{task}/Dockerfile" )
    if ( !file.exists( dockerfile ) ) next
    parent_task <- make.names( unlist(
      strsplit( readLines( dockerfile, n = 1 ), split = '/' ) )[2] )
    parent_task <- unlist( strsplit( parent_task, split = '\\.' ) )[1]
    if ( parent_task == "NA" ){
      parent_task <- make.names( unlist( strsplit(
        readLines( dockerfile, n = 1 ), split = ' ' ) )[2] )
    }
    if ( !( parent_task %in% names( map ) ) )
      map[[parent_task]] <- c( task )
    else map[[parent_task]] <- c( map[[parent_task]], task )
  }
  map.file <- glue( "{task.path}/map.txt" )
  blurb <- c()
  for ( parent in names( map ) )
    blurb <- c( blurb,
                glue( "={parent}\n;{paste0( map[[parent]],
                      collapse = '\n;' )}" ) )
  writeLines( blurb, map.file, sep = '\n')
  return( map )
}

recurse_targets <- function( node, map, targets ) {
  if ( !( node %in% names( map ) ) ) return( c() )
  children <- map[[node]]
  #print( glue( "{node}->{paste0(children, collapse=';')}" ) )
  for ( child in children ) {
    targets <- c( targets, child, recurse_targets( child, map, c() ) )
    #print( glue( "{child} || {paste0( targets, collapse='; ' )}" ) )
  }
  return( targets )
}

get_targets <- function( map ) {
  targets <- recurse_targets( opt$task, map, c() )
  dtarget <- glue( "{opt$panoply}/hydrant/tasks/targets" )
  dir.create( dtarget )
  ftarget <- glue( "{dtarget}/{opt$task}-targets.txt" )
  writeLines( paste0( targets, collapse = ';' ), ftarget )
}

main <- function(){
  set_arguments()
  opt$task <<- make.names( opt$task )
  map <- create_map()
  get_targets( map )
}

if ( !interactive() ) {
  main()
}
