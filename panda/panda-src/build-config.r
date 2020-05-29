## Loading libraries needed to run this notebook
library( cmapR );
library( pacman );
library( RColorBrewer )
p_load( scales );
p_load( glue );
p_load( yaml );

### ====
### Section 0. Setting global parameters
### ====

### separators
sep <- "============================"
DONE <- glue( "\n{sep}\n.. DONE.\n{sep}\n" )
WAIT <- glue( "\n{sep}\n.. Please wait..\n{sep}\n" )

cat.map <- c(
  "proteome",
  "phosphoproteome",
  "acetylome",
  "rna",
  "cna",
  "annotation",
  "groups" )

required.cols <- c(
  "Sample.ID",
  "Participant",
  "Experiment",
  "Channel",
  "Type" )

### "globals" variable
globals <- list()
globals$project <- "broad-firecloud-cptac"
globals$group <- "GROUP_Broad_CPTAC@firecloud.org"
globals$meth_space <- "broadcptac"

### small global functions
add_space <- function( arr ){
  for ( id in 1:length( arr ) )
    if ( id < 10 )
      arr[id] <- glue( " {as.character( arr[id] )}" )
   return( arr )
}

trim <- function( x )
  gsub( "^\\s+|\\s+$", "", x )

y2true <- function( choice ){
  if ( choice == "y" )
    flag <- T else flag <- F
    return( flag )
}

### ====
### Section. Inputs
### ====

map_my_files <- function( ext ){
  mapped <- list()
  for ( file in ext ){
    category <- as.numeric( readline( prompt = glue( ".. {file}: " ) ) )
    if ( category == 0 ){
      cat.map <- c( cat.map, readline( prompt = "\n$$ Enter category name:" ) )
      category <- length( cat.map )
    } else if ( category < 0 || category > 7 )
      stop( "\n\n{sep}\n.. ERROR. Invalid entry. Try again.." )
    mapped[[cat.map[category]]] <- file
  }
  return ( mapped )
}

load_unzipped_files <- function(){
  home <- glue( "/home/jupyter-user/notebooks/{terra.wkspace}/edit/" )
  setwd( home )
  ## locate uploaded file
  input.zip.name <- readline(
    prompt = "\n$$ Enter uploaded zip file name (test.zip): " )
  input.zip <<- glue( "{google.bucket}/{input.zip.name}" )
  zip.name <- tail( unlist( strsplit( input.zip, split = '/' ) ), 1 )
  system( glue( "gsutil cp {input.zip} {getwd()}/." ) )
  if( dir.exists( "input" ) )
    unlink( "input", recursive = T )
  dir.create( "input" )
  system( glue( "unzip {zip.name} -d input/" ) )
  setwd( glue( "{getwd()}/input" ) )
  gcts <- list.files( pattern = "*.gct" )
  csvs <- list.files( pattern = "*.csv" )
  typemap.gct <<- map_my_files( gcts )
  typemap.csv <<- map_my_files( csvs )
  #return( list( gct = typemap.gct, csv = typemap.csv ) )
}

panda_input <- function(){
  load_unzipped_files()
  print( DONE )
}

### ====
### Section. Groups
### ====

get_all_groups <- function( typemap.csv ){
  annot <- read.csv( typemap.csv$annotation, header = T,
                     stringsAsFactors = F, quote = '"' )
  all.groups <- colnames( annot )
  missing.cols <- setdiff( required.cols, all.groups )
  error <- glue(
    "Missing Columns. \n{paste0( missing.cols, collapse = ', ' )}" )
  if( length( missing.cols ) > 0 )
    stop( glue( "\n\n{sep}\n.. ERROR. {error}.\n{sep}\n" ) )
  else
    print( glue( "\n\n{sep}\n.. INFO. No missing columns.\n{sep}\n" ) )
  return( all.groups )
}

display_all_groups <- function( typemap.csv ){
  all.groups <- get_all_groups( typemap.csv )
  print( glue( "\n\n{sep}\n.. Annotations:\n\n" ) )
  print( glue( "{add_space( 1:length( all.groups ) )}: {all.groups}" ) )
  print( glue( "\n{sep}\n" ) )
  return( all.groups )
}

select_groups_case <- function( all.groups, groups.cols ){
  groups.cols <- c()
  stmt <- glue( "\n$$ Enter an index or a range of indices " )
  stmt <- glue( "{stmt}(5 or 23:34): " )
  range <- unlist( strsplit(
      readline( prompt = stmt ), split = ':' ) )
  lwend <- as.numeric( range[1] )
  if( length( range ) == 1 )
    upend <- lwend else
      upend <- as.numeric( range[2] )
    groups.cols <- c( groups.cols, all.groups[lwend:upend] )
  return( groups.cols )
}

select_groups <- function( all.groups ){
  groups.flag <- T
  groups.cols <- c()
  while ( groups.flag ){
    groups.cols <- c( groups.cols, select_groups_case(
      all.groups, groups.cols ) )
    groups.flag <- y2true( readline(
      prompt = "$$ Continue adding groups? (y/n): "))
  }
  return( groups.cols )
}

verify_group_validity <- function( groups.cols, typemap.csv ){
  drop <- c()
  annot <- read.csv( typemap.csv$annotation, header = T,
                     stringsAsFactors = F, quote = '"' )
  for ( group.idx in 1:length( groups.cols ) ){
    groups.vals <- unique( annot[[groups.cols[group.idx]]] )
    if ( length( groups.vals ) > 5 || length( groups.vals ) == 1 ){
      drop <- c( drop, group.idx )
    }
  }
  if ( length( drop ) > 0 ) groups.cols <- groups.cols[-drop]
  return( groups.cols )
}


display_validated_groups <- function( groups.cols ){
  warning <- "Annotations with <2 or 10+ categories were not added."
  print( glue( "\n\n{sep}\n.. WARNING. {warning}" ) )
  print( glue( "\n.. Your groups:\n\n" ) )
  print( glue( "{add_space( 1:length( groups.cols ) )}: {groups.cols}" ) )
}

panda_display_groups <- function(){
  all.groups  <<- display_all_groups( typemap.csv )
}
panda_select_groups <- function(){
  print( WAIT )
  groups.cols <<- select_groups( all.groups )
  print( WAIT )
  groups.cols <<- verify_group_validity( groups.cols, typemap.csv )
  print( WAIT )
  display_validated_groups( groups.cols )
  print( DONE )
}

### ====
### Section. Colors
### ====

display_default_colors <- function( groups.cols, groups.colors ){
  par( cex.axis = 1.2 )
  par( mar = c( 0, 0, 0, 0 ) )
  par( pin = c( 3, 0.2 ) )
  par( las = 3 )
  par( mfrow = c( 2, 2 ) )
  for ( group in names( groups.colors ) ){
    values <- groups.colors[[group]]$vals
    colors <- groups.colors[[group]]$colors
    values[is.na( values )] <- "NA"
    group.pal <- colors
    names( group.pal ) <- values
    image( 1:length( colors ), 1, as.matrix( 1:length( colors ) ),
           col = colors,
           xlab = "",
           ylab = "",
           xaxt = "n",
           yaxt = "n",
           bty = "n" )
    title( group, line = -5 )
    axis( 3, at = 1:length( colors ), labels = values, tick = F )
  }
}

set_default_colors <- function( groups.cols,typemap.csv ){
  annot <- read.csv( typemap.csv$annotation, header = T,
                     stringsAsFactors = F, quote = '"' )
  pair.colors <- brewer.pal( n = 12, name = "Paired" )
  qual.pals <- c( "Set1", "Dark2", "Set2", "Set3" )
  qual.maxc <- c( 9, 8, 8, 12 )
  qual.colors <- c()
  for ( pal.idx in 1:length( qual.pals ) ) {
    qual.colors <- c( qual.colors,
                      brewer.pal( n = qual.maxc[pal.idx],
                                  name = qual.pals[pal.idx] ) )
  }
  pair.count <- 0
  qual.count <- 1
  groups.colors <- list()
  for ( group in groups.cols  ) {
    groups.vals <- unique( annot[[group]] )
    groups.vals[is.na( groups.vals )] <- "NA"
    if( length( groups.vals ) == 2 || ( length( groups.vals ) == 3
                                        && "NA" %in% groups.vals ) ){
      pair.idx <- ( pair.count * 2 ) + 1
      colors <- pair.colors[ c( pair.idx, pair.idx + 1 ) ]
      pair.count <- ( pair.count + 1 ) %% 6
      if ( length( groups.vals ) == 3 ){
        groups.vals <- c( groups.vals[-which( groups.vals == "NA" )], "NA" )
        colors <- c( colors, "#EEEEEE" )
      }
    } else{
      upend <- qual.count + length( groups.vals ) - 1
      colors <- qual.colors[qual.count:upend]
      qual.count <- ( upend %% length( qual.colors ) ) + 1
    }
    groups.colors[[group]]$vals <- groups.vals
    groups.colors[[group]]$colors <- colors
  }
  return( groups.colors )
}

panda_defaults <- function(){
  groups.colors <<- set_default_colors( groups.cols, typemap.csv )
  print( WAIT )
  display_default_colors( groups.cols, groups.colors )
  print( DONE )
}

panda_colors_edit_reference <- function( groups.colors ){
  for ( group.idx in 1:length( names( groups.colors ) ) ){
    group <- names( groups.colors )[group.idx]
    print( glue( "{sep}\n{add_space( group.idx )}] {group} =\n" ) )
    values <- groups.colors[[group]]$vals
    print( glue( "\t{add_space( 1:length( values ) )}: {values}" ) )
  }
  print( DONE )
}

change_default_colors <- function( groups.cols, all.groups, groups.colors ){
  groups.flag <- TRUE
  while ( groups.flag ){
    stmt <- glue( "\n$$ Enter group index for color modification: " )
    change.group <- groups.cols[as.numeric( readline( prompt = stmt ) )]
    ## iterate over values for group
    value.flag = T
    while( value.flag ){
      prompt.stmt <- glue( "\n$$ Enter value index: " )
      change.value <- as.numeric( readline( prompt = prompt.stmt ) )
      change.color <- trim( as.character(
        readline( prompt = "\n$$ Enter color hex: " ) ) )
      if ( substring( change.color, 1, 1 ) != '#'
           || nchar( change.color ) != 7 )
        stop( glue( "\n\n{sep}\n.. ERROR. Invalid color. Try Again.." ) )
      groups.colors[[change.group]]$colors[change.value] <- change.color
      value.flag <- y2true( readline( prompt = glue(
        "\n$$ Continue editing colors for {change.group}? (y/n): " ) ) )
    }
    groups.flag <- y2true(
      readline(
        prompt = "\n$$ Continue modifying defaults for other groups? (y/n): " ) )
  }
  display_default_colors( groups.cols, groups.colors )
  return( groups.colors )
}

panda_colors_edit <- function(){
  groups.colors <<- change_default_colors(
    groups.cols, all.groups, groups.colors )
  print( DONE )
}


### ===
### Section. Sample sets
### ===

sample_set_sanity_check <- function( sample.sets ){
  command <- glue( "fissfc sset_list -w {terra.wkspace} -p {globals$project}" )
  exist.ss <- system( command, intern = T )
  if ( "all" %in% exist.ss )
    sample.sets$all <- NULL
  if ( !( identical( character(), exist.ss ) ) ){
    overlap.flags <- names( sample.sets ) %in% exist.ss
    if ( any( overlap.flags ) ){
      overlap.ssets <- names( sample.sets )[overlap.flags]
      warning <- "WARNING. Some set names you entered already exist."
      print( glue( "\n\n{sep}\n.. {warning} Duplicates created." ) )
      names( sample.sets )[overlap.flags] <- paste0( overlap.ssets, "_copy" )
    }
  }
  return( sample.sets )
}

define_sample_sets <- function( all.groups, typemap.csv ){
  annot <- read.csv( typemap.csv$annotation, header = T,
                     stringsAsFactors = F, quote = '"' )
  sample.sets <- list()
  sample.sets$all$fil_col <- "Type"
  sample.sets$all$fil_val <- glue(
    "{paste0( unique( annot[[sample.sets$all$fil_col]] ), collapse = ';' )}" )
  set.flag <- y2true( readline( prompt = "\n$$ Add sample subsets? (y/n) : " ) )
  while( set.flag ){
    prompt.stmt <- glue( "{sep}\n$$ Enter set name: " )
    name <- as.character( readline( prompt = prompt.stmt ) )
    fil_col <- as.character( readline(
      prompt = "$$ Filter samples based on this column. Enter name: " ) )
    if ( !( fil_col %in% all.groups ) )
      stop( "\n\n{sep}\n.. ERROR. Invalid column entry. Try Again." )
    values <- unique( annot[[fil_col]] )

    select.value.msg <- "$$ Enter one or more values from above."
    prompt.stmt <- glue( "{paste0( values, collapse = ', ' )}" )
    prompt.stmt <- glue( ".. Values: {prompt.stmt}\n{select.value.msg}" )
    prompt.stmt <- glue( "{prompt.stmt} (Tumor or LumA;LumB): " )
    fil_val <- trim( as.character( readline( prompt = prompt.stmt ) ) )
    if ( !( fil_val %in% values ) )
      stop( glue( "\n\n{sep}\n.. ERROR. Value entered is incorrect." ) )
    sample.sets[[name]]$fil_col <- fil_col
    sample.sets[[name]]$fil_val <- fil_val

    add.param.flag <- y2true(
      readline( prompt = "$$ Add set-specific parameters? (y/n): " ) )
    while( add.param.flag ){
      param.name <- as.character( readline(
        prompt = "$$ Enter Parameter Name: " ) )
      param.value <- as.character( readline(
        prompt = "$$ Enter parameter value: " ) )
      sample.sets[[name]][[param.name]] <- param.value
      add.param.flag <- y2true( readline(
        prompt = "$$ Continue adding more set parameters? (y/n): " ) )
    }
    set.flag <- y2true(
      readline( prompt = "\n$$ Continue adding more sets? (y/n): " ) )
  }
  sample.sets <- sample_set_sanity_check( sample.sets )
  if ( length( names( sample.sets ) ) == 0 )
    print( glue( "\n{sep}\n.. No sample sets will be added." ) ) else{
      print(
        glue( "\n{sep}\n.. Sample sets to be added to Terra Workspace: " ) )
      print( glue( "\n.. {paste0( names( sample.sets ), collapse = ', ' )}" ) )
  }
  return( sample.sets )
}

panda_create_sets <- function(){
  sample.sets <<- define_sample_sets( all.groups, typemap.csv )
  print( DONE )
}

### ===
### Section. Build YAML
### ===

write_config_yaml <- function() {
  lines <- list()
  lines$TEDMINT <- "/panda/bin"
  lines$freeze.path <-
    glue( "/home/jupyter-user/notebooks/{terra.wkspace}/input/" )
  lines$wkspace <- terra.wkspace
  lines$globals <- globals
  lines$typemap.gct <- typemap.gct
  lines$typemap.csv <- typemap.csv
  lines$groups.cols <- groups.cols
  lines$groups.colors <- list()
  for ( group in names( groups.colors ) ) {
    lines$groups.colors[[group]] <- list()
    for ( idx in 1:length(groups.colors[[group]]$vals ) ) {
      value <- groups.colors[[group]]$vals[idx]
      color <- groups.colors[[group]]$colors[idx]
      lines$groups.colors[[group]][[value]] <- color
    }
  }
  lines$sample.sets <- sample.sets
  setwd( glue( "{getwd()}/../" ) )
  write_yaml( x = lines, file = "config.yaml" )
  system( glue( "gsutil cp {getwd()}/config.yaml {google.bucket}/config.yaml" ) )
  print( DONE )
}

### ===
### Section. PANDA tool
### ===

run_panda <- function(){
  runspace <- "runspace"
  if ( dir.exists( runspace ) )
    unlink( runspace, recursive = T )
  dir.create( runspace )
  print( WAIT )
  invisible( file.copy( "/panda/bin/Makefile", glue( "{runspace}/Makefile" ) ) )
  config.addr <- glue( "{google.bucket}/config.yaml" )
  system( glue( "gsutil cp {config.addr} {runspace}/." ) )
  setwd( "runspace" )
  system( glue( "make restart && make init" ) )
  setwd( "pipeline-input" )
  for ( type in names( typemap.gct ) ) {
    system( glue(
      "ln -s ../../input/{typemap.gct[[type]]} {terra.wkspace}-{type}.gct" ) )
  }
  for ( type in names( typemap.csv ) ) {
    system( glue(
      "ln -s ../../input/{typemap.csv[[type]]} {terra.wkspace}-{type}.csv" ) )
  }
  setwd( "../" )
  system( glue( "make groups" ) )
  system( glue( "make terrainit" ) )
  run.all.flag <- y2true( readline(
    prompt = "$$ (Re)Run PANDA on samples too? (y/n): " ) )
  if( run.all.flag )
    system( glue( "make panda-samples" ) )
  print( WAIT )
  system( glue( "make panda-sets" ) )
  print( DONE )
}


