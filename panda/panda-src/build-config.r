#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


##
## NOTE: To run this script manually (outside of the Terra Jupyter notebook)
##       set the following env variables (in bash)
##       % export PANDA=<$PANOPLY/panda/panda-src>
##       % export WORKSPACE_NAME=<Terra workspace>
##       % export WORKSPACE_BUCKET=<Bucket id for Terra workspace>
##       % export WORKSPACE_NAMESPACE=<Billing project for Terra workspace>
##       Souce this script from R, and then run R functions to replicate notebook
##       (ie. R code in notebook code blocks, in sequence)
##       Working directory will be ./$WORKSPACE_NAME
## Additional NOTE: Google-Cloud-SDK, R-3.5, Perl-5.8 and Python-2.7 need to be loaded
##


## Loading libraries needed to run this notebook
library( cmapR );
library( pacman );
library( RColorBrewer )
p_load( scales );
p_load( glue );
p_load( yaml );
p_load( ids );

### ====
### Section 0. Setting global parameters
### ====
options (warn=-1)

### separators
sep <- "============================"
DONE <- glue( "\n{sep}\n.. DONE.\n{sep}\n" )
WAIT <- glue( "\n{sep}\n.. Please wait ..\n{sep}\n" )

cat.map <- c(
  "proteome",
  "phosphoproteome",
  "acetylome",
  "ubiquitylome",
  "rna",
  "cna",
  "annotation",
  "groups",
  "parameters",
  "ptmseaDB",
  "gseaDB")

proteome.types <- c(
  "proteome",
  "phosphoproteome",
  "acetylome",
  "ubiquitylome"
)

required.genomics.types <- c(
  "cna",
  "rna"
)

required.cols <- c(
  "Sample.ID",
  "Type" )

ignore.cols <- c(
  "Experiment",
  "Channel",
  "Participant",
  "QC.status" )



### "global" variable
globals <- list()
globals$meth_space <- "broadcptac"
# globals$project defined in panda_initialize
# globals$group unnecessary, excluded (workspace already exists)

### Home dir on Terra
terra_notebook_home <- "/home/jupyter"

### defaults (files to be located in docker image, or locally)
# for local use outside of a notebook in Terra, set PANDA to point to $PANOPLY/panda/panda-src
# if $PANDA is not found, default to /panda as the code base
if ( (panda <- Sys.getenv("PANDA")) == "" ) panda <- "/panda"
Sys.setenv (PANDA=panda)   # set in env so that other scripts can find the code base
defaults <- list()
defaults$parameters <- glue ("{panda}/defaults/master-parameters.yaml")
defaults$ptmsea_db <- glue ("{panda}/defaults/ptm.sig.db.all.uniprot.human.v1.9.0.gmt")
defaults$gsea_db <- glue ("{panda}/defaults/h.all.v6.2.symbols.gmt")
defaults$all_parameters_file_name <- "panoply-parameters.yaml"
defaults$panda_parameters_file_name <- "config.yaml"

### small global functions
process_file <- function (f) {
  # copy f to current directory and return name of file
  name <- basename(f)
  file.copy(from=f, to=name)
  return (name)
}

read_annot <- function () {
  # reads annotation file, if specified; else throws an error
  if (is.null (typemap.csv$annotation)) {
    stop( glue( "\n\n{sep}\n.. ERROR. Sample annotation file missing. Run panda_input().\n{sep}\n" ) )
  }
  annot <- read.csv( glue ("{home}/input/{typemap.csv$annotation}"), header = T,
                     stringsAsFactors = F, quote = '"' )
  return (annot)
}

add_space <- function( arr ){
  for ( id in 1:length( arr ) )
    if ( id < 10 )
      arr[id] <- glue( " {as.character( arr[id] )}" )
   return( arr )
}

trim <- function( x )
  gsub( "^\\s+|\\s+$", "", x )

y2true <- function( prompt ) {
  # displays prompt and reads in y/n response
  valid_choice <- function (ch) ch %in% c ('Y', 'y', 'N', 'n')
  while (! valid_choice (choice <- readline (prompt = glue ("{prompt} (y/n): "))) ) {}
  if ( choice == "y" || choice == "Y")
    flag <- TRUE else flag <- FALSE
  return( flag )
}

printX <- function (type, line1, line2=NULL, line3=NULL) {
  outX <- glue ("\n\n{sep}")
  outX <- glue ( "{outX}\n.. {type}. {line1}." )
  if (!is.null (line2)) outX <- glue ( "{outX}\n.. {type}. {line2}." )
  if (!is.null (line3)) outX <- glue ( "{outX}\n.. {type}. {line3}." )
  outX <- glue( "{outX}\n{sep}\n")
  if (type != "ERROR") print (outX)
  invisible (outX)
}

### ====
### Section. Inputs
### ====

panda_initialize <- function (workspace.type) {
  # identiry workspace type ('modules' or 'pipelines')
  modules.workspace <<- ifelse (workspace.type == 'modules', TRUE, FALSE)
  
  # initialize system variables
  terra.wkspace <<- Sys.getenv( 'WORKSPACE_NAME' )
  google.bucket <<- Sys.getenv( 'WORKSPACE_BUCKET' )
  globals$project <<- Sys.getenv( 'WORKSPACE_NAMESPACE' )

  # check to make sure workspace does not have spaces or other special characters
  if (grepl ( '[^[:alnum:]|_|-]', terra.wkspace )) {
    error <- printX ("ERROR", "Invalid Workspace name", 
                     "Only letters, numbers, dashes and underscores are allowed")
    stop (error)
  }
  
  sys_home <- Sys.getenv( 'HOME' )
  if (sys_home == terra_notebook_home) {
    # this is running in Terra
    home <<- glue( "{sys_home}/{terra.wkspace}/edit/" )
  } else {
    # running (manually) elsewhere
    home <<- glue( "{Sys.getenv('PWD')}/{terra.wkspace}" )
    if (! dir.exists (home)) dir.create (home)   # in case doesn't exist
  }
  setwd (home)
  
  new.config <<- TRUE
  # check if a config is available, and the working directory is intact
  cfg <- defaults$panda_parameters_file_name
  check.cfg <- system (glue ("gsutil -q stat {google.bucket}/{cfg}"))
  if ( check.cfg==0 ) {
    # get config file
    system( glue( "gsutil cp {google.bucket}/{cfg} {cfg}" ) )
    # read and set parameter values
    p <- read_yaml (cfg)
    # check working directory for files
    ok <- TRUE
    for (f in unlist (c (p$typemap.gct, p$typemap.csv, p$typemap.gmt, p$typemap.yml))) {
      if (! file.exists (glue ("input/{f}")) ) {
        ok <- FALSE
        printX ("WARNING", glue ("Previous {cfg} exists, but input files are missing"),
                glue ("Existing {cfg} will be ignored"))
        break
      }
    }
    if (ok) {
      printX ("INFO", glue("Previous {cfg} exists"),
              glue("Input file map, groups and colors restored"),
              glue("Modify using appropriate Sections"))

      typemap.gct <<- p$typemap.gct
      typemap.csv <<- p$typemap.csv
      typemap.gmt <<- p$typemap.gmt
      typemap.yml <<- p$typemap.yml
      groups.cols <<- p$groups.cols
      groups.continuous <<- p$groups.continuous
      groups.colors <<- list()
      for ( group in names( p$groups.colors ) ) {
        groups.colors[[group]] <<- list()
        groups.colors[[group]]$vals <<- names (p$groups.colors[[group]])
        groups.colors[[group]]$colors <<- unlist (p$groups.colors[[group]])
      }
      new.config <<- FALSE
    }
  }
  print (DONE)
}

panda_datatypes <- function () {
  # print list of available data types
  print( glue( "{add_space( 1:length( cat.map ) )}: {cat.map}" ) )
}

map_my_files <- function( ext ){
  mapped <- list()
  for ( file in ext ){
    while (is.na (category <- as.numeric( readline( prompt = glue( ".. {file}: " ))))) {}
    if ( category == 0 ){
      cat.map <- c( cat.map, readline( prompt = "\n$$ Enter category name:" ) )
      category <- length( cat.map )
    } else if ( category < 0 || category > length(cat.map) )
      stop( "\n\n{sep}\n.. ERROR. Invalid entry. Try again.." )
    mapped[[cat.map[category]]] <- file
  }
  return ( mapped )
}

load_unzipped_files <- function(){
  setwd( home )
  ## locate uploaded file
  input.zip.name <- readline(
    prompt = "\n$$ Enter uploaded zip file name (test.zip): " )
  input.zip <<- glue( "{google.bucket}/{input.zip.name}" )
  zip.name <- tail( unlist( strsplit( input.zip, split = '/' ) ), 1 )
  system( glue( "gsutil cp {input.zip} {home}/." ) )
  if (!file.exists( zip.name )) 
    stop ( printX ("ERROR", glue("Zip file {zip.name} not found in workspace bucket")) )
  
  if( dir.exists( "input" ) )
    unlink( "input", recursive = T )
  dir.create( "input" )
  system( glue( "unzip -j {zip.name} -d input/" ) )
  setwd( glue( "{home}/input" ) )
  gcts <- list.files( pattern = "*.gct" )
  csvs <- list.files( pattern = "*.csv" )
  ymls <- list.files( pattern = "*.yaml" )
  gmts <- list.files( pattern = "*.gmt" )
  typemap.gct <<- map_my_files( gcts )
  typemap.csv <<- map_my_files( csvs )
  typemap.yml <<- map_my_files( ymls )
  if (is.null (typemap.yml$parameters)) 
    typemap.yml$parameters <<- process_file (defaults$parameters)
  typemap.gmt <<- map_my_files( gmts )
  if (is.null (typemap.gmt$ptmseaDB)) 
    typemap.gmt$ptmseaDB <<- process_file (defaults$ptmsea_db)
  if (is.null (typemap.gmt$gseaDB)) 
    typemap.gmt$gseaDB <<- process_file (defaults$gsea_db)
  
  # check if proteomics data should be normalized
  normalize.prot <<- y2true ("Does proteomics data need normalization?")
}

validate_annotation_table <- function (typemap.csv) {
  annot <- read.csv( typemap.csv$annotation, header = T,
                     stringsAsFactors = F, quote = '"' )
  all.groups <- colnames( annot )
  missing.cols <- setdiff( required.cols, all.groups )
  if( length( missing.cols ) > 0 )
    stop( printX ("ERROR", 
                  glue("Missing Columns. \n{paste0( missing.cols, collapse = ', ' )}")) )
  else
    printX ("INFO", "Sample annotation file successfully validated")
  
  return (annot$Sample.ID)
}

validate_data_table <- function (data.type, table.samples, sample.list) {
  common <- intersect (table.samples, sample.list)
  if (length(common) == 0) {
    # no common samples between data.type and annotation table
    error <- printX ("ERROR", glue("No samples found in {data.type}"),
                     glue("Harmonize samples names in annotation table and data tables"))
    stop( error )
  }
  if (length(common) < length(sample.list)/2) {
    # less than half of samples in annotation table found in data table
    printX ("WARNING", glue("Only {length(common)} (out of {length(sample.list)}) found in {data.type}"))
  } else {
    print( glue( "\n.. INFO. {data.type} successfully validated.\n" ) )
  }
}

validate_input <- function () {
  ## check input data against annotation table
  # start by checking if annotation table is present, and has required columns
  if (is.null (typemap.csv$annotation)) {
    stop( printX ("ERROR", "Sample annotation file missing" ) )
  }
  samples <- validate_annotation_table (typemap.csv)
  
  # check for sample id matches to ensure proper harmonization
  if ( (!exists ("typemap.gct")) || length (typemap.gct) == 0 )
    stop( printX ("ERROR", "No gct files specified") )
  if (!modules.workspace) {
    if ( all (names (typemap.gct) %in% proteome.types == FALSE) )
      stop( printX ("ERROR", "No proteomics dataset specified") )
    if ( ! all (required.genomics.types %in% names (typemap.gct)) )
      stop( printX ("ERROR", glue ("Genomics data ({paste (required.genomics.types, collapse='/')}) missing")) )
  }
  print( glue( "\n.. INFO. Validating sample IDs in other files ..\n" ) )
  flush.console ()  # without this, the display shows up later
  sapply (names (typemap.gct),   # gct files
          function (n) {
            d <- suppressMessages (parse.gctx (typemap.gct[[n]]))
            validate_data_table (n, d@cid, samples)
          })
  sapply (setdiff (names (typemap.csv), c('annotation', 'groups')),   # csv files
          function (n) {
            d <- read.csv (typemap.csv[[n]], header = T,
                           stringsAsFactors = F, quote = '"')
            validate_data_table (n, d$Sample.ID, samples)
          })
}

panda_input <- function(){
  new.config <<- TRUE
  load_unzipped_files()
  validate_input()
  print( DONE )
}

### ====
### Section. Groups
### ====

get_all_groups <- function( typemap.csv ){
  annot <- read_annot()
  all.groups <- colnames( annot )
  return( all.groups )
}

display_all_groups <- function( typemap.csv ){
  all.groups <- get_all_groups( typemap.csv )
  print( glue( "\n\n{sep}\n.. Annotations:\n\n" ) )
  print( glue( "{add_space( 1:length( all.groups ) )}: {all.groups}" ) )
  print( glue( "\n{sep}\n" ) )
  return( all.groups )
}

select_groups_case <- function( all.groups ){
  groups.cols <- c()
  stmt <- glue( "\n$$ Enter comma separated indexes or ranges " )
  stmt <- glue( "{stmt}(1:2, 5, 9, 12:14): " )
  selection <- readline( prompt = stmt )
  for (s in unlist (strsplit (selection, split=','))) {
    range <- unlist (strsplit (s, split = ':' ))
    lwend <- as.numeric( range[1] )
    if( length( range ) == 1 )
      upend <- lwend else
        upend <- as.numeric( range[2] )
    if (any (is.na (c (lwend, upend))) || 
        lwend > length(all.groups) || upend > length(all.groups)) {
      print ( glue( "\n.. Index {s} out of range, skipping ... \n"))
      next
    }
    groups.cols <- c( groups.cols, all.groups[lwend:upend] )
  }
  return( groups.cols )
}

select_groups <- function( all.groups ){
  groups.flag <- T
  groups.cols <- c()
  while ( groups.flag ){
    groups.cols <- c( groups.cols, select_groups_case(all.groups) )
    groups.flag <- y2true("$$ Continue adding groups?")
  }
  return( groups.cols )
}

verify_group_validity <- function( groups.cols, typemap.csv ){
  if (length (groups.cols) == 0) return (groups.cols)
  drop <- c()
  cont <- groups.continuous <- c()
  annot <- read_annot()
  for ( group.idx in 1:length( groups.cols ) ){
    groups.vals <- unique( annot[[groups.cols[group.idx]]] )
    if ( length( groups.vals ) > max.categories || length( groups.vals ) <= 1 ){  
      # drop if column has <= 1 unique values (ie column not present, or is identical for all samples)
      # if column has > max.categories, drop from groups.cols, 
      #   but, if numeric, treat as a continuous column and add to groups.continuous
      drop <- c( drop, group.idx )
      if (length( groups.vals ) > max.categories && is.numeric (groups.vals)) 
        cont <- c (cont, group.idx)
    }
  }
  if ( length( cont ) > 0 ) groups.continuous <- groups.cols[cont]
  if ( length( drop ) > 0 ) groups.cols <- groups.cols[-drop]     # do this last
  return( list (groups.cols=groups.cols, groups.continuous=groups.continuous) )
}

display_validated_groups <- function( groups.cols, groups.continuous, groups=FALSE ){
  if (!groups ) {
    warning <- glue( "Annotations (if any) with <2 categories were excluded." )
    warning <- glue( "\nAnnotations with > {max.categories} are considered continuous. ")
    printX ("WARNING", warning)
  }
  print( glue( "\n.. Selected groups:\n" ) )
  print( glue( "{add_space( 1:length( groups.cols ) )}: {groups.cols}" ) )
  if (length (groups.continuous) > 0) {
    print( glue( "\n\n.. Continuous groups:\n" ) )
    print( glue( "{add_space( 1:length( groups.continuous ) )}: {groups.continuous}" ) )
  }
}

panda_select_groups <- function(){
  groups.cols <<- select_groups( all.groups )
  if (length (groups.cols) == 0) {
    # no groups were selected -- choose all valid columns from annotation
    printX ("INFO", "No groups selected. Using all valid annotations")
    groups.cols <<- setdiff (all.groups, c (required.cols, ignore.cols))
  }
  final.cols <- verify_group_validity( groups.cols, typemap.csv )
  groups.cols <<- final.cols$groups.cols
  groups.continuous <<- final.cols$groups.continuous
}

panda_groups <- function(){
  new.config <<- TRUE
  keep <- FALSE
  if ( !is.null (typemap.csv$groups) ) {
    # groups file already specified
    keep <- y2true("Groups file already present. Keep it?")
    if (keep) {
      gr <- read.csv( typemap.csv$groups, header = T,
                      stringsAsFactors = F, quote = '"' )
      final.cols <- verify_group_validity( as.vector (unlist (gr)), typemap.csv )
      groups.cols <<- final.cols$groups.cols
      groups.continuous <<- final.cols$groups.continuous
    }
  } 
  
  if (keep==FALSE) {
    all.groups  <<- display_all_groups( typemap.csv )
    flush.console ()  # without this, the display shows up after the next function
    panda_select_groups ()
  }
  
  display_validated_groups( groups.cols, groups.continuous, groups=keep )
  panda_colors_defaults (display=FALSE)  ## always set default colors automatically
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

set_default_colors <- function( groups.cols, typemap.csv ){
  annot <- read_annot()
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

panda_colors_defaults <- function(display=TRUE){
  groups.colors <<- set_default_colors( groups.cols, typemap.csv )
  if (display) display_default_colors( groups.cols, groups.colors )
}

panda_colors_edit_index <- function( groups.colors ){
  for ( group.idx in 1:length( names( groups.colors ) ) ){
    group <- names( groups.colors )[group.idx]
    print( glue( "{sep}\n{add_space( group.idx )}] {group} =\n" ) )
    values <- groups.colors[[group]]$vals
    print( glue( "\t{add_space( 1:length( values ) )}: {values}" ) )
  }
}

change_default_colors <- function( groups.cols, all.groups, groups.colors ){
  new.config <<- TRUE
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
      value.flag <- y2true( glue("\n$$ Continue editing colors for {change.group}?") )
    }
    groups.flag <- y2true("\n$$ Continue modifying defaults for other groups?")
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
### Section. COSMO labels
### ===

# initialize this
cosmo.params <- list(run_cosmo = FALSE)

# function to get used input for cosmo
select_COSMO_attributes <- function() { 
  
  # initialize again
  cosmo.params <<- list(run_cosmo = FALSE)
  
  annot <- read_annot()
  potential_cols <- setdiff(names(annot), ignore.cols)
  
  # get valid columns
  valid_columns <- c()
  for (col in potential_cols) {
    if (!(col %in% names(annot))) {
      warning(paste0("Sample label '", col, "' is not in sample annotation file"))
    } else if (min(base::table(annot[, col])) < min(10, dim(annot)[1] / 5)) {
      warning(paste0("Sample label '", col, "' is not well-balanced. It is being excluded."))
    } else if (any(is.na(annot[, col]))) {
      warning(paste0("Sample label '", col, "' has NAs. It is being excluded."))
    } else if (length(unique(annot[, col])) < 2) {
      warning(paste0("Sample label '", col, "' only has one level. It is being excluded."))
    } else if (length(unique(annot[, col])) > 2) {
      warning(paste0("Sample label '", col, "' being excluded because COSMO does not handle classes with more than 2 levels."))
    } else {
      valid_columns <- c(valid_columns, col)
    }
  }
  
  # make prompt
  if (length(valid_columns > 0)) {
    prompt <- paste("VALID ATTRIBUTES:",
                    paste(paste(' *', valid_columns), collapse = '\n'), 
                    sep,
                    "Select sample label column(s) for COSMO: ",
                    sep = '\n')
  } else {
    cat("NO VALID ATTRIBUTES FOUND. COSMO will not be run.\n", DONE)
    return()
  }
  
  # get user input
  user_columns <- readline(prompt)
  
  # validate column selection
  user_columns <- strsplit(gsub(' ', '', user_columns), ',')[[1]]
  
  invalid_user_columns <- setdiff(user_columns, valid_columns)
  if (length(invalid_user_columns > 0)) {
    message(paste("Invalid attribute:", 
                  paste(invalid_user_columns, collapse = ', ')))
  }
  
  valid_user_columns <- intersect(user_columns, valid_columns)
  
  
  cat(sep, '\n')
  
  if (length(valid_user_columns) > 0) {
    cat("ATTRIBUTE SELECTION:\n", paste(paste(' *', valid_user_columns), collapse = '\n'), sep = '')
    cat("\n\nCOSMO WILL BE RUN.\n")
    run_cosmo <- TRUE
  } else {
    cat("NO ATTRIBUTES SELECTED. COSMO will not be run.\n")
    run_cosmo <- FALSE
  }
  
  cosmo.params <<- list(run_cosmo = run_cosmo,
                        sample_label = paste(valid_user_columns, collapse = ','))
  
  cat(DONE)
}

### ===
### Section. Sample sets
### ===

sample_set_sanity_check <- function( sample.sets ){
  command <- glue( "fissfc sset_list -w {terra.wkspace} -p {globals$project}" )
  exist.ss <- system( command, intern = T )
  
  # delete the 'all' sample set so that it can be re-created with updated parameters
  # but don't delete yet, in case sample set definition is aborted
  delete.all <- ifelse ("all" %in% exist.ss, TRUE, FALSE)
    
  exist.ss <- setdiff (exist.ss, 'all')  # always include 'all' -- keeps samples synchronized
  if ( !( identical( character(), exist.ss ) ) ){
    overlap.flags <- names( sample.sets ) %in% exist.ss
    if ( any( overlap.flags ) ){
      overlap.ssets <- names( sample.sets )[overlap.flags]
      warning <- glue ("WARNING. Sample sets -- {overlap.ssets} -- already exist.")
      print( glue( "\n\n{sep}\n.. {warning}\n" ) )
      flush.console()
      overwrite <- y2true ("Overwrite sample sets?")
      if (overwrite) {
        # delete existing sample sets, including 'all'
        # "echo y" automates the "y(es)/no" prompt needed to confirm deletion
        sapply (overlap.ssets, 
                function (s) 
                  system ( glue ( "echo y | fissfc sset_delete -w {terra.wkspace} -p {globals$project} -e {s}" ))
        )
      } else {
        print( glue( "\n\n.. Aborting sample subset definition\n") )
        sample.sets <- list()
        return (sample.sets)
      }
    }
  }
  # delete 'all' sample set, if it exists (and sample set definition was not aborted)
  if (delete.all) 
    system ( glue ( "echo y | fissfc sset_delete -w {terra.wkspace} -p {globals$project} -e all" ))
  
  return( sample.sets )
}

define_sample_sets <- function( all.groups, typemap.csv ){
  annot <- read_annot()
  all.groups <- colnames (annot)
  sample.sets <- list()
  sample.sets$all$fil_col <- "Type"
  sample.sets$all$fil_val <- glue(
    "{paste0( unique( annot[[sample.sets$all$fil_col]] ), collapse = ';' )}" )
  # add pathway databases and parameters file as set-specific parameters
  for (x in names (typemap.gmt))  
    sample.sets$all[[x]] <- typemap.gmt[[x]]
  # this will be the concatenated parameter file
  for (x in names (typemap.yml))  
    sample.sets$all[[x]] <- typemap.yml[[x]]
  
  
  set.flag <- y2true( "\n$$ Add additional sample subsets?" )
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
    fil_val <- unlist (strsplit (fil_val, split=';'))
    if ( !( all (fil_val %in% values) ) )
      stop( glue( "\n\n{sep}\n.. ERROR. Value entered is incorrect." ) )
    sample.sets[[name]]$fil_col <- fil_col
    sample.sets[[name]]$fil_val <- fil_val

    # add pathway databases and parameters as set-specific parameters for every set
    for (x in names (typemap.gmt))  
      sample.sets[[name]][[x]] <- typemap.gmt[[x]]
    # this will be the concatenated parameter file
    for (x in names (typemap.yml))  
      sample.sets[[name]][[x]] <- typemap.yml[[x]]

    # add.param.flag <- y2true("$$ Add set-specific parameters?")
    # while( add.param.flag ){
    #   param.name <- as.character( readline(
    #     prompt = "$$ Enter Parameter Name: " ) )
    #   param.value <- as.character( readline(
    #     prompt = "$$ Enter parameter value: " ) )
    #   sample.sets[[name]][[param.name]] <- param.value
    #   add.param.flag <- y2true( "$$ Continue adding more set parameters?" )
    # }
    set.flag <- y2true("\n$$ Continue adding more sets?")
  }
  sample.sets <- sample_set_sanity_check( sample.sets )
  if ( length( names( sample.sets ) ) == 0 )
    print( glue( "\n{sep}\n.. No new sample sets will be added." ) ) else{
      print(
        glue( "\n{sep}\n.. Sample sets to be added to Terra Workspace: " ) )
      print( glue( "\n.. {paste0( names( sample.sets ), collapse = ', ' )}" ) )
  }
  return( sample.sets )
}

panda_sample_subsets <- function(){
  new.config <<- TRUE
  sample.sets <<- define_sample_sets( all.groups, typemap.csv )
  print( DONE )
}

### ===
### Section. Build YAML
### ===

panda_finalize <- function (internal=FALSE) {
  # perform sanity checks (in case this function is called without going through other steps)
  if (!exists ("typemap.csv") || is.null (typemap.csv$annotation))
    stop( glue( "\n\n{sep}\n.. ERROR. Sample annotation file missing. Run panda_input().\n{sep}\n" ) )
  if (!exists ("groups.cols") || length (groups.cols) == 0) 
    stop( glue( "\n\n{sep}\n.. ERROR. No groups selected. Run panda_groups()." ) )
  if (!exists ("sample.sets") || length (sample.sets) == 0) 
    stop( glue( "\n\n{sep}\n.. ERROR. No sample subsets. Run panda_sample_subsets()." ) )

  # write panda config file
  # also creates a combined parameter file for use in the PANOPLY pipeline
  # created yaml files are transferred to the workspace google bucket
  lines <- list()
  lines$TEDMINT <- glue ("{panda}/bin")
  lines$freeze.path <-
    glue( "{home}/input/" )
  lines$wkspace <- terra.wkspace
  lines$globals <- globals
  lines$typemap.gct <- typemap.gct
  lines$typemap.csv <- typemap.csv
  lines$typemap.gmt <- typemap.gmt
  lines$typemap.yml <- typemap.yml
  lines$groups.cols <- groups.cols
  lines$groups.cols.continuous <- groups.continuous
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
  lines$normalize.proteomics <- normalize.prot
  lines$cosmo.params <- cosmo.params
  
  
  # output config for panda and copy to google bucket
  setwd( home )
  cfg <- defaults$panda_parameters_file_name
  write_yaml( x = lines, file = cfg,
              handlers = list (logical=function (x) {
                result <- ifelse (x, "TRUE", "FALSE")
                class (result) <- "verbatim"
                return (result)
              }))
  system( glue( "gsutil cp {home}/{cfg} {google.bucket}/{cfg}" ) )
  new.config <<- FALSE
  # concatenate master parameter (or input parameter) file -- will be included in sample_set
  original.params <- glue ("{home}/input/{typemap.yml$parameters}.orig")
  params <- glue ("{home}/input/{typemap.yml$parameters}")
  if (! file.exists (original.params)) file.copy (params, original.params)
  catcmd <- glue ("cat {cfg} {original.params} > {params}")
  system (catcmd)
  if (!internal) print( DONE )
}

### ===
### Section. PANDA tool
### ===

run_panda <- function(){
  # run the PANDA tool which transfers data to the google bucket, 
  # sets up metadata table in the workspace and creates sample 
  # subsets and associated data files
  
  # if new.config is TRUE, run panda_write_config to create (or update) config file
  if (new.config) panda_finalize (internal=TRUE)
  
  setwd (home)
  runspace <- "runspace"
  if ( dir.exists( runspace ) )
    unlink( runspace, recursive = T )
  dir.create( runspace )
  invisible( file.copy( glue ("{panda}/bin/Makefile"), glue( "{runspace}/Makefile" ) ) )
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
  for ( f in c (typemap.yml, typemap.gmt) ) {
    system( glue (
      "ln -s ../../input/{f} {f}" ))
  }
  setwd( "../" )
  print( WAIT )
  flush.console()  # without this, the display shows up later
  system( glue( "make groups" ) )
  system( glue( "make terrainit" ) )
  system( glue( "make panda-samples" ) )  # always call -- will upload only missing samples
  system( glue( "make panda-sets" ) )
  print( DONE )
}


