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
suppressMessages(library( dplyr )); # load quietly
library( pacman );
library( RColorBrewer )
p_load( scales );
p_load( glue );
p_load( yaml );
p_load( ids );

sink(); source('/prot/proteomics/Projects/R-utilities/map-to-genes.r')

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
defaults$ptmsea_db <- glue ("{panda}/defaults/ptm.sig.db.all.flanking.human.v1.9.0.gmt")
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

trim <- function( x ) { gsub( "^\\s+|\\s+$", "", x ) }

y2true <- function( prompt ) {
  # displays prompt and reads in y/n response
  valid_choice <- function (ch) ch %in% c ('Y', 'y', 'N', 'n')
  while (! valid_choice (choice <- readline (prompt = glue ("{prompt} (y/n): "))) ) {}
  if ( choice == "y" || choice == "Y")
    flag <- TRUE else flag <- FALSE
  return( flag )
}

### valid_choice() assesses whether an input (choice) is valid, according to default and a custom
valid_choice <- function (choice, FUN=NULL, exit_commands = c("q", "Q", "quit", "QUIT", "exit", "EXIT"),
                          allow_na=FALSE, allow_empty=FALSE, allow_zero=FALSE) {
  flush.console() # automatically flush console
  
  if (is.null(FUN)) { # if no custom_function was provided
    FUN = function(...) { return(TRUE)} # set FUN to automatically return TRUE
  } else { FUN <- match.fun(FUN) } # coerce FUN input into function
  
  if (choice %in% exit_commands) return(1) # if quit command was triggered, return a NUMERIC true
  else # otherwise, check that choice fulfils ALL of the following conditions:
    return(all(!is.na(choice) || allow_na, # choice is not NA
               choice!="" || allow_empty, # choice is not an empty string
               if(!is.na(suppressWarnings(as.numeric(choice)))) {as.numeric(choice)!=0} || allow_zero, # if choice can be coerced into number, check if it's nonzero
               FUN(choice))) # choice fulfills custom_condition function
}


# readline(), but only accepts valid inputs & recognizes exit-commands
smart_readline <- function (prompt="$$ Input: ", trim_input=TRUE,
                            custom_condition=NULL, custom_warning="Invalid input, try again.",
                            exit_commands = c("q", "Q", "quit", "QUIT", "exit", "EXIT", "cancel", "CANCEL"),
                            exit_message = "Quit condition triggered.", exit_function = NULL,
                            allow_na=FALSE, allow_empty=FALSE, allow_zero=FALSE) {
  if (!exists('trim')) { trim <- function( x ) { gsub( "^\\s+|\\s+$", "", x ) } } # define trim(), if missing
  while (! valid_choice ( choice <- readline(prompt) %>% # prompt until the choice is valid
                            { ifelse(trim_input, trim(.), . ) } , #trim choice, if toggled on
                          FUN = custom_condition, # apply custom condition; this function should return a single TRUE/FALSE
                          exit_commands = exit_commands,
                          allow_na=allow_na, allow_empty=allow_empty, allow_zero=allow_zero)) { cat(custom_warning) } 
  if ( is.numeric(valid_choice(choice)) ) { # if we've triggered the quit condition
    if (!is.null(exit_message)) { cat(exit_message) } # print exit message
    if (!is.null(exit_function)) { try(exit_function) } # run exit function
    return(NULL) } # return NULL
  
  return(choice) # otherwise, return choice
}

printX <- function (type, line1, line2=NULL, line3=NULL, invisible=FALSE) {
  outX <- glue ("\n\n{sep}")
  outX <- glue ( "{outX}\n.. {type}. {line1}." )
  if (!is.null (line2)) outX <- glue ( "{outX}\n.. {type}. {line2}." )
  if (!is.null (line3)) outX <- glue ( "{outX}\n.. {type}. {line3}." )
  outX <- glue( "{outX}\n{sep}\n")
  if (type != "ERROR" && !invisible) print (outX)
  invisible (outX)
}

stnd_custom_warning<-printX("ERROR", "Invalid input. Please try again") # standard exit message, for use with smart_readline()
stnd_exit_message<-printX("CANCELLED","Previous inputs not saved", invisible = TRUE) # standard exit message, for use with smart_readline()

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
      groups.cols.continuous <<- p$groups.cols.continuous
      groups.colors <<- list()
      for ( group in names( p$groups.colors ) ) {
        groups.colors[[group]] <<- list()
        groups.colors[[group]]$vals <<- names (p$groups.colors[[group]])
        groups.colors[[group]]$colors <<- unlist (p$groups.colors[[group]])
      }
      new.config <<- FALSE
    }
    
    if (!is.null(p$groups.colors)) { # if we have colors from the previous dataset
      # printX ("INFO", glue ("Colorscheme found in previous {cfg}"))
      # save old colors, in case they should be restored
      groups.colors.old <<- list()
      for ( group in names( p$groups.colors ) ) {
        groups.colors.old[[group]] <<- list()
        groups.colors.old[[group]]$vals <<- names (p$groups.colors[[group]])
        groups.colors.old[[group]]$colors <<- unlist (p$groups.colors[[group]])
      }
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
    choice = smart_readline(prompt = glue( ".. {file}: " ), allow_zero=TRUE,
                            custom_condition = function(category) { !is.na(as.numeric(category)) && between(as.numeric(category), 0, length(cat.map)) }, # check if input is numerical, and in the appropriate range
                            custom_warning = printX("ERROR", glue("Invalid category index. Please enter an index in the range {paste(0, 'to', length(cat.map))}"), invisible = TRUE),
                            exit_message = stnd_exit_message)
    if (is.null(choice)) stop() # if we triggered a quit condition, stop
    category <- as.numeric(choice)
    if ( category == 0 ){
      cat.map <- c( cat.map, readline( prompt = "\n$$ Enter category name:" ) )
      category <- length( cat.map )
    }
    mapped[[cat.map[category]]] <- file
  }
  return ( mapped )
}

load_unzipped_files <- function(){
  setwd( home )
  ## locate uploaded file
  process_zip <- function(input.zip.name) {
    input.zip <<- glue( "{google.bucket}/{input.zip.name}" )
    system( glue( "gsutil cp {input.zip} {home}/." ) )
    zip.name <- tail( unlist( strsplit( input.zip, split = '/' ) ), 1 )
    return(zip.name)
  }
  input.zip.name <- smart_readline( prompt = "\n$$ Enter uploaded zip file name (test.zip): ",
                                    custom_condition = function(input.zip.name) { file.exists( process_zip(input.zip.name) ) },
                                    custom_warning = printX ("ERROR", glue("Zip file not found in workspace bucket")),
                                    exit_message = stnd_exit_message) # ensure that zip file exists
  if (is.null(input.zip.name)) stop() # if we triggered a quit condition, stop
  zip.name = process_zip(input.zip.name)
  
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

validate_sample_id <- function (data.type, table.samples, sample.list) {
  common <- intersect (table.samples, sample.list)
  if (length(common) == 0) {
    # no common samples between data.type and annotation table
    error <- printX ("ERROR", glue("No samples found in {toupper(data.type)}"),
                     glue("Harmonize samples names in annotation table and data tables"))
    stop( error )
  }
  if (length(common) < length(sample.list)/2) {
    # less than half of samples in annotation table found in data table
    printX ("WARNING", glue("Only {length(common)} (out of {length(sample.list)}) found in {toupper(data.type)}"))
  } else {
    print( glue( "\n.. INFO. {toupper(data.type)} successfully validated.\n" ) )
  }
  flush.console ()  # without this, the display shows up later
}

gene_id_column_select <- function(ome, gct, gene.id.col.default) {
  rdesc = gct@rdesc
  rdesc_names = names(rdesc)
  #cat(paste0(sep, '\n'))
  cat(paste0('\n'))
  cat(paste(glue("{toupper(ome)} Annotation Columns:\t"),
            paste(glue("{rdesc_names}" ), collapse = ', '), '\n'))
  
  # Prompt user for the relevant column
  valid_id = FALSE # initialize
  while ( !valid_id ) { # until we select a valid ID
    gene.id.col = smart_readline( prompt = glue("\n$$ Enter the column from the {toupper(ome)} data that has HUGO Gene Symbols:"),
                                  custom_condition = function(column) { column %in% rdesc_names },
                                  custom_warning = printX ("ERROR", glue("Invalid column name, please try again")),
                                  exit_message = "\n.. RETURNING\n" )
    if(is.null(gene.id.col)) break() # go back to id.col.type prompt
    
    # check if IDs are valid
    valid_id = TRUE
    tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = rdesc[[ gene.id.col ]] ,
                                   keytype="SYMBOL", columns = "SYMBOL"),
             error = function(cond) {
               cat(printX("ERROR", glue("IDs in column '{gene.id.col}' are not valid HUGO Gene Symbols, please try again")))
               valid_id <<- FALSE }) # set valid to FALSE
  }
  
  # If we got a successful column selection
  if (valid_id) { #select column w/ gene IDs
    printX("INFO", glue("Using column '{gene.id.col}' as {gene.id.col.default} column for {toupper(ome)} data"))
    gct@rdesc[[gene.id.col.default]] <- gct@rdesc[[ gene.id.col ]] # overwrite default gene.id.col with new gene.id.col
    write.gct(gct, typemap.gct[[ome]], appenddim=FALSE) # overwrite GCT file
  }
  
  flush.console()
  return(valid_id)
}

gene_id_column_create <- function(ome, gct, gene.id.col.default, protein.id.col) {
  rdesc = gct@rdesc
  rdesc_names = names(rdesc)
  # check that we have protein IDs at all
  if (!( protein.id.col %in% rdesc_names )) { # if the protein.id.col is missing
    cat(printX("ERROR", glue("Protein ID column '{protein.id.col}' missing from {ome} GCT file '{typemap.gct[[ome]]}'
                               ..        Please ensure that all proteomics datasets have a protein ID column in the row-description")))
    stop() # quit
  }
  # print valid protein ID types
  # keytype.opt <- AnnotationDbi::keytypes(org.Hs.eg.db) # full list of ID types
  keytype.opt <- c('ENSEMBLPROT', 'IPI', 'REFSEQ', 'PFAM', 'PROSITE', 'UNIPROT') # ONLY protein ID types
  #cat(paste0(sep, '\n'))
  cat(paste0('\n'))
  cat(paste("Valid Protein Keytypes:",
            paste(glue( "{add_space( 1:length( keytype.opt ) )}: {keytype.opt}" ), collapse = '\n'), '\n',
            sep = '\n'))
  
  # Prompt user for the Protein ID Keytype
  valid_id = FALSE # initialize value
  while (!valid_id) { # until we select a valid ID column
    protein.id.type.index = smart_readline( prompt = glue("\n$$ Enter the keytype (from above) that corresponds to protein IDs in column '{protein.id.col}' (format '{rdesc[['id']][1]}'):"),
                                            custom_condition = function(keytype) { keytype %in% 1:length(keytype.opt) },
                                            custom_warning = printX ("ERROR", glue("Invalid keytype index, please enter one of the keytypes indexes above")),
                                            exit_message = "\n.. RETURNING\n" )
    if(is.null(protein.id.type.index)) break()  # go back to id.col.type prompt
    protein.id.type = keytype.opt[[as.numeric(protein.id.type.index)]]
    
    protein.ids = rdesc[[ protein.id.col ]]
    if ( y2true("Should version numbers be trimmed?") ) {
      protein.ids = sub ("(.*?)\\..*?$", "\\1", protein.ids) # assumes an ID with a version number will be of format <ID>.<version_number> (e.g. ENSP00000301740.8)
      cat(glue("Format of trimmed IDs: '{protein.ids[[1]]}'\n"))
    }
    
    
    valid_id = TRUE # initialize value
    gene.ids = tryCatch(
      map_id(protein.ids, keytype_from = protein.id.type, keytype_to = "SYMBOL"), # attempt to map to symbols
      error = function(cond) {
        cat(printX("ERROR", glue("Protein IDs in column '{protein.id.col}' are not valid {protein.id.type} IDs, please try again")))
        valid_id <<- FALSE })
  }
  
  if (valid_id) { # convert protein IDs to gene IDs
    printX("INFO", glue("Protein IDs in column '{protein.id.col}' was successfully converted from {protein.id.type} IDs to Hugo Gene Symbols for {toupper(ome)} data"))
    gct@rdesc[[gene.id.col.default]] <- gene.ids # overwrite default gene.id.col with new gene.id.col
    write.gct(gct, typemap.gct[[ome]], appenddim=FALSE) # overwrite GCT file
  }
  
  flush.console()
  return(valid_id)
}

validate_gene_id <- function(ome) {
  suppressMessages(require(org.Hs.eg.db)) # read in org.Hs.eg.db, to let us check that gene id column has Hugo Gene Symbols
  gene.id.col.default = read_yaml(typemap.yml$parameters)$global_parameters$gene_mapping$gene_id_col
  # gene.id.col.default = "testing"
  protein.id.col = read_yaml(typemap.yml$parameters)$global_parameters$gene_mapping$protein_id_col
  
  gct = suppressMessages (parse.gctx (typemap.gct[[ome]]))
  rdesc = gct@rdesc
  rdesc_names = names(rdesc)
  
  # Check for default gene.id.col
  valid_id = TRUE # initialize value
  if ( gene.id.col.default %in% rdesc_names ) { # if the default gene.id.col exists
    print( glue( "\n.. INFO. Default Gene ID column '{gene.id.col.default}' detected in {toupper(ome)} data.\n" ) )
    tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = rdesc[[ gene.id.col.default ]] , # check if they're valid symbols
                                   keytype="SYMBOL", columns = "SYMBOL"),
             error = function(cond) {
               cat(printX("WARNING", glue("IDs in column '{gene.id.col.default}' are not valid HUGO Gene Symbols")))
               valid_id <<- FALSE })
  } else {
    printX("WARNING", glue("Default Gene ID column '{gene.id.col.default}' was NOT detected in {toupper(ome)} data"))
    valid_id = FALSE
  }
  
  # If the default gene.id.col wasn't found / wasn't valid
  if ( !valid_id ) {
    while ( !valid_id ) { # until we select a valid ID
      # check whether user is inputting a gene.id.col, or a prot.id.col and prot.id.type
      id.col.type = smart_readline( prompt = paste(glue("\n$$ To create a Gene ID column for {toupper(ome)}, please choose to either:"),
                                                   "\t1) Select an existing annotation column with HUGO Gene Symbols",
                                                   "\t2) Convert protein IDs to HUGO Gene Symbol \n",
                                                   sep = "\n"),
                                    custom_condition = function(input) { input %in% c(1,2) },
                                    custom_warning = stnd_custom_warning,
                                    exit_message = stnd_exit_message )
      if(is.null(id.col.type)) stop()
      
      # attempt to create Gene ID Column
      if (id.col.type==1)
        valid_id = gene_id_column_select(ome, gct, gene.id.col.default) #select column w/ gene IDs
      else if (id.col.type==2) 
        valid_id = gene_id_column_create(ome, gct, gene.id.col.default, protein.id.col) # convert protein IDs to gene IDs
      else { printX("ERROR", "This shouldn't have happened!"); stop() }
    }
  }
    
  flush.console()
  # nothing to return
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

  printX ("INFO", "Validating sample IDs in all files")
  flush.console ()  # without this, the display shows up later
  sapply (names (typemap.gct),   # gct filenames
          function (n) {
            d <- suppressMessages (parse.gctx (typemap.gct[[n]]))
            validate_sample_id (n, d@cid, samples)
          })
  sapply (setdiff (names (typemap.csv), c('annotation', 'groups')),   # csv files
          function (n) {
            d <- read.csv (typemap.csv[[n]], header = T,
                           stringsAsFactors = F, quote = '"')
            validate_sample_id (n, d$Sample.ID, samples)
          })
  
  # check for gene id columns in all GCT files
  printX ("INFO", "Validating gene IDs in GCT files")
  flush.console ()  # without this, the display shows up later
  sapply (names (typemap.gct),   # gct filenames
          function (ome) {
            validate_gene_id (ome)
          })
}

panda_input <- function(){
  new.config <<- TRUE
  load_unzipped_files()
  validate_input()
  print( DONE )
}


### ====
### Section. Toggle Processing
### ====

validate_flanking_sequence <- function() {
  flush.console()
  ome = 'phosphoproteome'
  valid_sequence_pattern="[A-z-]{7}[sty][A-z-]{7}"
  
  gct = suppressMessages(cmapR::parse.gctx(typemap.gct[[ome]]))
  rdesc = gct@rdesc
  rdesc_names = names(rdesc)
  
  flanking.seq.col.default = read_yaml(typemap.yml$parameters)$panoply_preprocess_gct$seqwin_column
  # flanking.seq.col.default = "testing"
  
  
  # Check for default gene.id.col
  valid_id = TRUE # initialize value
  if ( flanking.seq.col.default %in% rdesc_names ) { # if the default flanking.seq.col.default exists
    print( glue( "\n.. INFO. Default Flanking Sequence column '{flanking.seq.col.default}' detected in {toupper(ome)} data.\n" ) )
    if ( ! any( grepl(valid_sequence_pattern, rdesc[[flanking.seq.col.default]]) ) ) { # check that it has valid flanking sequences
      printX("WARNING", glue("IDs in column '{flanking.seq.col.default}' are not valid Flanking Sequences"))
      valid_id = FALSE
    } 
  } else {
    printX("WARNING", glue("Default Flanking Sequence column '{flanking.seq.col.default}' was NOT detected in {toupper(ome)} data"))
    valid_id = FALSE
  }
  flush.console()
  
  if (!valid_id) {
    # print column options
    cat(paste(glue("{toupper(ome)} COLUMNS:"),
              paste( glue( "{add_space( 1:length( rdesc_names ) )}: {rdesc_names}" ), collapse = '\n'), 
              sep, '\n',
              sep = '\n'))
    
    while (!valid_id) { # until we select a valid flanking sequence column
      flush.console()
      flanking.seq.col.index = smart_readline( prompt = glue("\n$$ Enter the column from the {toupper(ome)} data that has Flanking Sequences:"),
                                               custom_condition = function(column) { column %in% 1:length(rdesc_names) },
                                               custom_warning = printX ("ERROR", glue("Invalid column index, please try again")),
                                               exit_message = stnd_exit_message )
      if( is.null(flanking.seq.col.index) ) stop()
      flanking.seq.col = rdesc_names[ as.numeric(flanking.seq.col.index) ] # set global flanking.seq.col parameter
      
      # check if IDs are valid
      if ( any( grepl(valid_sequence_pattern, rdesc[[flanking.seq.col]]) ) ) { # check if ANY entries are valid flanking sequences
        valid_id = TRUE
      } else {
        cat(printX("ERROR", glue("Entries in column '{flanking.seq.col}' are not valid flanking-sequences, please try again")))
      }
    }
    
    # If we got a successful column selection
    if (valid_id) { #select column w/ gene IDs
      printX("INFO", glue("Using column '{flanking.seq.col}' as {flanking.seq.col.default} column for {toupper(ome)} data"))
      gct@rdesc[[flanking.seq.col.default]] <- gct@rdesc[[ flanking.seq.col ]] # overwrite default gene.id.col with new gene.id.col
      write.gct(gct, typemap.gct[[ome]], appenddim=FALSE) # overwrite GCT file
    } else { stop(printX("ERROR", "Something has gone terribly wrong.")) }
    
  }
  
  flush.console()
}

panda_preprocessing <- function() {
  new.config <<- TRUE
  
  ### check if proteomics data should be normalized / filtered
  normalize.prot <<- y2true ("Does proteomics data need normalization?")
  filter.prot <<- y2true ("Does proteomics data need filtering?")
  
  cat(paste0(sep, '\n'))
  
  ### PTM-SEA / Flanking Sequence Column
  if ( 'phosphoproteome' %in% names(typemap.gct) && !is.null(typemap.gmt) ) { # if we have phosphoproteome data
    run.ptmsea <<- y2true("Phosphoproteome data detected. Should PTM-SEA be run?")
    
    if (run.ptmsea) validate_flanking_sequence()
    else run.ptmsea <<- FALSE
  }
  
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
  
  selection <- smart_readline( prompt = stmt,
                               custom_condition = function(selection) {
                                 unlist (strsplit (selection, split=',')) %>% trim() %>% # strsplit ranges and trim
                                   grepl("^[0-9]+?(:[0-9]+?)?$", .) %>% all() }, # check that all ranges are in a valid format, i.e. number(s), or number(s):number(s)
                               custom_warning = printX ("ERROR", glue("Invalid range formatting, please try again")),
                               exit_message = stnd_exit_message) # ensure that zip file exists
  if (is.null(selection)) stop() # if we triggered a quit condition, stop
  
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
  cont <- groups.cols.continuous <- c()
  annot <- read_annot()
  for ( group.idx in 1:length( groups.cols ) ){
    groups.vals <- sort(unique( annot[[groups.cols[group.idx]]] ))
    if ( length( groups.vals ) > max.categories || length( groups.vals ) <= 1 ){  
      # drop if column has <= 1 unique values (ie column not present, or is identical for all samples)
      # if column has > max.categories, drop from groups.cols, 
      #   but, if numeric, treat as a continuous column and add to groups.cols.continuous
      drop <- c( drop, group.idx )
      if (length( groups.vals ) > max.categories && is.numeric (groups.vals)) 
        cont <- c (cont, group.idx)
    }
  }
  if ( length( cont ) > 0 ) groups.cols.continuous <- groups.cols[cont]
  if ( length( drop ) > 0 ) groups.cols <- groups.cols[-drop]     # do this last
  return( list (groups.cols=groups.cols, groups.cols.continuous=groups.cols.continuous) )
}

display_validated_groups <- function( groups.cols, groups.cols.continuous, groups=FALSE ){
  if (!groups ) {
    warning <- glue( "Annotations (if any) with <2 categories were excluded." )
    warning <- glue( "\nAnnotations with > {max.categories} are considered continuous. ")
    printX ("WARNING", warning)
  }
  print( glue( "\n.. Selected groups:\n" ) )
  print( glue( "{add_space( 1:length( groups.cols ) )}: {groups.cols}" ) )
  if (length (groups.cols.continuous) > 0) {
    print( glue( "\n\n.. Continuous groups:\n" ) )
    print( glue( "{add_space( 1:length( groups.cols.continuous ) )}: {groups.cols.continuous}" ) )
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
  groups.cols.continuous <<- final.cols$groups.cols.continuous
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
      groups.cols.continuous <<- final.cols$groups.cols.continuous
    }
  } else {
    all.groups  <<- display_all_groups( typemap.csv )
    flush.console ()  # without this, the display shows up after the next function
    panda_select_groups ()
  }
  
  display_validated_groups( groups.cols, groups.cols.continuous, groups=keep )
  
  # set groups colors
  if ( exists("groups.colors.old") && #old color scheme exists
       y2true ("Color scheme found in previous config-file. Attempt to restore previous color scheme?") ) { #user would like to restore color scheme
    groups.colors <<- restore_config_colors ( groups.colors.old, typemap.csv ) # restore old colors, to whatever extent is possible
  } else {
    groups.colors <<- set_default_colors( groups.cols, typemap.csv )  # set default colors automatically
  }
  
  print( DONE )
}


### ====
### Section. Colors
### ====

display_colors <- function( groups.cols, groups.colors ){
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

panda_display_current_colors <- function() {
  if (is.null(groups.colors)) {
    printX ("WARNING", glue ("Color scheme is missing. Please run panda_set_default_colors() below, before running this block."))
  } else {
    display_colors( groups.cols, groups.colors )
  }
}

set_default_colors <- function( groups.cols, typemap.csv ){
  annot <- read_annot()
  
  # utility function, allows pseudorandom choice based on a provided string, which is used as a temporary seed
  eval_with_seed <- function(expression, seed_string="randomstring") {
    old_seed = .Random.seed # pull current seed
    on.exit({.Random.seed <<- old_seed}) # reset seed back to normal on exit
    
    set.seed( sum(as.numeric(charToRaw(seed_string))) ) # set seed, based on string passed in
    eval( expression ) # match function
  }
  
  # Paul Tol Color Palattes-- colorblind safe! (https://personal.sron.nl/~pault/)
  qual.pals <- list("Bright" = c('#4477AA', '#66CCEE', '#27B13E', '#CCBB44', '#EE6677', '#AA3377'), # green (originally '#228833') was modified to be more distinct from blue under tritanopia colorblindness
                    "Vibrant" = c('#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377'),
                    "Muted" = c( '#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499'),
                    "HighContrast" = c('#004488', '#DDAA33', '#BB5566'))

  # if adding a new palate, please add LIGHTER colors SECOND, to ensure consistent formatting
  pair.pals <- list("BrightPaired" = c('#4477AA', "#A8DBFF", '#66CCEE', "#CAFFFF", '#27B13E', "#8BFFA2", '#CCBB44', "#FFFFA8", '#EE6677', "#FFCADB", '#AA3377', "#FF97DB"),
                    "VibrantPaired" = c('#0077BB', "#64DBFF", '#009988', "#64FDEC", '#EE7733', "#FFDB97", '#EE3377', "#FF97DB" ),
                    "HighContrastPaired" = c('#004488', '#6699CC', '#997700', '#EECC66', '#994455', '#EE99AA'))
  pair.pals <- eval_with_seed( { lapply(pair.pals, # for each palate
                                        function(color_vec) {
                                          reorder_index = rep(sample( 1:(length(color_vec)/2) *2) , each=2) + rep(c(-1,0), length(color_vec)/2) # shuffle the pairs (TOGEHTER)
                                          color_vec[reorder_index]} ) } , # reorder palate to be shuffled
                               seed_string = "paultolcolorblindsafe")
  pair.colors = unlist(pair.pals) # unlist, to allow for easy color-selection
  
  pair.count <- 0 # indexed at 0 bc pairs
  qual.count <- 1 # indexed at 1 bc singles
  groups.colors <- list()
  for ( group in groups.cols  ) {
    groups.vals <- sort(unique( annot[[group]] )) # sort unique annotations (sorting prevents different color assignments based on order-of-appearanace)
    groups.vals[is.na( groups.vals )] <- "NA" # set NA to a string
    # groups.vals[groups.vals=="n.a." | groups.vals=="na"] <- "NA" # set n.a. or na to a NA (this would... probably mess with annot values. leaving alone for now)
    
    # assign paired or qualitative colors, depending on vector length
    if( length( groups.vals ) == 2 || ( length( groups.vals ) == 3 # if we have a pair, or a group of three with an NA val
                                        && "NA" %in% groups.vals ) ){
      pair.idx <- ( pair.count * 2 ) + 1 # choose color-pair index
      colors <- pair.colors[ c( pair.idx, pair.idx + 1 ) ] # assign colors
      pair.count <- ( pair.count + 1 ) %% length(pair.colors) # increment pair.count, or restart if we've exhausted the pairs list
      
      # Assign "normal" annotations the lighter color
      normal_annots = c("nat", "wt", "unmut","normal","0") # possible "normal" annotations (not case sensitive)
      norm_index = unlist(sapply(normal_annots, function(annot) { grep(annot, groups.vals, ignore.case=TRUE, fixed = TRUE) }, USE.NAMES = FALSE)) # identify index of groups.vals that is "normal"
      if ( length(norm_index)!=0 ) { # if any of the "normal" annots can be found
        groups.vals <- c( groups.vals[-norm_index], groups.vals[norm_index]) # move norm value to end. this will assign it the lighter color.
      }
      
    } else { # pick a qualitative colorscheme
      qual.colors = qual.pals[[qual.count]] # choose color palate
      n_colors = ifelse("NA" %in% groups.vals, length(groups.vals)-1, length(groups.vals)) # choose n-1 or n colors, depending on whether there were / were not NA vals
      
      while ( n_colors > length(qual.colors) ) # if there aren't enough colors in this palate
        qual.colors = c(qual.colors, qual.pals[[qual.count+1]]) # keep tacking on new palates until we have enough colors. 
        # NOTE: if someone had a ridiculously large number of annotations, this WILL recycle colors. but dear god I hope someone doesn't try to run an annot with 28+ values....
      
      colors <- eval_with_seed( { sample(qual.colors, n_colors) } , seed_string = group);  # sample colors pseudo-randomly (using group name as seed, to make sure we always choose the "same" random colors)
      qual.count <- qual.count + 1 %% length(qual.pals) # increment qual.count, or restart if we've exhausted the pairs list
    }
    
    # assign NA values grey
    if ( "NA" %in% groups.vals ){ # if there was an NA value
      groups.vals <- c( groups.vals[-which( groups.vals == "NA" )], "NA" ) # move NA to end
      colors <- c( colors, "#BBBBBB" ) # set color to grey
    }
    
    groups.colors[[group]]$vals <- groups.vals
    groups.colors[[group]]$colors <- colors
  }
  return( groups.colors )
}

restore_config_colors <- function( groups.colors.old, typemap.csv ){
  annot <- read_annot() # read in annotations
  
  groups.colors <- list() # initialize groups.colors
  for ( group in groups.cols ) {
    groups.vals <- sort(unique( annot[[group]] )) # get / sort unique values
    groups.vals[is.na( groups.vals )] <- "NA"
    
    # assign default colors to that group // initialize
    groups.colors[[group]] <- set_default_colors( group , typemap.csv)[[group]]
    
    # if the group already existed, and some values already had colors assigned
    if ( (group %in% names(groups.colors.old))  &&  any(groups.vals %in% groups.colors.old[[group]]$vals) ) {
      # override those values 
      ind = match(groups.vals, groups.colors.old[[group]]$vals) # indices of groups.colors.old values that appear in current data-- ordered as they should appear in the current data
      groups.colors[[group]]$colors[ind] = unname(groups.colors.old[[group]]$colors[ind]) # use corresponding colors
    }
  }
  
  return( groups.colors )
}


panda_set_default_colors <- function(display=TRUE) {
  groups.colors <<- set_default_colors( groups.cols, typemap.csv )
  if (display) display_colors( groups.cols, groups.colors )
}

panda_colors_edit_index <- function( groups.colors ){
  for ( group.idx in 1:length( names( groups.colors ) ) ){
    group <- names( groups.colors )[group.idx]
    print( glue( "{sep}\n{add_space( group.idx )}] {group} =\n" ) )
    values <- groups.colors[[group]]$vals
    print( glue( "\t{add_space( 1:length( values ) )}: {values}" ) )
  }
}

change_current_colors <- function( groups.cols, all.groups, groups.colors, byIndex = TRUE ){
  color_exit_message = printX("CANCELLED","Prior changes have been saved", invisible = TRUE)
  new.config <<- TRUE
  groups.flag <- TRUE
  while ( groups.flag ){
    ### Prompt for Group Index
    choice = smart_readline(prompt = "\n$$ Enter group index for color modification: ",
                            custom_condition = function(choice) {!is.na(groups.cols[as.numeric(choice)])},  # check if valid group index
                            custom_warning = printX("ERROR", glue("Invalid group index. Please enter an index in the range {paste(1, 'to', length(groups.colors))}")),
                            exit_message = color_exit_message, exit_function = break)
    change.group <- groups.cols[as.numeric(choice)] # set change.group according to choice

    if (byIndex) { # if we are inputting colors by index
      value.flag = T # initialize value.flag
      while( value.flag ){
        ### Prompt for Value Index
        choice = smart_readline(prompt = "\n$$ Enter value index: ",
                                custom_condition = function (choice) {!is.na(groups.colors[[change.group]]$colors[as.numeric(choice)])}, # check if valid value index, within change.group
                                custom_warning = printX("ERROR", glue("Invalid value index for {change.group}. Please enter an index in the range {paste(1, 'to', length(groups.colors[[change.group]]))}")),
                                exit_message = NULL, exit_function = break)
        change.value <- as.numeric(choice)
        
        ### Prompt for Color Hexcode
        choice = smart_readline(prompt = "\n$$ Enter color hex: ",
                                custom_condition = function(choice) {grepl("^#?[A-Fa-f0-9]{6}$", choice)}, # check if entry is valid HEX code
                                custom_warning = printX("ERROR", "Invalid hex color. Please try again"),
                                exit_message = NULL, exit_function = break)
        change.color <- choice
        # if change.color is a valid hex code, but has no # symbol, add it.
        if ( grepl("^[A-Fa-f0-9]{6}$", change.color ) )  { change.color <- paste0("#",change.color) }
        
        ### Modify Color
        groups.colors[[change.group]]$colors[change.value] <- change.color 
        
        ### Check whether another value in this group should be changed
        value.flag <- y2true( glue("\n$$ Continue editing colors for {change.group}?") )
      }
    } else {
      ### Prompt for Color-Hexcode Vector
      nvals = length(groups.colors[[change.group]]$colors) # number of elements required for groups.colors vector
      rnd_hex_codes <- function(n) { # generates n random hex code. not important. only used to generate example-syntax for the custom_warning below
        hexdigits = "abcdefABCDEF0123456789" %>% strsplit(., "") %>% unlist()
        hexcodes = replicate(n, paste0("#", paste0(sample(hexdigits, 6), collapse="")))
        return(hexcodes)
      } 
      choice = smart_readline(prompt = paste0("\n$$ ", glue("Updating all colors for {change.group}. Values are {paste(groups.colors[[change.group]]$vals, collapse=', ')}."),
                                              "\n$$ ", glue("Enter {nvals} hex colors, separated by commas: " )),
                              custom_condition = function(choice) {
                                choice %>% gsub(" ","", .) %>% # remove whitespace
                                  strsplit(., split=",") %>% .[[1]] %>% # split into vector
                                  lapply(., function(ch) grepl("^#[A-Fa-f0-9]{6}$", ch)) %>% # every input is a valid hex code (must have #)
                                  all(length(.)==nvals, .) },
                              custom_warning = printX("ERROR", glue("Invalid formatting. Please enter {nvals} hex codes, separated by commas (e.g. {paste(rnd_hex_codes(nvals), collapse=', ')})")),
                              exit_message = color_exit_message, exit_function = break)
      change.color <- choice %>% gsub(" ","", .) %>% # remove whitespace
        strsplit(., split=",") %>% .[[1]]
      
      groups.colors[[change.group]]$colors <- change.color
    }
    groups.flag <- y2true("\n$$ Continue modifying defaults for other groups?")
  }
  
  display_colors( groups.cols, groups.colors ) # display colors
  return( groups.colors )
}

panda_colors_edit <- function( byIndex = TRUE ){
  groups.colors <<- change_current_colors(
    groups.cols, all.groups, groups.colors, byIndex = byIndex )
  print( DONE )
}

### ===
### Section. COSMO labels
### ===

# initialize cosmo parameters in case user never runs COSMO cell
cosmo.params <- list(run_cosmo = FALSE)

# function to get user input for cosmo
select_COSMO_attributes <- function() { 
  
  # initialize again in case user runs cell multiple times
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
    "{paste0( sort(unique( annot[[sample.sets$all$fil_col]] )), collapse = ';' )}" )
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
    
    fil_col <- smart_readline(prompt = "$$ Filter samples based on this column. Enter name: ",
                              custom_condition = function(fil_col) { ( fil_col %in% all.groups ) }, # check if input is a valid column
                              custom_warning = printX("ERROR", glue("Invalid column entry, please try Again"), invisible = TRUE),
                              exit_message = NULL,
                              exit_function = ifelse(y2true("\n$$ Continue adding more sets?"),
                                                     {flush.console(); next},
                                                     {flush.console(); printX("QUITTING", "Prior changes have been saved"); break}) )

    # fil_col <- as.character( readline(
    #   prompt = "$$ Filter samples based on this column. Enter name: " ) )
    # if ( !( fil_col %in% all.groups ) )
    #   stop( "\n\n{sep}\n.. ERROR. Invalid column entry. Try Again." )
    values <- sort(unique( annot[[fil_col]] ))

    select.value.msg <- "$$ Enter one or more values from above, separated by semicolons"
    prompt.stmt <- glue( "{paste0( values, collapse = ', ' )}" )
    prompt.stmt <- glue( ".. Values: {prompt.stmt}\n{select.value.msg}" )
    prompt.stmt <- glue( "{prompt.stmt} (e.g. {sample(values,1)}")
    prompt.stmt <- paste0(prompt.stmt, ifelse( length(values)>1 ,
                                               paste0( ', or ', paste0(sample(values,2), collapse=';'), ")."),
                                               ').')  )
    
    
    fil_val <- smart_readline(prompt = prompt.stmt,
                              custom_condition = function(fil_val) { unlist (strsplit (fil_val, split=';')) %in% values %>% all() }, # check if all fil_val are in values
                              custom_warning = printX("ERROR", glue("Value entered is incorrect, please try again"), invisible = TRUE),
                              exit_commands = c("quit","QUIT","exit","EXIT","exit","EXIT"), exit_message = NULL, # only
                              exit_function = ifelse(y2true("\n$$ Continue adding more sets?"),
                                                     {flush.console(); next},
                                                     {flush.console(); printX("QUITTING", "Prior changes have been saved"); break}) )
    fil_val <- unlist (strsplit (fil_val, split=';'))
    
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
  if (!exists ("normalize.prot") || !exists ("filter.prot") || !exists ("run.ptmsea") ) 
    stop( glue( "\n\n{sep}\n.. ERROR. Preprocessing has not been selected. Run panda_preprocessing()." ) )
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
  lines$groups.cols.continuous <- groups.cols.continuous
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
  lines$filter.proteomics <- filter.prot
  lines$run.ptmsea <- run.ptmsea # needs to be integrated ?
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
  catcmd <- glue ("cat {cfg} {original.params} > {params}") # system command to concatenate config.yaml file, with typemap.yml$parameters file
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


