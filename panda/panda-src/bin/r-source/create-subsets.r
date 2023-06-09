#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
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
      c( "-i", "--inputs.dir"  ),
      dest   = "inputs.dir",
      action = "store",
      type   = 'character'),
    
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


read_sets <- function( location ){
  
  meta.file <- read.csv(glue( "{location}/sample_set_membership.tsv" ), sep="\t") # sample-set membership file
  set.names = unique(meta.file[[1]]) # get list of sets from the first column
  
  #  convert the meta.file into a list
  sets_list = sapply( set.names, function(set) {
    set.members = meta.file[meta.file[[1]]==set,2] #for each set, return a vector listing the samples in that set
    return(make.names(set.members)) # return make.names() version of Sample.IDs
  } )
  
  return( sets_list )
}

load_files <- function( gct.types, csv.types, inputs.dir, wkspace ){
  #initialize lists, to house the three file-types
  data <- list(gcts=list(),
               csvs=list(),
               groups=NULL)
  
  for (type in gct.types)
    data[['gcts']][[type]] <- parse_gctx(glue("{inputs.dir}/{wkspace}-{type}.gct")) # read in GCT object
  for (type in csv.types)
    data[['csvs']][[type]] <- read.csv(glue("{inputs.dir}/{wkspace}-{type}.csv")) # read in CSV object
  
  return(data)
}


set_subset_gct <- function( gct, type, set.name, set.members, set.path ){
  
  gct.subset <- cmapR::subset_gct(gct, cid = set.members) # subset GCT to set.members
  
  out.file <- glue( "{set.path}/{type}-subset.gct" )
  write.gct( gct.subset, out.file, ver = 3, appenddim = FALSE )
}

set_subset_csv <- function( csv, type, set.name, set.members, set.path ){

  csv.subset = csv[csv$Sample.ID %in% set.members, ] # subset full CSV file, based on when Sample.ID is in set.members
  

  out.file <- glue( "{set.path}/{type}-subset.csv" )
  write.csv( csv.subset, out.file, row.names = F, quote = T, na = "" )
}


subset_all_files <- function( gct.types, csv.types, location, inputs.dir, wkspace ){
  sets <- read_sets( location ) # read in sets, using the sample_set_membership.tsv file in opt$location
  data <- load_files( gct.types, csv.types, inputs.dir, wkspace ) # read in GCT/CSV/group data
  
  for ( set.name in names( sets ) ){ # for each set
    set.members = sets[[set.name]]
    
    # set up directory
    set.path <- glue( "{location}/subsets/{set.name}" )
    dir.create( set.path, showWarnings = T )
    
    # subset GCT files
    invisible( lapply( gct.types, function( type ) # for each GCT object
      set_subset_gct( data[['gcts']][[type]], type, set.name, set.members, set.path ) ) ) # subset to the sample_set and write file
    
    # subset CSV files
    invisible( lapply( csv.types, function( type )
      set_subset_csv( data[['csvs']][[type]], type, set.name, set.members, set.path ) ) ) # subset to the sample_set and write file
  }
}

main <- function(){
  set_arguments()
  gct.types <- unlist( strsplit( opt$gct.types, split = ';' ) )
  csv.types <- unlist( strsplit( opt$csv.types, split = ';' ) )
  subset_all_files( gct.types=gct.types,
                    csv.types=csv.types,
                    location=opt$location,
                    inputs.dir=opt$inputs.dir,
                    wkspace=opt$wkspace)
}

if ( !interactive() ){
  main()
}
