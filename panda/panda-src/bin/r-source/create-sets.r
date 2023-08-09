#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load( optparse )
p_load( cmapR )
p_load( glue )

set_arguments <- function(){
  optionList <- list(
    make_option(
      c( "-s", "--set" ),
      dest = "set",
      action = "store",
      type = 'character' ),

    make_option(
      c( "-l", "--location" ),
      dest = "location",
      action = "store",
      type = 'character' ),

    make_option(
      c( "-w", "--wkspace" ),
      dest = "wkspace",
      action = "store",
      type = 'character' ) )

  opt <<- parse_args( OptionParser( option_list=optionList ) )
}

set_arguments()
set_subs <- unlist( strsplit( opt$set, split = ',' ) )
set_name <- set_subs[1]
set_fcol <- set_subs[2]
set_fval <- set_subs[3]

oa.print <- ""
other_attributes <- set_subs[-c( 1, 2, 3 )]
for ( attr in other_attributes ) {
  pair <- unlist( strsplit( attr, split = '@' ) )
  oa.print <- glue( "{oa.print}{pair[1]}=\"{pair[2]}\"\n\n" )
}

if ( !( identical( oa.print, "" ) ) ) {
  oa.file.name <- glue( "{opt$location}/other-attributes/{set_name}-other.txt" )
  oa.file <- file( oa.file.name )
  writeLines( oa.print, oa.file )
  close( oa.file )
}

## annotation file
annot.file <- glue( "{opt$location}/pipeline-input" )
annot.file <- glue( "{annot.file}/{opt$wkspace}-annotation.csv" )
annot.data <- read.csv( annot.file, header = T, stringsAsFactors = F )

s.ids <- annot.data$Sample.ID
# check that Sample.IDs are unique-- throw error if not
if( length(s.ids)!=length(unique(s.ids)) ) { stop("Sample.IDs are not unique in Sample Annotations file. Please ensure that Sample.IDs are unique.") }
set_fval <- unlist( strsplit( set_fval, split = ';' ) )
set.list  <- c()
for ( value in set_fval ){
  set.list <- c( set.list,  s.ids[which( annot.data[, set_fcol] == value )] )
}

set.list <- data.frame( paste( set.list, collapse = '\n' ) )
write.table( set.list,
             glue( "{opt$location}/sample-set-members/{set_name}.txt" ),
             col.names = F,
             row.names = F,
             quote = F )
