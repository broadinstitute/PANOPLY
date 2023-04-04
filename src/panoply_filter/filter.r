#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source ('config.r')
Source ('gct-io.r')
Source ('map-to-genes.r')
library(dplyr)  # transition to dplyr to make data wranging more consistent and bug-free

# Locate or Create (gene.id.col) and (qc.col)
# Separate QC types into distinct files
preprocess.dataset <- function (ds, out.prefix=NULL, protein.id.type="REFSEQ", 
                                separate.qc.types=FALSE) {
  
  # if no output prefix has been set, use original file prefix
  if (is.null (out.prefix)) out.prefix <- file.prefix
  
  ### Locate (or Create) gene.id.col in ds@rdesc
  if ( any (grepl (gene.id.col, colnames (ds@rdesc), ignore.case=TRUE)) ) {
    # gene.id.col is already present as an annotation column
    genesym.col <- grep (gene.id.col, colnames (ds@rdesc), ignore.case=TRUE)
    if (length (genesym.col) > 1) {
      genesym.col <- genesym.col[1]
      warning ( paste ('Identified multiple gene symbol columns. Using', genesym.col) )
    }
    ds@rdesc [,gene.id.col] <- as.character (ds@rdesc [,genesym.col])
  } else if ( any (grepl ('gene[.-_]?(name|id|symbol)s?$', colnames (ds@rdesc), ignore.case=TRUE)) ) {
    # gene symbol is already present as an annotation column, with a generic name
    genesym.col <- grep ('gene[.-_]?(name|id|symbol)s?$', colnames (ds@rdesc), ignore.case=TRUE)
    if (length (genesym.col) > 1) {
      genesym.col <- genesym.col[1]
      warning ( paste ('Identified multiple gene symbol columns. Using', genesym.col) )
    }
    ds@rdesc [,gene.id.col] <- as.character (ds@rdesc [,genesym.col])
  } else {
    # gene sumbol is NOT present; map protein id to gene symbols
    ds@rdesc [,gene.id.col] <- sub ("(.*?)\\..?", "\\1", ds@rdesc[,'id']) %>% # prune suffix
      map_id(., keytype_from=protein.id.type, keytype_to="SYMBOL")  %>%  # map to SYMBOLS; map_id() sourced from map-to-genes.R
      as.character(.) # convert to character datatype, if its not already
  }  
  
  
  ### Locate (or Create) qc.col in ds@cdesc
  # Check if qc.col column exists
  if ( is.null(ds@cdesc[[qc.col]]) ) { # if qc.col doesn't exists in ds@cdesc
    warning("No QC column detected. Assuming all samples passed QC.")
    ds@cdesc[[qc.col]] <- rep (qc.pass.label, ncol (ds@mat)) # assign qc.pass.label to every sample
  }
  # Separate QC.Pass samples 
  qc.types = unique ( ds@cdesc[[qc.col]] ) 
  if ( !(qc.pass.label %in% qc.types) ) { # if NO samples have the qc.pass.label, throw an error
    stop(paste0("No samples matched qc.pass.label \'",qc.pass.label,"\' in ds@cdesc[[qc.col]] \'",qc.col,"\'."))
  } else {
    # if separate.qc.types is TRUE, write GCT for each qc.type
    if ( separate.qc.types && length (qc.types) > 1 ) {
      for (qc.type in qc.types) { # for each QC type
        ds.tmp <- col.subset.gct (ds, ds@cdesc[[qc.col]]==qc.type) # subset GCT by column-index; col.subset.gct() sourced from gct-io.r)
        if (qc.type == qc.pass.label) { # QC passes cases
          ds.qc.pass = ds.tmp  # save gct-subset
          file.name <- paste (out.prefix, '.gct', sep='')
        }
        else file.name <- paste (out.prefix, '-', qc.type, '.gct', sep='')
        write.gct (ds.tmp, file.name, ver=3, precision=ndigits, appenddim=FALSE)
      }
    } else { # write GCT for JUST the QC.pass samples (qc.type = qc.pass.label()
      ds.qc.pass <- col.subset.gct (ds, ds@cdesc[[qc.col]]==qc.pass.label) # subset GCT by column-index; col.subset.gct() sourced from gct-io.r)
      file.name <- paste (out.prefix, '.gct', sep='')
      write.gct (ds.qc.pass, file.name, ver=3, precision=ndigits, appenddim=FALSE)
    }
  }
  
  return(ds.qc.pass)
}

filter.dataset <- function (ds, out.prefix=NULL, silent=FALSE,
                            sd.threshold=NULL,
                            combine.replicates=NULL,
                            na.max=NULL,
                            no.na=FALSE) {
  
  # tags = "" #initialize tags to add to output filename (tags added according to applied filters)
  
  # if no output prefix has been set, use original file prefix
  if (is.null (out.prefix)) out.prefix <- file.prefix
  
  
  ### Filters
  
  # filter rows with SD less than sd.threshold across samples
  if (!is.null (sd.threshold)) {
    # tags = paste(tags, "sdfilter", sep="-") # add sdFilt tag to output filename
    rows_before=nrow(ds@mat) # record starting dimensions
    
    ds <- row.subset.gct (ds, index=apply (ds@mat, 1, sd, na.rm=TRUE) > sd.threshold) # subset GCT by row-index; row.subset.gct() sourced from gct-io.r)
    # print filtering results
    if (!silent) print(paste0( as.character(nrow(ds@mat))," rows remaining after SD filtering. ",
                               as.character(rows_before - nrow(ds@mat)),
                               " rows removed for having a standard deviation below ", as.character(sd.threshold), " across all samples"))
    }

  # filter replicates
  if (!is.null (combine.replicates)) {
    # tags = paste(tags, "combRep", sep="-") # add sdFilt tag to output filename # NO TAG NECESSARY
    cols_before=ncol(ds@mat) # record starting dimensions
    
    # combine replicate samples; possible options -- first, mean, median
    # replicates are identified by identical (Participant, Type, and Timepoint) 
    # or (Participant) only if (Type) and (Timepoint) do not exist
    avail.cols <- intersect (colnames (ds@cdesc), c ('Participant', 'Type', 'Timepoint'))
    if (length (avail.cols) > 0 && 'Participant' %in% avail.cols) {
      # identify duplicates
      dup.table <- ds@cdesc %>% mutate (sample.num=row_number()) %>%
        group_by ( across(all_of(avail.cols)) ) %>%  # group by content of avail.cols vector
        mutate(group.num = cur_group_id()) %>%       # group_indices() is deprecated
        filter (n() > 1)
      # for each group of replicates, override the first replicate with the "combined" replicate value 
      for (g in unlist (unique (dup.table [,'group.num']))) {
        rows <- unlist (dup.table [ dup.table[,'group.num']==g, 'sample.num'])
        # replace first occurance of replicate with combined replicate value (use combine.replicates method)
        ds@mat[,rows[1]] <- apply (ds@mat[,rows], 1, function (x) { 
          switch (combine.replicates,
                  "first" = x[which(!is.na(x))][1], # first non-NA value (unless all are NA)
                  "mean" = mean (x, na.rm=T),
                  "median" = median (x, na.rm=T),
                  stop(paste0("Invalid combine.replicates method '",combine.replicates,"'."))) })
      }
      # remove the remaining replicates
      keep.samples <- !duplicated (ds@cdesc[, avail.cols])
      ds <- col.subset.gct (ds, keep.samples)   
    } else { stop("Replicates could not be combined; sample-grouping columns are not present in ds@cdesc.
                  At minimum, the 'Participant' column must be present.") } 
    
    if (!silent) print(paste0( as.character(ncol(ds@mat))," samples remaining after replicate combination via '",combine.replicates,"' method. ",
                               as.character(cols_before - ncol(ds@mat)),
                               " columns removed."))
  }
  
  # filter rows with too many NAs across samples (ideally should happen LAST before write.gct(); other filters impact % missing value)
  if (!is.null (na.max)) {
    # tags = paste(tags, "NArm", sep="-") # add NArm tag to output filename
    
    rows_before=nrow(ds@mat) # record starting dimensions
    
    if (na.max < 1) na.max <- ceiling (ncol(ds@mat) * na.max) # convert na.max to integer, if a fraction 
    # keep rows if observed in at least ( ncol(ds@mat) - na.max ) samples
    ds <- row.subset.gct (ds, index=apply (ds@mat, 1, function (x) sum (!is.finite (x)) < na.max))
    
    # print filtering results
    if (!silent) print(paste0( as.character(nrow(ds@mat))," rows remaining after NA filtering. ",
                               as.character(rows_before - nrow(ds@mat)),
                               " rows removed for having >", round(na.max/ncol(ds@mat)*100), "% missing values across all samples"))
  } 
  
  
  ### Write Filtered Dataset
  write.gct (ds, paste0 (out.prefix, "-filt", ".gct"), ver=3, precision=ndigits, appenddim=FALSE)
  
  
  ### Additional filters (create additional GCT files)
  
  # filter rows with ANY NA readings; creates additional GCT (must happen AFTER na.max filter)
  if (no.na) {
    # tags = gsub("-NArm|$", "-noNA", tags) # replace NArm with noNA (or append noNA)
    rows_before=nrow(ds@mat) # record starting dimensions
    
    # keep rows if observed in ALL samples; exclude any rows with missing values
    ds_noNA <- row.subset.gct (ds, index=apply (ds@mat, 1, function (x) all (is.finite (x))))

    # print filtering results
    if (!silent) print(paste0( as.character(nrow(ds_noNA@mat))," rows remaining after NA filtering. ",
                               as.character(rows_before - nrow(ds_noNA@mat)),
                               " rows removed for having a missing value in at least one samples"))
    # write a noNA GCT file
    write.gct (ds_noNA, paste0 (out.prefix, "-filt-noNA", '.gct'), ver=3, precision=ndigits, appenddim=FALSE)
  }

}





file.prefix = paste (type, '-ratio-norm', sep='')

### Read in Dataset (must be GCT v1.3)
ds <- parse.gctx ( file.path( norm.dir, paste (file.prefix, '.gct', sep='') ) )
if (grepl("1\\.2", ds@version)) {
  stop("GCT v1.2 no longer supported. Please convert to GCT v1.3 or higher.") } # throw error if input is GCT v1.2 (upgrade your file!)


# preprocess dataset, select qc-pass samples
ds_qcpass = preprocess.dataset(ds, out.prefix=file.prefix,
                               protein.id.type=protein.id.type, 
                               separate.qc.types=separate.qc.types) 

if (filter.proteomics) {
  print(paste0(ncol(ds_qcpass@mat), " of ", ncol(ds@mat), " samples are ", qc.pass.label,
               ". Filtering will be applied to these samples."))
  
  ## Change default to No SD filter -- filter removes very few items
  filter.dataset (ds_qcpass, out.prefix=NULL, silent=FALSE,
                  sd.threshold=sd.filter.threshold,
                  combine.replicates=combine.replicates,
                  na.max=na.max,
                  no.na=no.na)
} else {
  print(paste0(ncol(ds_qcpass@mat), " of ", ncol(ds@mat), " samples are ", qc.pass.label,
               ". Filtering will not be applied."))
  file.copy( paste0(file.prefix, '.gct', sep=''), paste0(file.prefix, '-filt', '.gct', sep='') )
}
