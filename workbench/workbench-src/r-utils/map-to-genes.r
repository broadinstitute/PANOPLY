

# $Id$
print ("$Id$")


Rutil.path <- switch (Sys.info()[['sysname']],
                      Windows = {'//flynn-cifs/prot_proteomics/Projects/R-utilities'},
                      Darwin = {'/Volumes/prot_proteomics/Projects/R-utilities'},
                      Linux = {ifelse (file.exists('/prot/proteomics/Projects/R-utilities'),
                                       '/prot/proteomics/Projects/R-utilities',   # on Broad systems
                                       '~/proteomics/R-utilities')})              # elsewhere

Source <- function (f) {
  # adaptive source file location -- from current directory or manidr's R-utiliites
  if (file.exists (f)) {source (f)}
  else {
    f.os <- file.path (Rutil.path, f)
    if (file.exists (f.os)) source (f.os)
    else stop (paste ("Can't find file", f.os, '-- mount prot_proteomics?'))
  }
}


Source ('io.r')
Source ('misc.r')
## NB: needs "repeated" function from misc.r



# Map Gene/Protein IDs, between keytypes available with AnnotationDbi (uses org.Hs.eg.db (HUMAN) database)
map_id = function(id_col, keytype_from, keytype_to="REFSEQ", remove_duplicate_keys = TRUE, ome_lookup_key="", silent=FALSE) {
  if (!remove_duplicate_keys) {
    warning("haven't added support for keeping duplicates-- this argument doesn't do anything! sorry!")
  }
  
  if (keytype_from==keytype_to){
    if (!silent) print("No ID conversion necessary.")
    return(id_col) 
  } else {
    # load in main AnnotationDbi Database
    require(org.Hs.eg.db, include.only = 'org.Hs.eg.db')
    db = org.Hs.eg.db
    
    # Check if keytypes are valid ID options
    keytype_options <- AnnotationDbi::keytypes(db) #show key options
    if (!(keytype_from %in% keytype_options)) {
      stop(paste0("The specified input keytype, ", keytype_from, ", is not a valid keytype. Run AnnotationDbi::keytypes(", db$packageName, ") to get a list of valid keytypes."))
    } else if (!(keytype_to %in% keytype_options)) {
      stop(paste0("The specified output keytype, ", keytype_to, ", is not a valid keytype. Run AnnotationDbi::keytypes(", db$packageName, ") to get a list of valid keytypes."))
    }
    
    # Get list of unique IDs (omitting NA, NULL, or "" values)
    id_col_noMissing = na.omit(id_col %>% na_if(c("", NULL)))
    unique_ids = unique(id_col_noMissing)
    
    # if ENSEMBLPROT is involved, use the EnsDb.Hsapiens.v79 database
    # convert to ENTREZID as an intermediate step, for a higher ID-conversion rate
    if (keytype_from=="ENSEMBLPROT" | keytype_to=="ENSEMBLPROT") {
      require(EnsDb.Hsapiens.v79, include.only = 'EnsDb.Hsapiens.v79') #load relevant DB
      if (!silent) print(paste0("To get a higher ID-match rate for ENSEMBLPROT IDs, ", keytype_from,
                                " is first being mapped to ENTREZID, before being mapped to ",
                                keytype_to, "; uses the EnsDb.Hsapiens.v79 database."))
      
      ### Initial Conversion to ENTREZID; overwrite keytype_from
      db_tmp = db # use default gene-id-database for initial conversion to ENTREZID
      if (keytype_from=="ENSEMBLPROT") { # if we're converting ENSEMBLPROT->ENTREZID
        db_tmp = EnsDb.Hsapiens.v79; keytype_from="PROTEINID" } # use EnsDb.Hsapiens.v79 for initial conversion to ENTREZID
      unique_ids = AnnotationDbi::mapIds(db_tmp, keys=unique_ids, column=c("ENTREZID"), keytype=keytype_from)
      
      # print how many IDs were converted successfully
      percent_converted = sum(!is.na(unique_ids)) / sum(!is.na(unique(id_col_noMissing)))
      if (!silent) print(paste0(round(percent_converted*100,1),"% of IDs were successfully converted from ", keytype_from, " to ENTREZID."))
      
      # replace NA w/ unique-placeholder-value, to avoid mismatched-length-errors
      unique_ids[which(is.na(unique_ids))] = paste0("nomatch.", which(is.na(unique_ids))) 
      
      # update keytypes for downstream use
      keytype_from="ENTREZID" #set keytype_from ID to ENTREZID
      if (keytype_to=="ENSEMBLPROT") { # if we're converting ENTREZID->ENSEMBLEPROT
        db = EnsDb.Hsapiens.v79; keytype_to="PROTEINID" } # use EnsDb.Hsapiens.v79 for downstream conversions
      
      # sanity check to ensure that original IDs still match up
      # if (sum(names(unique_ids)!=unique(na.omit(id_col)))>0) {
      #   stop("Something went wrong with the ID conversion-- IDs no longer match up!")
      # }
    }
    
    # Create key to match keytype_old->keytype_new, pulling the FIRST ID that contains ome_lookup_key
    key_vec=AnnotationDbi::mapIds(db, keys=as.character(unique_ids), column=keytype_to, keytype=keytype_from, multiVals = 'list') %>%
      lapply(., function(e) {grep(ome_lookup_key, e, value=TRUE) %>% .[1]}) # grab FIRST ID that contains ome_lookup_key
    # put data into dataframe format, for use with left_join()
    key = data.frame(old_id = unique(id_col_noMissing), #use original ID as old_id (omitting missing values)
                     new_id = unlist(key_vec)) #use final ID as new_id
    
    # print how many IDs were converted successfully
    percent_converted = sum(!is.na(key$new_id)) / sum(!is.na(key$old_id))
    if (!silent) print(paste0(round(percent_converted*100,1),"% of IDs were successfully converted from ", keytype_from, " to ", keytype_to, "."))
    
    # Use key to map old id_col to new IDs (id_col_converted)
    tmp = data.frame(old_id = id_col) # tmp dataframe
    id_col_converted <- left_join(tmp, key, by='old_id') %>% #match old IDs to new IDs
      .[,'new_id'] #pull column with new IDs
    
    return(as.character(id_col_converted)) #return as characters so OTHER functions don't decide this is FACTORS for some reason
  }
}


map.to.unique.genes <- function (input.file, output.file, duplicate.gene.policy='median',
                                 gene.map=file.path (Rutil.path, '..', 'gene-mapping', 'RefSeq.20170701-RefSeq-GeneName-map.txt'),
                                 map.id.col='refseq_protein', map.genename.col='gene_name', official.genenames=FALSE,
                                 official.genenames.file=file.path (Rutil.path, '..', 'gene-mapping', 'std-gene-map', 'gene-symbol-map.csv')) {
  # map protein ids to gene symbols using the default RefSeq -> gene symbol map
  # if required convert gene symbols to "offical" gene symbols (disabled by default)
  # also eliminate duplicate gene ids by combining all rows with the same gene
  # this function is a convenient implementation of other basic functions in this file
  
  d <- read.gct2 (input.file, check.names=FALSE)
  d.info <- d[, 1:2]
  d.data <- d[, 3:ncol(d)]
  
  d.genes <- map.proteins.to.genes (cbind (d.info, d.data), id.col=1, gi.gene.map.file=gene.map,
                                    protein.id.col=map.id.col, gene.id.col=map.genename.col,
                                    keep.cols=c(map.id.col, colnames (d.data)))
  if (nrow (d.genes) == 0)
    # gene mapping failed -- try using a different id.col
    # (in some cases, the protein id is the second column)
    d.genes <- map.proteins.to.genes (cbind (d.info, d.data), id.col=2, gi.gene.map.file=gene.map,
                                      protein.id.col=map.id.col, gene.id.col=map.genename.col,
                                      keep.cols=c(map.id.col, colnames (d.data)))
  
  if (nrow (d.genes) == 0) 
    # nothing worked -- give up
    stop ('Gene mapping failed ... giving up.')
  
  d.output <- process.duplicate.genes (d.genes, genesym.col=1, data.cols=3:ncol(d.genes),
                                         official.symbols.file=official.genenames.file, map.genes=official.genenames,
                                         policy=duplicate.gene.policy)
  
  write.gct (d.output, output.file)
}

combine.duplicate.genes <- function (input.file, output.file, policy='median') {
  # roll-up duplicate genes in input.file (gct format, gene symbols in first column) 
  # using specified policy and write resulting dataset to output.file (gct format)
  # if input data has duplicated genes, combine rows with appropriate "policy":
  #   maxvar: select row with largest variance
  #   union: union of binary (0/1) values in all rows (e.g for mutation status)
  #   median: median of values in all rows (for each column/sample)
  #   mean: mean of values in all rows (for each column/sample)
  #   min: minimum of values in all rows (for each column/sample)
  
  data <- read.gct2 (input.file, check.names=FALSE)
  data.rollup <- process.duplicate.genes (data, genesym.col=1, data.cols=3:ncol(data),
                                          official.symbols.file=NULL, map.genes=FALSE, policy=policy)
  write.gct (data.rollup, output.file)
}



process.duplicate.genes <- function (data, genesym.col, data.cols, official.symbols.file, 
                                     map.genes=TRUE, policy="maxvar") {
  ## map gene symbols to official gene symbols in map.file
  ##   map.file must have 2 columns -- symbol and official.symbol
  ##   data has (unofficial) gene symbols in genesym.col, with data in data.cols
  ##   for gct, genesym.col=1 (Name) and data.cols=3:ncol(data)
  ##   for netgestalt, genesym.col=1 (GeneSymbol) and data.cols=2:ncol(data)
  ##   n.b: data.cols must the last set of columns
  ## if data has duplicated genes, combine rows with appropriate "policy"
  ##   maxvar: select row with largest variance
  ##   union: union of binary (0/1) values in all rows (e.g for mutation status)
  ##   median: median of values in all rows (for each column/sample)
  ##   mean: mean of values in all rows (for each column/sample)
  ##   min: minimum of values in all rows (for each column/sample)
  
  
  map.genesym <- function (d, map.file=official.symbols.file) {
    # gene map for mapping to official symbols created in create-gene-map.r
    gene.map <- read.csv (map.file, as.is=TRUE)  
    
    # keep only those genes that map to an official symbol, and perform mapping
    d <- d [ d[,genesym.col] %in% gene.map[,'symbol'], ]
    d[,genesym.col] <- unlist (lapply (d[,genesym.col],
                                        function (x) gene.map [ which (gene.map[,'symbol'] %in% x), 'official.symbol']))
    invisible (d)
  }
  
  
  data [,genesym.col] <- toupper (as.character (data [,genesym.col]))  # convert all gene symbols to uppercase
  
  # keep only those genes that map to an official symbol
  if (map.genes) data <- map.genesym (data)
  
  # deal with duplicate symbols (if any)
  repeated.index <- repeated (data[,genesym.col])
  data.unique <- data [ !repeated.index, ]
  
  data.rep <- data [ repeated.index, ]
  data.rep.chosen <- NULL
  if (nrow (data.rep) > 0) {
    rep.genes <- data.rep [,genesym.col]
    for (g in (unique (rep.genes))) {
      data.subset <- data.rep [g==rep.genes,]
      if (policy == "union") {
        # union of "mutations" in all rows
        all.cols <- colnames (data)
        if (is.numeric (data.cols)) all.cols <- 1:ncol(data)
        row <- c ( sapply (data.subset[1, setdiff (all.cols, data.cols)], as.character), 
                   apply (data.subset[, data.cols], 2, function (x) ifelse (sum(x)>0, 1, 0)))
        data.rep.chosen <- rbind (data.rep.chosen, row)
      }
      
      if (policy == "median") {
        # median of row values
        all.cols <- colnames (data)
        if (is.numeric (data.cols)) all.cols <- 1:ncol(data)
        row <- c ( sapply (data.subset[1, setdiff (all.cols, data.cols)], as.character), 
                   apply (data.subset[, data.cols], 2, median, na.rm=TRUE) )
        data.rep.chosen <- rbind (data.rep.chosen, row)
      }

      if (policy == "mean") {
        # mean of row values
        all.cols <- colnames (data)
        if (is.numeric (data.cols)) all.cols <- 1:ncol(data)
        row <- c ( sapply (data.subset[1, setdiff (all.cols, data.cols)], as.character), 
                   apply (data.subset[, data.cols], 2, mean, na.rm=TRUE) )
        data.rep.chosen <- rbind (data.rep.chosen, row)
      }
      
      if (policy == "maxvar") {
        # select row with largest variance
        choose <- which.max (apply (data.subset[, data.cols], 1, sd, na.rm=TRUE))
        data.rep.chosen <- rbind (data.rep.chosen, data.subset [choose, ])
      }
      
      if (policy == "min") {
        # minimum of row values -- useful for p-values and the like
        all.cols <- colnames (data)
        if (is.numeric (data.cols)) all.cols <- 1:ncol(data)
        row <- c ( sapply (data.subset[1, setdiff (all.cols, data.cols)], as.character), 
                   apply (data.subset[, data.cols], 2, min, na.rm=TRUE) )
        data.rep.chosen <- rbind (data.rep.chosen, row)
      }
    }
    colnames (data.rep.chosen) <- colnames (data.rep)
  }
  
  data.final <- rbind (data.unique, data.rep.chosen)
  for (i in data.cols)
    if ( class (data.final[,i]) != "numeric" )
      data.final[,i] <- as.numeric (as.character (data.final[,i]))

  invisible (data.final)
}



map.proteins.to.genes <- function (data, id.col, gi.gene.map.file,
                                   protein.id.col='refseq_protein', gene.id.col='gene_name', keep.cols=NULL) {
  ## maps protein ids in data (in id.col) to gene ids
  ## gi.gene.map.file (tab delimited) must contain both protein.id.col and gene.id.col (may have extra columns)
  gi.gene.map <- read.delim (gi.gene.map.file, sep='\t')
  colnames (data)[id.col] <- protein.id.col
  mapped.data <- merge (data, gi.gene.map, by=protein.id.col)
  if (is.null (keep.cols)) data.final <- mapped.data [, c (gene.id.col, protein.id.col, colnames(data)[-id.col])]
  else data.final <- mapped.data [, c (gene.id.col, keep.cols)]
  
  invisible (data.final)  
}
