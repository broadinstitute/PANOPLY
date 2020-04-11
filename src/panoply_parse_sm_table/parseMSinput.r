

## Code for preprocessing TMT/iTRAQ data
## Converts Spectrum Mill output to preprocessed data files
## Additional files are needed to run this code.


source ('config.r')


process.dataset <- function (dataset, out.prefix, id.col, proteome=FALSE, exclude=NULL,
                             additional.cols=NULL, species.filter=TRUE, expt.design=NULL) {
  # the input dataset is an ssv file that is output from Spectrum Mill
  # that has a 2-line header -- the header is parsed to get sample and run # info
  #
  # exclude is a vector of TRUE/FALSE values that indicate which data columns (in order, ala cls)
  # should be excluded from the data processing
  # additional.cols specifies extra columns to keep in the output table (and will result
  # in a gct3 file, since gct2 cannot support >2 info columns)
  # species.filter if TRUE will discard all proteins/peptides not from homo sapiens
  #
  # expt.design file should contain at least Sample.ID, Experiment and Channel columns
  # (Sample.IDs MUST be unique, valid R names; duplicate samples should have the same sample names, 
  #  but include a replicate.indicator, eg .REP, followed by a unique suffix: <name>.REP1)
  # additional columns in the file will be treated as sample annotations; could be
  # derived from Spectrum Mill reporter_sample_template (converted to long from)
  
  
  parse.info <- function (str) {
    # extract sample id and run number from the str
    # used only when the experiment design file is not specified
    if (project.name == 'cptac2.tcga') {
      # specific to the CPTAC2 TCGA (retrospective) sample analysis
      pieces <- strsplit (toString (str), split='/')[[1]]
      parts <- strsplit (pieces[length(pieces)], split='_')[[1]]
      run <- parts[1]
      samples <- parts[2:(n.channel-1)]
      
      samples.run <- paste (samples, run, sep='.')
      return (samples.run)
    }
    
    # for other project.names, write appropriate sample/run extraction code
    return (str)
  }

  get.sample.names <- function (run, exptdesign, channels=plex.channels) {
    # get sample ids for experiment specified by run, using experiment design file
    subset <- exptdesign[ exptdesign[,'Experiment']==run, ]
    rownames (subset) <- as.character (subset[,'Channel'])
    return (as.vector (unlist (subset [channels,'Sample.ID'])))
  }
  
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)    # remove leading and trailing blanks
  
  process.header <- function (dataset) {
    # collect info from the first 2 rows of a Spectrum Mill ssv file
    # and generate cols/colnames for ratios
    header <- unlist (lapply (scan (dataset, what="char", nlines=1, sep=';'), trim))
    if (!proteome) {
      header.1 <- header
      header.2 <- unlist (lapply (scan (dataset, what="char", nlines=1, skip=1, sep=';'), trim))
      header.2 <- header.2 [ sapply (header.2, function (x) nchar (x) > 0) ]
    } else {
      # protein tables have different formatting and an empty column at the start
      header.1 <- unlist (lapply (header,
                                  function (x) {
                                      s <- strsplit (x, split=' ')[[1]]
                                      retval <-  s [length(s)]	 # keep the string after the first space
                                      if (length (s) < 2 && grepl (header.pat, s[1]))
                                          retval <- paste (s[1], 'ignore', sep='')
                                      return (retval)
                                  }))
      header.2 <- unlist (lapply (header,
                                  function (x) {
                                      s <- strsplit (x, split=' ')[[1]]
                                      return ( s[1] ) 	 # keep the string before the first space
                                  }))
    }
    n <- n.channels - 1
    ratio.fields <- grep (ratio.pat, header.1 [1:length(header.2)])
    stddev.fields <- grep (stddev.pat, header.1 [1:length(header.2)])
    intensity.fields <- grep (intensity.pat, header.1 [1:length(header.2)])
    numratio.fields <- grep (numratio.pat, header.1[1:length(header.2)])
    numspectra.fields <- rep (grep (numspectra.pat, header.1[1:length(header.2)]), each=n)
    totalint.fields <- rep (grep (totalint.pat, header.1[1:length(header.2)]), each=n)
    unique_pep.fields <- rep (grep (unique_pep.pat, header.1[1:length(header.2)]), each=n)
    refint.fields <- rep (grep (refint.pat, header.1[1:length(header.2)]), each=n)
    
    first.in.run <- ratio.fields [ seq (1, length(ratio.fields), n) ]  # X-plex labels with (X-1) ratios
    if (is.null (expt.design)) {
      col.names <- unlist (lapply (header.2 [first.in.run], parse.info))
    } else {
      ed <- read.csv (expt.design)
      col.names <- unlist (lapply (1:length(first.in.run), get.sample.names, ed))
    }
        
    if (project.name == 'default') {
      # use this section to perform alternative processing based on sample.source
    }
    
    return (list (col.numbers=list (ratio.fields=ratio.fields,
                                    stddev.fields=stddev.fields,
                                    intensity.fields=intensity.fields,
                                    numratio.fields=numratio.fields, 
                                    numspectra.fields=numspectra.fields,
                                    totalint.fields=totalint.fields,
                                    unique_pep.fields=unique_pep.fields,
                                    refint.fields=refint.fields),
                  cols.list=c('ratio', 'stddev', 'intensity', 'num-ratio', 'num-spectra',
                              'precursor-intensity', 'unique-peptides', 'reference-intensity'),
                  col.names=col.names, col.names.all=header.1))
  }

  
  write.out <- function (info, data, file) {
    # if there are no additional.cols or sample annotations, write standard gct(2); 
    # else output gct3 format
    sample.annotations <- data.frame ()
    if (!is.null (expt.design)) {
      # check if there are any sample annotations
      d <- read.csv (expt.design)
      annot.cols <- setdiff (colnames (d), 'Sample.ID')
      sample.annotations <- d[, c('Sample.ID', annot.cols)]
      # add sample QC info (pass/fail)
      if (! qc.col %in% sample.annotations) {
        if (!is.null(sampleQC.cls)) qc <- read.cls (sampleQC.cls)
        else qc <- rep (qc.pass.label, nrow(sample.annotations))
        # add to annotation table
        sample.annotations <- cbind (sample.annotations, qc)
        colnames (sample.annotations)[ncol(sample.annotations)] <- qc.col
      }
    }
    # set up GCT v1.3 object
    gct <- new ('GCT')
    gct@mat <- as.matrix (data)
    gct@rid <- as.character (info[,1])
    gct@cid <- as.character (colnames (data))
    gct@rdesc <- fix.datatypes (info)
    gct@cdesc <- fix.datatypes (sample.annotations)
    gct@version <- "#1.3"
    gct@src <- file
    # write output file
    write.gct (gct, file, ver=3, precision=ndigits, appenddim=FALSE)
  }


  species.col <- function (col.names.all) {
    # species is indicated in the species or species2 column
    # species2 supercedes species
    col.num <- grep ('^species2$', col.names.all)
    if (length(col.num) == 0) col.num <- grep ('^species$', col.names.all)
    if (length(col.num) < 1 || length(col.num) > 1) stop ('Having trouble finding the species column ...')

    return (col.num)
  }
  
  
  d <- read.delim (dataset, sep=';', skip=ifelse(proteome,1,2), header=FALSE)
  if (proteome) d <- d[,-1]   # remove first column which is blank
  header.info <- process.header (dataset)
  colnames (d) <- header.info$col.names.all

  # filter to remove non-human entries (rows)
  if (species.filter) {
    species <- sapply (d[, species.col(header.info$col.names.all)],
                       function (x) tolower (gsub ("\\s", "", x)))
    # keep <- species %in% c ('human', 'homosapiens')
    # change to below to support ORFs or other "species" names that include 'human' as a substring
    keep <- sapply (species, function (x) grepl ('human', x) || grepl ('homosapiens',x))
    d <- d [keep, ]
  }
    
  
  # info columns (Name and Description for gct)
  # if additional.cols are requested, include those that are present in the input
  info.col1 <- as.character (d[,id.col])
  info.col2 <- as.character (d[,'entry_name'])
  if (is.null (additional.cols)) info.data <- cbind (Name=info.col1, Description=info.col2)
  else info.data <- cbind (id=info.col1, id.description=info.col2, 
                           d[, intersect (additional.cols, colnames(d))])
  
  for (i in 1:length (header.info$cols.list)) {
    fields <- header.info$col.numbers[[i]]
    if (length(fields) == 0) next
    
    if (! is.null (exclude)) {
      # exclude requested columns from the rest of the data processing
      keep <- !exclude
      fields <- fields [keep]
      header.info$col.names <- header.info$col.names [keep]
    }
    
    data <- d [ , fields]
    colnames (data) <-  header.info$col.names
    
    # write out raw data file with no normalization and no filtering    
    write.out (info.data, data, paste (out.prefix, '-', header.info$cols.list[i], '.gct', sep=''))
  }
  
}


###
### generate gct files for analysis
###
if (type == "proteome") {
  info.cols <- c ('geneSymbol', 
                  'numColumnsProteinObserved', 
                  'numSpectraProteinObserved',
                  'protein_mw', 
                  'percentCoverage', 
                  'numPepsUnique',
                  'scoreUnique', 
                  'numPepsUniqueSubgroupSpecificCI',
                  'scoreUniqueSubgroupSpecificCI',
                  'species', 
                  'orfCategory',
                  'accession_number', 
                  'accession_numbers',
                  'subgroupNum', 
                  'entry_name',
                  'scoreUniqueSubgroupSpecificCI',
                  'numPepsUniqueSubgroupSpecificCI')
  process.dataset (dataset = input.data.file,
                   out.prefix = 'proteome', proteome=TRUE, 
                   id.col = 'accession_number', additional.cols=info.cols,
                   expt.design=expt.design.file)
} else { 
  # phosphoproteome, acetylome, or other ptm-ome
  info.cols <- c ('geneSymbol', 
                  'numColumnsVMsiteObserved',
                  'bestScore',
                  'bestDeltaForwardReverseScore',
                  'Best_scoreVML',
                  'Best_numActualVMSites_sty',
                  'Best_numLocalizedVMsites_sty',
                  'Best_numAmbiguousVMsites_sty',
                  'StartAA',
                  'VMsiteFlanks',
                  'Best_numActualVMSites_k',
                  'Best_numLocalizedVMsites_k',
                  'Best_numAmbiguousVMsites_k',
                  'variableSites',
                  'sequence',
                  'sequenceVML',
                  'accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA',
                  'protein_mw',
                  'species',
                  'speciesMulti',
                  'orfCategory',
                  'accession_number',
                  'accession_numbers',
                  'protein_group_num',
                  'entry_name')
  process.dataset (dataset = input.data.file,
                   out.prefix = type, proteome=FALSE, additional.cols=info.cols,
                   id.col = 'accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA',
                   expt.design=expt.design.file)
}

