
# $Id$
print ("$Id$")


## in preparation for dockerization, use this file to specify all customizable inputs


Rutil.path <- switch (Sys.info()[['sysname']],
                      Windows = {'//flynn-cifs/prot_proteomics/Projects/R-utilities'},
                      Darwin = {'/Volumes/prot_proteomics/Projects/R-utilities'},
                      Linux = {'/prot/proteomics/Projects/R-utilities'})

Source <- function (f) {
  # adaptive source file location -- from current directory or manidr's R-utiliites
  if (file.exists (f)) {source (f)}
  else {
    f.os <- file.path (Rutil.path, f)
    if (file.exists (f.os)) source (f.os)
    else stop (paste ("Can't find file", f.os, '-- mount prot_proteomics?'))
  }
}


Source ('gct-io.r')
Source ('io.r')

# output precision for gct tables
ndigits <- 5    


# determine type (proteome/phosphoproteome) based on command line args; defaults to proteome
# data type defaults to unimodal, unless specified; if unimodal, subset.str is NULL
args <- commandArgs (trailingOnly=TRUE)
type <- ifelse (length(args) > 0, toString (args[1]), "proteome")
data.subset <- ifelse (length(args) > 1, toString (args[2]), "")
subset.str <- ifelse (data.subset=="", data.subset, paste ('-', data.subset, sep='')) 



# directories with raw, pre-processed and normalized data
data.dir <- '../data'
pre.dir <- '../data-with-intensity'
norm.dir <- '../2comp-normalization'
mrna.dir <- '../mrna'
rnaseq.dir <- '../rna-seq'
netgestalt.dir <- '../netgestalt'
mut.dir <- '../mutation'
cna.dir <- data.dir


# QC status (bimodal samples) can be indiated using a separate cls file,
# or included in the experiment design file with column name qc.col;
# if neither is found, all samples are marked "unimodal"
# if both exist, the experiment design file supercedes
# NB: QC fail == bimodal; QC pass == unimodal
bimodal.cls <- NULL
qc.col <- 'QC.status'

# CNA & RNA data
cna.version <- 1
cna.data.file <- 'cna-data.gct'
rna.version <- 1
rna.data.file <- 'rnaseq-data.gct'


# missing values
na.max <- 27                  # must be present in at least 2 out of 5 TMT 10-plexes
nmiss.factor <- 0.5           # for some situations, a more stringent condition is needed
sd.filter.threshold <- 0.5    # SD threshold for SD filtering


# normalization
norm.method <- 'median'        # options: 2comp (default), median, mean
alt.method <- '2comp'          # alt.method for comparison -- filtered datasets not generated
                               # ignored if alt.method is NULL

## data subset -- whether to process the entire data set or only unimodal samples
#  (applicable for analysis that uses a specific dataset)
master.prefix <- paste (type, '-ratio-norm-NAmax', na.max, subset.str, sep='')
master.file <- file.path (norm.dir, paste (master.prefix, '.gct', sep=''))



## gene mapping not needed -- use SM geneSymbol (but map to official symbols for CNA analysis)
official.genesyms <- 'gene-symbol-map.csv'
# policy for combining/collapsing duplicate gene names -- uncomment appropriate line to use
# duplicate.gene.policy <- ifelse (type == 'phosphoproteome', 'median', 'maxvar')  
duplicate.gene.policy <- 'maxvar'


# mutation and correlation analysis -- max number of LSF jobs
LSF.mut.jid.max <- 400



# id column (with unique row ids) in the gct file
# and other differences between proteome and phosphoproteome
if (type == 'proteome') {
  # for proteome, perform marker selection at the protein level
  id.col <- 1   # protein ID
  desc.col <- 2 # gene symbol
  min.numratio <- 2
} else {
  # for phosphoproteome, perform marker selection at the phosphosite level
  id.col <- 1    # phosphosite id
  desc.col <- 2  # gene symbol
  min.numratio <- 1
}


gene.id.col <- 'geneSymbol'
protein.gene.map <- 'RefSeq.20170701-RefSeq-GeneName-map.txt'


# data source -- for managing some operations (esp related to sample IDs and names)
#  [all current options listed below -- uncomment only one]
# use project.name to manage project specific processing and options
project.name <- 'default'
# project.name <- 'cptac2.tcga'



# definitions for label type, used in parsing input SM table
# set label.type to appropriate value -- all other variables will be assigned as shown below
# use this section to define new label.types

label.type <- 'TMT10'   # alternatives: iTRAQ4

## TMT-10
if (label.type == 'TMT10') {
  n.channels <- 10   # number of columns per experiment
  # match the following channel names in the experiment design file to find sample names
  plex.channels <- c('126','127N','127C','128N','128C','129N','129C','130N','130C')
  # the above does not contain the reference channel:
  ref.channel <- '131'
  ## regex patterns for various matching
  header.pat <- '.*TMT*'
  ratio.pat <- '.*median$'
  intensity.pat <- '^TMT_1[23][67890][NC]*_total$'
  numratio.pat <- '.*numRatios.*_131$'
  numspectra.pat <-'^num_?Spectra$'
  totalint.pat <- '^totalIntensity$'
  unique_pep.pat <- '^unique_peptides$'
  refint.pat <- '^TMT_131_total$'
}

## iTRAQ-4
if (label.type == 'iTRAQ4') {
  n.channels <- 4   # number of columns per experiment
  # match the following channel names in the experiment design file to find sample names
  plex.channels <- c('114', '115', '116')
  # the above does not contain the reference channel:
  ref.channel <- '117'
  ## regex patterns for various matching
  header.pat <- '.*iTRAQ*'
  ratio.pat <- '.*median$'
  intensity.pat <- '^iTRAQ_11[456]_total$'
  numratio.pat <- '.*numRatios.*_117$'
  numspectra.pat <-'^num_?Spectra$'
  totalint.pat <- '^totalIntensity$'
  unique_pep.pat <- '^unique_peptides$'
  refint.pat <- '^iTRAQ_117_total$'
}


