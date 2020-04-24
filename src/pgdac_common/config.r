

### Broad Institute PGDAC PIPELINE
### configuration and default parameters  

## Many parameters here are used and/or assumed to have specific values in
## run-pipeline.sh. If editing parameters before the "Pipeline Parameters" section,
## cross check run-pipeline.sh to make corresponding changes


## Path for R-utilites (for I/O and other misc functions)
Rutil.path <- switch (Sys.info()[['sysname']],
                      Windows = {'//flynn-cifs/prot_proteomics/Projects/R-utilities'},
                      Darwin = {'/Volumes/prot_proteomics/Projects/R-utilities'},
                      Linux = {'/prot/proteomics/Projects/R-utilities'})

Source <- function (f) {
  # adaptive source file location -- from current directory or from R-utiliites
  if (file.exists (f)) {source (f)}
  else {
    f.os <- file.path (Rutil.path, f)
    if (file.exists (f.os)) source (f.os)
    else stop (paste ("Can't find file", f.os, '-- R-utilities missing?'))
  }
}


Source ('gct-io.r')
Source ('io.r')

# Load 'yaml' package with pacman to read in master parameter file:
if( !suppressMessages( require( "pacman" ) ) ) install.packages( "pacman" )
p_load('yaml')
#yml_dir <- getwd() # find the parameters file directory
#yml_dir <- file.path("master_parameters.yml") # path to file
master.parameters <- read_yaml("master-parameters.yml")

## Directory structure
# directories with raw, pre-processed and normalized data
# also includes other analysis directories
# directory structure is created by appropriate options in run-pipeline.sh
# (changes to this will require corresponding changes in run-pipeline.sh)
data.dir <- master.parameters$config_r_parameters$data_dir #'../data'
pre.dir <- master.parameters$config_r_parameters$pre_dir #'../parsed-data'
norm.dir <- master.parameters$config_r_parameters$norm_dir #'../normalized-data'
harmonize.dir <- master.parameters$config_r_parameters$harmonize_dir #'../harmonized-data'


## Command line arguments
# determine type (proteome/phosphoproteome) based on command line args; defaults to proteome
# data subset defaults to QC.pass, unless specified; if QC.pass, subset.str is NULL
args <- commandArgs (trailingOnly=TRUE)
type <- ifelse (length(args) > 0, toString (args[1]), "proteome")
data.subset <- ifelse (length(args) > 1, toString (args[2]), "")
subset.str <- ifelse (data.subset=="", data.subset, paste ('-', data.subset, sep='')) 

## data subset -- whether to process the entire data set or only QC.pass samples
#  (applicable for analysis that uses a specific dataset)
master.prefix <- paste (type, '-ratio-norm-NArm', subset.str, sep='')
master.file <- file.path (norm.dir, paste (master.prefix, '.gct', sep=''))

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



## Label type definitions
# definitions for label type, used in parsing input SM table
# set label.type to appropriate value -- all other variables will be assigned as shown below
# use this section to define new label.types


## Label type for MS experiment
set.label.type <- function (label.type) {
  ## function to set variables based on label.type
  ## assign is used to set the value globally
  ## options: TMT10 and iTRAQ4
  
  ## TMT-10 (default reference channel in 131)
  if (label.type == 'TMT10') {
    assign ("n.channels",  10, envir = .GlobalEnv)   # number of columns per experiment
    # match the following channel names in the experiment design file to find sample names
    assign ("plex.channels",  c('126','127N','127C','128N','128C','129N','129C','130N','130C'), envir = .GlobalEnv)
    # the above does not contain the reference channel:
    assign ("ref.channel",  '131', envir = .GlobalEnv)
    ## regex patterns for various matching
    assign ("header.pat",  '.*TMT*', envir = .GlobalEnv)
    assign ("ratio.pat",  '.*median$', envir = .GlobalEnv)
    assign ("stddev.pat",  '.*stddev$', envir = .GlobalEnv)
    assign ("intensity.pat",  '^TMT_1[23][67890][NC]*_total$', envir = .GlobalEnv)
    assign ("numratio.pat",  '.*numRatios.*_131$', envir = .GlobalEnv)
    assign ("numspectra.pat", '^num_?Spectra$', envir = .GlobalEnv)
    assign ("totalint.pat",  '^totalIntensity$', envir = .GlobalEnv)
    assign ("unique_pep.pat",  '^unique_peptides$', envir = .GlobalEnv)
    assign ("refint.pat",  '^TMT_131_total$', envir = .GlobalEnv)
  }

  ## TMT-11 (default reference channel in 131C)
  if (label.type == 'TMT11') {
    assign ("n.channels",  11, envir = .GlobalEnv)   # number of columns per experiment
    # match the following channel names in the experiment design file to find sample names
    assign ("plex.channels",  c('126C','127N','127C','128N','128C','129N','129C','130N','130C','131N'), envir = .GlobalEnv)
    # the above does not contain the reference channel:
    assign ("ref.channel",  '131C', envir = .GlobalEnv)
    ## regex patterns for various matching
    assign ("header.pat",  '.*TMT*', envir = .GlobalEnv)
    assign ("ratio.pat",  '.*median$', envir = .GlobalEnv)
    assign ("stddev.pat",  '.*stddev$', envir = .GlobalEnv)
    assign ("intensity.pat",  '^TMT_1[23][678901][NC]*_total$', envir = .GlobalEnv)
    assign ("numratio.pat",  '.*numRatios.*_131C$', envir = .GlobalEnv)
    assign ("numspectra.pat", '^num_?Spectra$', envir = .GlobalEnv)
    assign ("totalint.pat",  '^totalIntensity$', envir = .GlobalEnv)
    assign ("unique_pep.pat",  '^unique_peptides$', envir = .GlobalEnv)
    assign ("refint.pat",  '^TMT_131C_total$', envir = .GlobalEnv)
  }

  ## TMT-10 with 126 as reference channel
  if (label.type == 'TMT10.126') {
    assign ("n.channels",  10, envir = .GlobalEnv)   # number of columns per experiment
    # match the following channel names in the experiment design file to find sample names
    assign ("plex.channels",  c('127N','127C','128N','128C','129N','129C','130N','130C', '131'), envir = .GlobalEnv)
    # the above does not contain the reference channel:
    assign ("ref.channel",  '126', envir = .GlobalEnv)
    ## regex patterns for various matching
    assign ("header.pat",  '.*TMT*', envir = .GlobalEnv)
    assign ("ratio.pat",  '.*median$', envir = .GlobalEnv)
    assign ("stddev.pat",  '.*stddev$', envir = .GlobalEnv)
    assign ("intensity.pat",  '^TMT_1[23][17890][NC]*_total$', envir = .GlobalEnv)
    assign ("numratio.pat",  '.*numRatios.*_126$', envir = .GlobalEnv)
    assign ("numspectra.pat", '^num_?Spectra$', envir = .GlobalEnv)
    assign ("totalint.pat",  '^totalIntensity$', envir = .GlobalEnv)
    assign ("unique_pep.pat",  '^unique_peptides$', envir = .GlobalEnv)
    assign ("refint.pat",  '^TMT_126_total$', envir = .GlobalEnv)
  }
  
  ## iTRAQ-4
  if (label.type == 'iTRAQ4') {
    assign ("n.channels",  4, envir = .GlobalEnv)   # number of columns per experiment
    # match the following channel names in the experiment design file to find sample names
    assign ("plex.channels",  c('114', '115', '116'), envir = .GlobalEnv)
    # the above does not contain the reference channel:
    assign ("ref.channel",  '117', envir = .GlobalEnv)
    ## regex patterns for various matching
    assign ("header.pat",  '.*iTRAQ*', envir = .GlobalEnv)
    assign ("ratio.pat",  '.*median$', envir = .GlobalEnv)
    assign ("stddev.pat",  '.*stddev$', envir = .GlobalEnv)
    assign ("intensity.pat",  '^iTRAQ_11[456]_total$', envir = .GlobalEnv)
    assign ("numratio.pat",  '.*numRatios.*_117$', envir = .GlobalEnv)
    assign ("numspectra.pat", '^num_?Spectra$', envir = .GlobalEnv)
    assign ("totalint.pat",  '^totalIntensity$', envir = .GlobalEnv)
    assign ("unique_pep.pat",  '^unique_peptides$', envir = .GlobalEnv)
    assign ("refint.pat",  '^iTRAQ_117_total$', envir = .GlobalEnv)
  }
}


## Input data location
# input files are copies to this location so that the tarball has all the data
input.data.file <- file.path (data.dir, paste (type, master.parameters$config_r_parameters$input_data_file, sep=''))
expt.design.file <- file.path (data.dir, master.parameters$config_r_parameters$expt_design_file)
rna.data.file <- file.path (data.dir, master.parameters$config_r_parameters$rna_data_file)
cna.data.file <- file.path (data.dir, master.parameters$config_r_parameters$cna_data_file)



### Pipeline parameters
# defaults shown below
# can be overridden by redefining the appropriate variable in the
# parameter redefinition section
###


## Label type for MS experiment (set.label.type must be called to initialize variables)
label.type <- master.parameters$config_r_parameters$label_type_for_MS_exp$label_type   # default TMT10, alternatives: iTRAQ4, TMT10.126, TMT11
set.label.type (label.type) 

## Sample replicate indicator
# Sample.IDs MUST be unique in the expt.design.file; duplicate samples should have the same 
# sample names, but include this replicate.indicator, followed by a unique suffix: <name>REP1)
replicate.indicator <- master.parameters$config_r_parameters$sample_replicate_indicator$replicate_indicator

## QC
# QC status can be indiated using a separate cls file,
# or included in the experiment design file with column name qc.col;
# if neither is found, all samples are marked as qc.pass.label
# if both exist, the experiment design file supercedes
sampleQC.cls <- master.parameters$config_r_parameters$QC$sampleQC_cls
qc.col <- master.parameters$config_r_parameters$QC$qc_col
qc.pass.label <- master.parameters$config_r_parameters$QC$qc_pass_label


## Output precision for gct tables
ndigits <- master.parameters$config_r_parameters$output_precision$ndigits   


## Missing values and filtering
na.max <- master.parameters$config_r_parameters$missing_values_and_filtering$na_max                 # maximum allowed NA values (per protein/site/row), can be fraction or integer number of samples
min.numratio.fraction <- master.parameters$config_r_parameters$missing_values_and_filtering$min_numratio_fraction  # fraction of samples in which min. numratio should be present to retain protein/phosphosite
sample.na.max <- master.parameters$config_r_parameters$missing_values_and_filtering$sample_na_max           # maximum allowed fraction of NA values per sample/column; pipeline error if violated
nmiss.factor <- master.parameters$config_r_parameters$missing_values_and_filtering$nmiss_factor            # for some situations, a more stringent condition is needed
sd.filter.threshold <-master.parameters$config_r_parameters$missing_values_and_filtering$sd_filter_threshold     # SD threshold for SD filtering
clustering.sd.threshold <- master.parameters$config_r_parameters$missing_values_and_filtering$clustering_sd_threshold   # threshold for filtering data before consensus clustering
clustering.na.threshold <- master.parameters$config_r_parameters$missing_values_and_filtering$clustering_na_threshold # max fraction of missing values for clustering; rest are imputed
apply.SM.filter <- master.parameters$config_r_parameters$missing_values_and_filtering$apply_SM_filter        # if TRUE, apply numRatio based filter (use TRUE if input is SM ssv)


## Normalization
norm.method <- master.parameters$config_r_parameters$normalization$norm_method         # options: 2comp (default), median, mean
alt.method <- master.parameters$config_r_parameters$normalization$alt_method         # alt.method for comparison -- filtered datasets not generated
if (norm.method == alt.method) alt.method <- NULL
                               # ignored if alt.method is NULL, or is identical to norm.method


## Gene mapping
# gene mapping not needed -- use SM geneSymbol (but map to official symbols for CNA analysis)
official.genesyms <- master.parameters$config_r_parameters$gene_mapping$official_genesyms
gene.id.col <- master.parameters$config_r_parameters$gene_mapping$gene_id_col
protein.gene.map <- master.parameters$config_r_parameters$gene_mapping$protein_gene_map
# policy for combining/collapsing duplicate gene names -- uncomment appropriate line to use
# duplicate.gene.policy <- ifelse (type == 'phosphoproteome', 'median', 'maxvar')  
duplicate.gene.policy <- master.parameters$config_r_parameters$gene_mapping$duplicate_gene_policy

## RNA related
rna.output.prefix <- master.parameters$config_r_parameters$rna$rna_output_prefix  # output prefix for tables creates during RNA analysis
rna.sd.threshold <- master.parameters$config_r_parameters$rna$rna_sd_threshold   # for variation filter (set to NA to disable)

## CNA/parallelism related
pe.max.default <- master.parameters$config_r_parameters$CNA_parallelism$pe_max_default     # default maximum processors/jobs

## Project
# data source -- for managing some operations (esp related to sample IDs and names)
#  [all current options listed below -- uncomment only one]
# use project.name to manage project specific processing and options
project.name <- master.parameters$config_r_parameters$project$project_name
# project.name <- 'cptac2.tcga'

## Disease
# disease setting is used to set disease specific options and
# run appropriate CNA subsets by creating run-cna-analysis-<disease>.r;
# by default disease is empty ('') this will check and only define it if the parameter is changed in the yaml file
if (!is.null(master.parameters$config_r_parameters$disease$disease)){
  disease <- master.parameters$config_r_parameters$disease$disease
}
