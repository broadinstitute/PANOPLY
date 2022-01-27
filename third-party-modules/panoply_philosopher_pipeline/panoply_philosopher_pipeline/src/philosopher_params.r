#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
# parameter_manager.r
### Commandline inputs:
library(optparse)
library(yaml)

option_list <- list(
    # Every use:
    make_option(c("--yaml"), type = "character", dest = 'yaml', help = "philosopher.yml file to provide defaults"),
    # make_option(c("--module"), type = "character", dest = 'module', help = "The module or workflow name."),

    # GLOBAL
    make_option(c("--memory"), type = "integer", dest = 'ndigits', help = "how much memory in GB to use"),
    make_option(c("--num_threads"), type = "integer", dest = 'num_threads', help = "number of CPU threads to use. 0=poll CPU to set num threads"),

    # Database Search (MSFragger):
    make_option(c("--protein_database"), type = "character", dest = 'protein_database', help = "path to the target-decoy protein database"),
    make_option(c("--decoy_tag"), type = "character", dest = 'decoy_tag', help = "prefix tag used added to decoy sequences"),
    make_option(c("--contaminant_tag"), type = "character", dest = 'contaminant_tag', help = "prefix tag used added to contaminant sequences"),

    make_option(c("--precursor_mass_lower"), type = "integer", dest = 'precursor_mass_lower', help = "lower bound of the precursor mass window"),
    make_option(c("--precursor_mass_upper"), type = "integer", dest = 'precursor_mass_upper', help = "upper bound of the precursor mass window"),

    make_option(c("--variable_mod_01"), type = "character", dest = 'variable_mod_01', help = "variable modification"),
    make_option(c("--variable_mod_02"), type = "character", dest = 'variable_mod_02', help = "variable modification"),
    make_option(c("--variable_mod_03"), type = "character", dest = 'variable_mod_03', help = "variable modification"),
    make_option(c("--variable_mod_04"), type = "character", dest = 'variable_mod_04', help = "variable modification"),
    make_option(c("--variable_mod_05"), type = "character", dest = 'variable_mod_05', help = "variable modification"),
    make_option(c("--variable_mod_06"), type = "character", dest = 'variable_mod_06', help = "variable modification"),
    make_option(c("--variable_mod_07"), type = "character", dest = 'variable_mod_07', help = "variable modification"),
    make_option(c("--max_variable_mods_per_peptide"), type = "integer", dest = 'max_variable_mods_per_peptide', help = "maximum number of variable modifications on a peptide"),

    make_option(c("--digest_min_length"), type = "integer", dest = 'digest_min_length', help = "minimum length of peptides to be generated during in-silico digestion"),
    make_option(c("--digest_max_length"), type = "integer", dest = 'digest_max_length', help = "maximum length of peptides to be generated during in-silico digestion"),
    make_option(c("--digest_mass_range"), type = "character", dest = 'digest_mass_range', help = "mass range of peptides to be generated during in-silico digestion in Daltons"),

    # Peptide Validation (PeptideProphet)

    # PTM Localization (PTMProphet)

    # Protein Inference (ProteinProphet)

    # FDR Filtering 
    make_option(c("--psmFDR"), type = "double", dest = 'psmFDR', help = "psm FDR level (default 0.01)"),
    make_option(c("--peptideFDR"), type = "double", dest = 'peptideFDR', help = "peptide FDR level (default 0.01)"),
    make_option(c("--ionFDR"), type = "double", dest = 'ionFDR', help = "peptide ion FDR level (default 0.01)"),
    make_option(c("--proteinFDR"), type = "double", dest = 'proteinFDR', help = "protein FDR level (default 0.01)")

    # Integrated Isobaric Quantification (TMT-Integrator)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# MODULE NAMES:
# global
# msfragger
# peptideprophet
# ptmprophet
# proteinprophet
# filter

### FUNCTIONS:

### HELPER FUNCTIONS:
# all values are optional other than help, module, and yaml so once those are removed
# opt will be of length 0 if no command line overwrites were given!
check_if_any_command_line <- function(opt) {
  opt_copy <- opt
  opt_copy$help <- NULL
  opt_copy$yaml <- NULL
  if (length(opt_copy) == 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Read in the yaml file:
read_yaml_file <- function(yaml) {
  open_yaml <- read_yaml(yaml)
  return(open_yaml)
}

# Write out final yaml file:
write_yaml_file <- function(yaml) {
  yaml_name <- file("philosopher.yml", "w")
  write_yaml(yaml, yaml_name)
  close(yaml_name)
}

######################################## CHECK PARAMS FUNCTIONS ################################
# checks if any of the global vals are being changed from commandline. if yes will change them in yaml if no returns org yaml
check_global_params <- function(opt, yaml) {
    if (!is.null(opt$memory)){
        yaml$`Database Search`$msfragger$memory <- opt$memory
        yaml$`Integrated Isobaric Quantification`$memory <- opt$memory
    }
    if (!is.null(opt$num_threads)){
        yaml$`Database Search`$msfragger$num_threads <- opt$num_threads
    }
    return(yaml)
}

# msfragger:
check_msfragger_params <- function(opt, yaml) {
    if (!is.null(opt$protein_database)){
        yaml$`Database Search`$protein_database <- opt$protein_database
    }
    if (!is.null(opt$decoy_tag)){
        yaml$`Database Search`$decoy_tag <- opt$decoy_tag
    }
    if (!is.null(opt$contaminant_tag)){
        yaml$`Database Search`$contaminant_tag <- opt$contaminant_tag
    }

    if (!is.null(opt$precursor_mass_lower)){
        yaml$`Database Search`$msfragger$precursor_mass_lower <- opt$precursor_mass_lower
    }
    if (!is.null(opt$variable_mod_01)) {
        yaml$`Database Search`$msfragger$variable_mod_01 <- opt$variable_mod_01
    }
    if (!is.null(opt$variable_mod_02)) {
        yaml$`Database Search`$msfragger$variable_mod_02 <- opt$variable_mod_02
    }
    if (!is.null(opt$variable_mod_03)) {
        yaml$`Database Search`$msfragger$variable_mod_03 <- opt$variable_mod_03
    }
    if (!is.null(opt$variable_mod_04)) {
        yaml$`Database Search`$msfragger$variable_mod_04 <- opt$variable_mod_04
    }
    if (!is.null(opt$variable_mod_05)) {
        yaml$`Database Search`$msfragger$variable_mod_05 <- opt$variable_mod_05
    }
    if (!is.null(opt$variable_mod_06)) {
        yaml$`Database Search`$msfragger$variable_mod_06 <- opt$variable_mod_06
    }
    if (!is.null(opt$variable_mod_07)) {
        yaml$`Database Search`$msfragger$variable_mod_07 <- opt$variable_mod_07
    }
    if (!is.null(opt$max_variable_mods_per_peptide)) {
        yaml$`Database Search`$msfragger$max_variable_mods_per_peptide <- opt$max_variable_mods_per_peptide
    }
    if (!is.null(opt$precursor_mass_upper)) {
        yaml$`Database Search`$msfragger$precursor_mass_upper <- opt$precursor_mass_upper
    }
    if (!is.null(opt$digest_min_length)) {
        yaml$`Database Search`$msfragger$digest_min_length <- opt$digest_min_length
    }
    if (!is.null(opt$digest_max_length)) {
        yaml$`Database Search`$msfragger$digest_max_length <- opt$digest_max_length
    }
    if (!is.null(opt$digest_mass_range)) {
        yaml$`Database Search`$msfragger$digest_mass_range <- opt$digest_mass_range
    }
    return(yaml)
}

# peptideprophet:
check_peptideprophet_params <- function(opt, yaml) {
    return(yaml)
}

# ptmprophet:
check_ptmprophet_params <- function(opt, yaml) {
    return(yaml)
}

# proteinprophet:
check_proteinprophet_params <- function(opt, yaml) {
    return(yaml)
}

# filter:
check_filter_params <- function(opt, yaml) {
    if (!is.null(opt$psmFDR)){
        yaml$`FDR Filtering`$msfragger$psmFDR <- opt$psmFDR
    }
    if (!is.null(opt$peptideFDR)){
        yaml$`FDR Filtering`$msfragger$peptideFDR <- opt$peptideFDR
    }
    if (!is.null(opt$ionFDR)){
        yaml$`FDR Filtering`$msfragger$ionFDR <- opt$ionFDR
    }
    if (!is.null(opt$proteinFDR)){
        yaml$`FDR Filtering`$msfragger$proteinFDR <- opt$proteinFDR
    }
    return(yaml)
}

# tmtintegrator:
check_tmtintegrator_params <- function(opt, yaml){
    return(yaml)
}

########################################################################################

# Main function and logic
# How to handle the commandline variables based on current module?
parse_command_line_parameters <- function(opt) {
    yaml <- read_yaml_file(opt$yaml)
    if (check_if_any_command_line(opt)) { # TRUE if command line changes are being made
        yaml <- check_global_params(opt, yaml)
        yaml <- check_msfragger_params(opt, yaml) 
        yaml <- check_peptideprophet_params(opt, yaml)   
        yaml <- check_ptmprophet_params(opt, yaml)
        yaml <- check_proteinprophet_params(opt, yaml)
        yaml <- check_filter_params(opt, yaml)
        yaml <- check_tmtintegrator_params(opt, yaml)
    }
    write_yaml_file(yaml)
}


# MAIN:
parse_command_line_parameters(opt)
