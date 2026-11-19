#!/usr/bin/env Rscript
#
# Copyright (c) 2026 The Broad Institute, Inc. All rights reserved.
#
## Command line interface to ssGSEA3 / PTM-SEA.
## Parameter defaults are read from the 'panoply_ssgsea' section of a YAML file;
## any parameter passed on the command line overrides its YAML value.
## Parameters set in neither place fall back to the ssGSEA3() defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(ssGSEA3)
})

## Parameters that can be set in the YAML file and on the command line (name = type).
## Names are the YAML keys and the CLI flags; they map to ssGSEA3() arguments by replacing
## '_' with '.' (except multi_core, which maps to 'par').
params <- c(
  output_prefix     = "character",
  sample_norm_type  = "character",  # rank, log, log.rank, none
  correl_type       = "character",  # rank, z.score, symm.rank
  statistic         = "character",  # area.under.RES, Kolmogorov-Smirnov
  output_score_type = "character",  # NES, ES
  weight            = "numeric",
  nperm             = "integer",
  min_overlap       = "integer",
  global_fdr        = "logical",
  extended_output   = "logical",
  export_signat_gct = "logical",
  multi_core        = "logical"
)

option_list <- c(
  list(
    make_option("--input", type = "character", help = "Path to input GCT file."),
    make_option("--db", type = "character", help = "Path to gene set database (GMT format)."),
    make_option("--yaml", type = "character", help = "Parameter file (.yaml)."),
    make_option("--tolerate_min_overlap_err", type = "character", default = "false",
                help = "TRUE/FALSE: if the data has < min_overlap features in common with every gene set, write a warning instead of failing.")
  ),
  lapply(names(params), function(p) make_option(paste0("--", p), type = "character"))
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input) || is.null(opt$db)) stop("--input and --db are required.")

as_flag <- function(x) {
  b <- as.logical(x)
  if (length(b) != 1 || is.na(b)) stop("Cannot interpret '", x, "' as TRUE/FALSE.")
  b
}
convert <- list(character = as.character, numeric = as.numeric, integer = as.integer, logical = as_flag)

## ---- merge parameters: YAML provides defaults, command line overrides ----
yml <- list()
if (!is.null(opt$yaml)) {
  yml <- read_yaml(opt$yaml)[["panoply_ssgsea"]]
  if (is.null(yml)) stop("No 'panoply_ssgsea' section found in ", opt$yaml)
}
yml <- yml[names(yml) %in% names(params)]
cli <- Filter(Negate(is.null), opt[names(params)])
val <- modifyList(yml, cli)
if (is.null(val$output_prefix)) stop("'output_prefix' must be set in the YAML file or with --output_prefix.")
val <- Map(function(v, type) convert[[type]](v), val, params[names(val)])

cat("ssGSEA3 parameters:\n")
cat(sprintf("  %-18s %-22s (%s)\n", names(val), sapply(val, paste, collapse = ","),
            ifelse(names(val) %in% names(cli), "command line", "yaml")), sep = "")

## ---- run ssGSEA3 ----
args <- c(
  list(input.ds = opt$input, gene.set.databases = opt$db,
       spare.cores = 0,  # use all available cpus when multi_core = TRUE
       log.file = paste0(val$output_prefix, "_ssgsea.log.txt")),
  setNames(val, sub("^multi\\.core$", "par", gsub("_", ".", names(val))))
)
tolerate <- as_flag(opt$tolerate_min_overlap_err)

res <- tryCatch(do.call(ssGSEA3, args), error = function(e) {
  if (tolerate && grepl("does not meet minimum-overlap", conditionMessage(e), fixed = TRUE)) {
    message("\n### WARNING\n", conditionMessage(e))
    return(NULL)
  }
  stop(e)
})

## save the full result object, so downstream modules (e.g. panoply_ssgsea_report) can plot from it
if (!is.null(res)) saveRDS(res, paste0(val$output_prefix, "-result.rds"))

## declare whether a tolerated min-overlap error occurred (read by the WDL)
writeLines(as.character(is.null(res)), "geneset_overlap_below_min.txt")
