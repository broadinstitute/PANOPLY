#!/usr/bin/env Rscript
#
# Copyright (c) 2026 The Broad Institute, Inc. All rights reserved.
#
## Command line interface to ssGSEA3::preprocessGCT(): converts a GCT file to the gene-centric
## or single-site-centric format used by ssGSEA / PTM-SEA.
## Parameter defaults are read from the 'panoply_preprocess_gct' section of a YAML file;
## any parameter passed on the command line overrides its YAML value.
## Parameters set in neither place fall back to the preprocessGCT() defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(ssGSEA3)
})

## Parameters that can be set in the YAML file and on the command line:
## names are the YAML keys / CLI flags, values the corresponding preprocessGCT() arguments.
params <- c(
  level               = "level",                # ssc, gc, gcr
  id_type             = "id.type",              # sm, ph
  id_type_out         = "id.type.out",          # uniprot, refseq, seqwin
  acc_type_in         = "acc.type",             # uniprot, refseq, symbol
  organism            = "organism",             # human, mouse, rat, zebrafish (needed to map to UniProt)
  seqwin_column       = "seqwin.col",
  gene_symbol_column  = "gene.col",
  humanize_gene_names = "humanize.gene.names",  # TRUE/FALSE
  sgt_column          = "SGT.col",
  localized           = "loc",                  # TRUE/FALSE
  mode                = "mode",                 # mean, median, sd, SGT, abs.max
  residue             = "mod.res",              # S|T|Y, K
  ptm                 = "mod.type"              # p, ac, ub
)
logicals <- c("humanize_gene_names", "localized")

option_list <- c(
  list(
    make_option("--input", type = "character", help = "Path to input GCT file."),
    make_option("--yaml", type = "character", help = "Parameter file (.yaml).")
  ),
  lapply(names(params), function(p) make_option(paste0("--", p), type = "character"))
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("--input is required.")

as_flag <- function(x) {
  b <- as.logical(x)
  if (length(b) != 1 || is.na(b)) stop("Cannot interpret '", x, "' as TRUE/FALSE.")
  b
}

## ---- merge parameters: YAML provides defaults, command line overrides ----
yml <- list()
if (!is.null(opt$yaml)) {
  yml.all <- read_yaml(opt$yaml)
  yml <- yml.all[["panoply_preprocess_gct"]]
  if (is.null(yml)) stop("No 'panoply_preprocess_gct' section found in ", opt$yaml)
  ## 'organism' is a project-wide setting: use global_parameters$organism unless set in this section
  if (is.null(yml$organism)) yml$organism <- yml.all$global_parameters$organism
}
yml <- yml[names(yml) %in% names(params)]
cli <- Filter(Negate(is.null), opt[names(params)])
val <- modifyList(yml, cli)
val <- Map(function(v, p) if (p %in% logicals) as_flag(v) else as.character(v), val, names(val))

cat("preprocessGCT parameters:\n")
cat(sprintf("  %-20s %-16s (%s)\n", names(val), sapply(val, paste, collapse = ","),
            ifelse(names(val) %in% names(cli), "command line", "yaml")), sep = "")

## ---- run preprocessGCT ----
res <- do.call(preprocessGCT, c(list(gct.str = opt$input, appenddim = FALSE),
                                setNames(val, params[names(val)])))

## name of the file written to the working directory (read by the WDL)
writeLines(attr(res, "path"), "fn.out")
