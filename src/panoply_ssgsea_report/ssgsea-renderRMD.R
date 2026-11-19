#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
## Render the ssGSEA / PTM-SEA HTML report from the results tarball of panoply_ssgsea.
## Heatmap parameters are read from the 'panoply_ssgsea_report' section of a YAML file;
## any parameter passed on the command line overrides its YAML value.
## Parameters set in neither place fall back to the ssGSEA3::plot_ssgsea_heatmap() defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(ssGSEA3)
  library(ComplexHeatmap)
  library(rmarkdown)
})

as_flag <- function(x) {
  b <- as.logical(x)
  if (length(b) != 1 || is.na(b)) stop("Cannot interpret '", x, "' as TRUE/FALSE.")
  b
}

## Parameters that can be set in the YAML file and on the command line:
## names are the YAML keys / CLI flags, values the plot_ssgsea_heatmap() arguments.
params <- c(
  fdr                 = "fdr.max",               # max. FDR; also marks significant cells with '*'
  top_n               = "top_n",                 # max. number of significant gene sets plotted per direction
  cluster_rows        = "cluster_rows",          # TRUE/FALSE
  split_by_prefix     = "split_by_db_category",  # TRUE/FALSE: split rows by '<prefix>-<name>' (e.g. PTMsigDB)
  geneset_groups_file = "custom_categories"      # 2-column file: gene set name -> category
)
convert <- list(fdr = as.numeric, top_n = as.integer, cluster_rows = as_flag,
                split_by_prefix = as_flag, geneset_groups_file = as.character)

option_list <- c(
  list(
    make_option("--tarball", type = "character", help = "Path to panoply_ssgsea results tarball."),
    make_option("--label", type = "character", help = "Label used in the name of the report."),
    make_option("--yaml", type = "character", help = "Parameter file (.yaml).")
  ),
  lapply(names(params), function(p) make_option(paste0("--", p), type = "character"))
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$tarball) || is.null(opt$label)) stop("--tarball and --label are required.")

## ---- merge parameters: YAML provides defaults, command line overrides ----
yml <- list()
if (!is.null(opt$yaml)) {
  yml <- read_yaml(opt$yaml)[["panoply_ssgsea_report"]]
  if (is.null(yml)) stop("No 'panoply_ssgsea_report' section found in ", opt$yaml)
}
yml <- Filter(Negate(is.null), yml[names(yml) %in% names(params)])  # drops unset ('NULL') and unknown keys
cli <- Filter(Negate(is.null), opt[names(params)])
val <- modifyList(yml, cli)
val <- Map(function(v, p) convert[[p]](v), val, names(val))

cat("Heatmap parameters:\n")
cat(sprintf("  %-20s %-16s (%s)\n", names(val), sapply(val, paste, collapse = ","),
            ifelse(names(val) %in% names(cli), "command line", "yaml")), sep = "")

## ---- import ssGSEA3 results from the tarball ----
tmp.dir <- tempdir()
untar(opt$tarball, exdir = tmp.dir)
rds <- dir(tmp.dir, pattern = "-result\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(rds) == 0)
  stop("No '*-result.rds' found in ", opt$tarball, ". Was it created with the ssGSEA3 version of panoply_ssgsea?")
res <- readRDS(rds[1])

param.file <- dir(tmp.dir, pattern = "parameters\\.txt$", full.names = TRUE, recursive = TRUE)
param <- if (length(param.file) > 0) readLines(param.file[1]) else "no parameter file found"

## ---- heatmap ----
ht <- do.call(plot_ssgsea_heatmap, c(list(result = res), setNames(val, params[names(val)])))

fn.png <- file.path(getwd(), paste0(opt$label, "_ssGSEA_heatmap.png"))
sz <- ht_size(ht)  # size of the heatmap (without legend); add room for legend and padding
png(fn.png, units = "in", res = 300,
    width = grid::convertWidth(sz$width, "in", valueOnly = TRUE) + 2,
    height = grid::convertHeight(sz$height, "in", valueOnly = TRUE) + 1)
draw(ht)
invisible(dev.off())

## ---- render report ----
## values shown in the report: those used for the heatmap, else the plot_ssgsea_heatmap() defaults
defaults <- formals(plot_ssgsea_heatmap)
script.dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
render(file.path(script.dir, "ssgsea_rmd.rmd"),
       params = list(title = paste0("ssGSEA Report - ", opt$label),
                     label = opt$label,
                     fdr = if (is.null(val$fdr)) defaults$fdr.max else val$fdr,
                     top_n = if (is.null(val$top_n)) defaults$top_n else val$top_n,
                     param = param, fn.png = fn.png),
       output_dir = getwd(), output_file = paste0(opt$label, "_ssGSEA_rmd.html"),
       intermediates_dir = tmp.dir, envir = new.env(), quiet = TRUE)
