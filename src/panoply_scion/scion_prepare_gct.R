#!/usr/bin/env Rscript
# Prepares a target or regulator matrix for SCION from a GCT file.
#
# SCION::read_scion_inputs(format = "gct") uses a GCT's @mat row names as-is,
# which is fine for a proteome-level regulator matrix with clean gene-symbol
# row IDs, but not for a PTM ome (phosphoproteome/acetylome/ubiquitylome),
# where the raw row ID is a SpectrumMill- or FragPipe-specific composite
# string that doesn't cleanly split into "gene symbol" + "site" the way
# SCION's own ptm_sep convention expects. This re-derives row names from the
# GCT's own rdesc metadata instead, adding the PTM site back on when needed.
#
# Ported from the pre-SCION-package panoply_scion task's
# scion_data_processing.R, generalized so the gene-symbol rdesc column is a
# parameter (gene_id_col) instead of a hardcoded "geneSymbol", matching the
# rest of PANOPLY's gene_id_col convention (see panoply_harmonize/harmonize.r).

#' @param gct_file path to the GCT(x) file.
#' @param role "target" (mRNA -- used as-is, no gene_id_col/PTM handling) or
#'   "reg" (proteomics -- row names always re-derived from gene_id_col).
#' @param gene_id_col name of the rdesc column holding the gene symbol.
#'   Ignored for role = "target".
#' @param is_ptm whether this is a PTM ome (phosphoproteome/acetylome/
#'   ubiquitylome). Only meaningful for role = "reg".
#' @param ptm_type "SM" (SpectrumMill) or "FP" (FragPipe): which convention
#'   the GCT's raw row IDs follow for encoding a PTM site. Only used when
#'   is_ptm = TRUE.
#' @param ptm_sep separator used to join gene symbol + site into a row name
#'   (must match the ptm_sep passed to SCION::run_scion() downstream, so
#'   SCION's own split_ptm_sites() can split it back apart later).
#' @return a data frame (genes/sites as rows, samples as columns).
prepare_scion_matrix <- function(gct_file, role = c("target", "reg"), gene_id_col = "geneSymbol",
                                  is_ptm = FALSE, ptm_type = c("SM", "FP"), ptm_sep = ".") {
  role <- match.arg(role)
  ptm_type <- match.arg(ptm_type)
  gct <- cmapR::parse_gctx(gct_file)
  # kept as a plain matrix, not a data frame, for as long as possible: a
  # matrix's row names don't need to be unique until the moment they're
  # actually assigned, so duplicate-gene/duplicate-site rows can be filtered
  # out first, in one pass, instead of erroring the instant a duplicated name
  # is assigned to a data frame.
  mat <- gct@mat

  if (role == "target") {
    return(as.data.frame(mat))
  }

  if (!gene_id_col %in% colnames(gct@rdesc)) {
    stop("gene_id_col '", gene_id_col, "' not found in ", gct_file, "'s row descriptors. ",
         "Available columns: ", paste(colnames(gct@rdesc), collapse = ", "))
  }
  gene_symbols <- as.character(gct@rdesc[[gene_id_col]])

  if (!is_ptm) {
    # non-PTM proteome: row names come from gene_id_col (not necessarily the
    # same as the raw GCT row ID); rows mapping to the same gene symbol are
    # averaged together
    keep <- !is.na(gene_symbols) & gene_symbols != ""
    mat <- mat[keep, , drop = FALSE]
    gene_symbols <- gene_symbols[keep]
    if (!anyDuplicated(gene_symbols)) {
      rownames(mat) <- gene_symbols
      return(as.data.frame(mat))
    }
    # aggregate() groups by gene_symbols directly, so mat's own row names
    # never need to be unique at any point
    averaged <- stats::aggregate(as.data.frame(mat), by = list(gene = gene_symbols),
                                  FUN = function(x) mean(x, na.rm = TRUE))
    row.names(averaged) <- averaged$gene
    averaged$gene <- NULL
    return(averaged)
  }

  # PTM ome: extract the site from the RAW row ID. This is always
  # underscore-delimited regardless of ptm_sep -- that's the ID format the
  # search engine itself writes (a fixed convention per tool), not something
  # PANOPLY exposes as configurable -- and rebuild "<gene><ptm_sep><site>",
  # which IS in SCION's own (configurable) convention.
  raw_ids <- rownames(mat)
  sites <- vapply(strsplit(raw_ids, "_", fixed = TRUE), function(tokens) {
    # trim whitespace before testing -- some (older?) GCTs have a trailing
    # space on the site token (e.g. "S20s " instead of "S20s"), which made the
    # "ends with a lowercase letter" check below fail to match at all if
    # tested on the raw, untrimmed token
    tokens <- trimws(tokens)
    if (ptm_type == "SM") {
      # SpectrumMill site tokens start with an uppercase letter and end with
      # a lowercase letter (e.g. "S123s") -- an empirical convention (not a
      # documented spec), ported as-is from the pre-package task.
      candidates <- tokens[grepl("^[[:upper:]]", tokens) & grepl("[[:lower:]]$", tokens)]
      if (length(candidates) == 0) NA_character_ else candidates[[1]]
    } else {
      tokens[[length(tokens)]]
    }
  }, character(1))

  keep <- !is.na(sites) & !is.na(gene_symbols) & gene_symbols != ""
  if (!any(keep)) {
    stop("No rows in ", gct_file, " had both a mappable gene symbol (column '", gene_id_col,
         "') and an extractable PTM site (ptm_type = '", ptm_type, "').")
  }
  mat <- mat[keep, , drop = FALSE]
  new_names <- paste(gene_symbols[keep], sites[keep], sep = ptm_sep)
  dup <- duplicated(new_names)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- new_names[!dup]
  as.data.frame(mat)
}
