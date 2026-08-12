# Input upload, category mapping, and validation.

wb_list_data_categories <- function() {
  cat("Data categories:\n")
  for (i in seq_along(CAT_MAP)) cat(sprintf("  %2d: %s\n", i, CAT_MAP[i]))
  cat("   0: (none of the above)\n")
  flush.console()
}

wb_default_asset <- function(pattern) {
  # GSEA/PTM-SEA gmt databases live in the repo's shared defaults/ (also used by panda/), not
  # committed under workbench-src/defaults/ -- that directory is a build-time staging target
  # populated by deploy-workbench.sh right before upload, so the deployed workbench/ tree is
  # still self-contained (only workbench/ gets synced to S3, not the whole repo). The second
  # candidate is a dev-mode fallback for running straight out of a full repo checkout without
  # having run the deploy script first. The pattern match (rather than a fixed filename) just
  # means a future version bump can drop in a new file and be picked up as "most recent"
  # without a code change.
  candidate_dirs <- c(file.path("workbench-src", "defaults"), file.path("..", "defaults"))
  source_path <- NULL
  for (dir in candidate_dirs) {
    # First directory with ANY match wins entirely (staged copy takes priority over the
    # dev-mode fallback) -- within a directory, prefer the alphabetically-last match, in case
    # more than one versioned file matches the pattern.
    hits <- list.files(dir, pattern = pattern, full.names = TRUE)
    if (length(hits) > 0) { source_path <- sort(hits, decreasing = TRUE)[1]; break }
  }
  if (is.null(source_path)) {
    stop(sprintf("No default asset found matching '%s' under: %s.", pattern,
                 paste(candidate_dirs, collapse = ", ")))
  }

  seeded_dir <- file.path(wb_workbench_root(), "defaults")
  dir.create(seeded_dir, showWarnings = FALSE, recursive = TRUE)
  seeded_path <- file.path(seeded_dir, basename(source_path))
  if (!file.exists(seeded_path)) file.copy(source_path, seeded_path)
  seeded_path
}

wb_load_and_map_inputs <- function(state, input_dir = file.path(wb_workbench_root(), "inputs"),
                                    zip_path = NULL) {
  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(zip_path)) utils::unzip(zip_path, exdir = input_dir, junkpaths = TRUE)

  files <- list.files(input_dir, pattern = "\\.(gct|csv|ya?ml|gmt)$", full.names = FALSE)
  if (length(files) == 0 && is.null(zip_path)) {
    stop(sprintf(
      "No .gct/.csv/.yaml/.gmt files found in %s. Place your input files there, or pass ",
      input_dir
    ), "zip_path= to unzip one, then re-run.")
  }

  already_mapped <- unlist(state$typemap, use.names = FALSE)
  files <- files[!file.path(input_dir, files) %in% already_mapped]

  if (length(files) > 0) wb_list_data_categories()

  for (f in files) {
    choice <- wb_smart_readline(
      sprintf("  %s -> category index: ", f),
      valid = function(ch) {
        n <- suppressWarnings(as.integer(ch))
        if (is.na(n) || n < 0 || n > length(CAT_MAP)) {
          sprintf("Invalid category number (0-%d).", length(CAT_MAP))
        } else TRUE
      }
    )
    if (is.null(choice)) { wb_msg("CANCELLED", "Stopped mapping remaining files."); break }
    choice <- as.integer(choice)
    if (choice == 0) next
    state$typemap[[CAT_MAP[choice]]] <- file.path(input_dir, f)
  }

  if (is.null(state$typemap$ptmseaDB)) {
    state$typemap$ptmseaDB <- wb_default_asset("^ptm\\.sig\\.db\\.all\\.flanking\\.human.*\\.gmt$")
  }
  if (is.null(state$typemap$gseaDB)) {
    state$typemap$gseaDB <- wb_default_asset("^h\\.all.*\\.symbols\\.gmt$")
  }

  wb_save_state(state)
}

wb_validate_annotation_table <- function(annot) {
  missing_cols <- setdiff(REQUIRED_COLS, colnames(annot))
  if (length(missing_cols) > 0) {
    stop("Annotation file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (any(duplicated(annot$Sample.ID))) {
    stop("Sample.ID values in the annotation file are not unique.")
  }
  wb_msg("INFO", "Sample annotation file successfully validated.")
  invisible(TRUE)
}

wb_validate_sample_ids <- function(data_type, data_ids, annot_ids) {
  common <- intersect(data_ids, annot_ids)
  if (length(common) == 0) {
    stop(sprintf("No overlapping Sample.IDs found between '%s' and the annotation table.", data_type))
  }
  if (length(common) < length(annot_ids) / 2) {
    wb_msg("WARNING", sprintf("Only %d of %d annotation samples found in '%s'.",
                              length(common), length(annot_ids), data_type))
  } else {
    wb_msg("INFO", sprintf("'%s' successfully validated (%d/%d samples matched).",
                           data_type, length(common), length(annot_ids)))
  }
  invisible(TRUE)
}

wb_validate_gene_id_column <- function(gct, gct_path, ome, params) {
  if (ome == "metabolome") {
    wb_msg("INFO", "Skipping gene-ID column check for METABOLOME data.")
    return(invisible(NULL))
  }

  gene_id_col_default <- params$global_parameters$gene_mapping$gene_id_col
  protein_id_col       <- params$global_parameters$gene_mapping$protein_id_col
  protein_id_type      <- params$global_parameters$gene_mapping$protein_id_type

  rdesc_names <- colnames(gct@rdesc)

  valid <- FALSE
  if (gene_id_col_default %in% rdesc_names) {
    valid <- tryCatch({
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("org.Hs.eg.db not installed")
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      ids <- as.character(stats::na.omit(unique(gct@rdesc[[gene_id_col_default]])))
      nrow(AnnotationDbi::select(org.Hs.eg.db, keys = ids, keytype = "SYMBOL", columns = "SYMBOL")) > 0
    }, error = function(e) FALSE)
    if (valid) {
      wb_msg("INFO", sprintf("Default gene-ID column '%s' detected and valid in %s data.", gene_id_col_default, toupper(ome)))
    } else {
      wb_msg("WARNING", sprintf("Column '%s' in %s data does not contain valid HUGO gene symbols.", gene_id_col_default, toupper(ome)))
    }
  } else {
    wb_msg("WARNING", sprintf("Default gene-ID column '%s' not found in %s data.", gene_id_col_default, toupper(ome)))
  }
  if (valid) return(invisible(gct_path))

  repeat {
    cat(sprintf("\n%s row-annotation columns: %s\n", toupper(ome), paste(rdesc_names, collapse = ", ")))
    flush.console()
    choice <- wb_smart_readline(
      paste0("To create a Gene ID column for ", toupper(ome), ", choose:\n",
             "  1) Select an existing column with HUGO gene symbols\n",
             "  2) Convert a protein-ID column ('", protein_id_col, "', type ", protein_id_type, ") to HUGO gene symbols\n",
             "  3) Skip (proceed without a validated gene-ID column)\n> "),
      valid = function(ch) if (ch %in% c("1", "2", "3")) TRUE else "Please enter 1, 2, or 3 (or 'quit' to cancel)."
    )
    if (is.null(choice) || choice == "3") {
      wb_msg("WARNING", sprintf("Skipping gene-ID validation for %s data. Many PANOPLY modules require this column.", toupper(ome)))
      break
    }
    if (choice == "1") {
      col <- wb_smart_readline("Column with HUGO gene symbols: ",
                               valid = function(ch) if (ch %in% rdesc_names) TRUE else "Column not found, try again.")
      if (is.null(col)) next
      gct@rdesc[[gene_id_col_default]] <- gct@rdesc[[col]]
      cmapR::write_gct(gct, gct_path, appenddim = FALSE)
      wb_msg("INFO", sprintf("Using column '%s' as '%s' for %s data.", col, gene_id_col_default, toupper(ome)))
      break
    } else {
      if (!(protein_id_col %in% rdesc_names)) {
        stop(sprintf("Configured protein-ID column '%s' not found in %s data.", protein_id_col, toupper(ome)))
      }
      if (!exists("map_id")) {
        stop("map_id() is unavailable (vendored proteomics-Rutil scripts failed to load) -- cannot convert protein IDs.")
      }
      gct@rdesc[[gene_id_col_default]] <- map_id(gct@rdesc[[protein_id_col]], keytype_from = protein_id_type, keytype_to = "SYMBOL")
      cmapR::write_gct(gct, gct_path, appenddim = FALSE)
      wb_msg("INFO", sprintf("Converted '%s' (%s) to gene symbols for %s data.", protein_id_col, protein_id_type, toupper(ome)))
      break
    }
  }
  invisible(gct_path)
}

wb_validate_flanking_sequence_column <- function(gct_path, params) {
  seqwin_default <- params$panoply_preprocess_gct$seqwin_column
  gct <- cmapR::parse_gctx(gct_path)
  rdesc_names <- colnames(gct@rdesc)
  pattern <- "[A-Za-z-]{7}[sty][A-Za-z-]{7}"

  valid <- seqwin_default %in% rdesc_names && any(grepl(pattern, gct@rdesc[[seqwin_default]]))
  if (valid) {
    wb_msg("INFO", sprintf("Default flanking-sequence column '%s' detected and valid.", seqwin_default))
    return(invisible(gct_path))
  }
  wb_msg("WARNING", sprintf("Default flanking-sequence column '%s' missing or invalid.", seqwin_default))
  cat(sprintf("PHOSPHOPROTEOME columns: %s\n", paste(rdesc_names, collapse = ", ")))
  flush.console()
  col <- wb_smart_readline(
    "Column with flanking sequences: ",
    valid = function(ch) {
      if (!(ch %in% rdesc_names)) return("Column not found, try again.")
      if (!any(grepl(pattern, gct@rdesc[[ch]]))) return("No valid flanking sequences found in that column, try again.")
      TRUE
    }
  )
  if (is.null(col)) {
    wb_msg("WARNING", "Skipped flanking-sequence setup. PTM-SEA requires this column.")
    return(invisible(gct_path))
  }
  gct@rdesc[[seqwin_default]] <- gct@rdesc[[col]]
  cmapR::write_gct(gct, gct_path, appenddim = FALSE)
  wb_msg("INFO", sprintf("Using column '%s' as '%s'.", col, seqwin_default))
  invisible(gct_path)
}

wb_metab_compound_db_path <- function(github_ref = GITHUB_REF) {
  # Same staged-first, live-fetch-as-backup pattern as wb_load_default_master_parameters():
  # the staged copy (populated by deploy-workbench.sh directly from the canonical
  # src/panoply_metaboanalyst/ right before upload) is the stable, version-pinned copy that
  # shipped with this deployment. Dev-mode fallback for running straight out of a full repo
  # checkout. Only if neither is present do we live-fetch from GitHub as a backup -- with a
  # warning, since a live fetch may not match the version actually pinned at deployment.
  staged_candidates <- c(
    file.path("workbench-src", "defaults", "master_compound_db.qs"),
    file.path("..", "src", "panoply_metaboanalyst", "pathway_db", "master_compound_db.qs")
  )
  staged <- staged_candidates[file.exists(staged_candidates)][1]
  if (!is.na(staged)) return(staged)

  wb_msg("WARNING", sprintf(
    paste("No staged master_compound_db.qs found (expected at %s); live-fetching from",
          "GitHub (%s@%s) as a backup -- this may not match the stable, version-pinned copy",
          "normally staged when the workbench is deployed."),
    paste(staged_candidates, collapse = " or "), GITHUB_REPO, github_ref
  ))
  out_path <- tempfile(fileext = ".qs")
  tryCatch({
    wb_gh_fetch_binary("src/panoply_metaboanalyst/pathway_db/master_compound_db.qs", out_path, ref = github_ref)
    out_path
  }, error = function(e) {
    stop("No staged master_compound_db.qs found, and the live GitHub fetch failed too (",
         conditionMessage(e), ").")
  })
}

wb_validate_metabolite_id_column <- function(gct_path, params, github_ref = GITHUB_REF) {
  metab_id_col_default  <- params$panoply_metaboanalyst$meta_id_col
  metab_id_type_default <- params$panoply_metaboanalyst$meta_id_type

  gct <- cmapR::parse_gctx(gct_path)
  rdesc_names <- colnames(gct@rdesc)
  compound_map <- qs::qread(wb_metab_compound_db_path(github_ref))

  if (!(metab_id_type_default %in% names(compound_map))) {
    stop(sprintf("Configured metabolite ID type '%s' is not one of the supported types (%s).",
                 metab_id_type_default, paste(names(compound_map), collapse = ", ")))
  }

  valid <- metab_id_col_default %in% rdesc_names &&
    any(!is.na(gct@rdesc[[metab_id_col_default]]) &
        gct@rdesc[[metab_id_col_default]] %in% compound_map[[metab_id_type_default]])
  if (valid) {
    wb_msg("INFO", sprintf("Default metabolite-ID column '%s' detected and valid.", metab_id_col_default))
    return(invisible(gct_path))
  }
  wb_msg("WARNING", sprintf("Default metabolite-ID column '%s' missing or invalid.", metab_id_col_default))
  cat(sprintf("METABOLOME columns: %s (or '0' to use GCT row IDs)\n", paste(rdesc_names, collapse = ", ")))
  flush.console()
  repeat {
    col <- wb_smart_readline(
      "Column with metabolite IDs (or 0 for row IDs): ",
      valid = function(ch) if (identical(ch, "0") || ch %in% rdesc_names) TRUE else "Column not found, try again."
    )
    if (is.null(col)) {
      wb_msg("WARNING", "Skipped metabolite-ID setup. panoply_metaboanalyst requires this column.")
      return(invisible(gct_path))
    }
    if (identical(col, "0")) { ids <- gct@rid; col_label <- "rid" }
    else { ids <- gct@rdesc[[col]]; col_label <- col }

    id_type <- metab_id_type_default
    if (!wb_confirm(sprintf("Does '%s' use %s IDs?", col_label, metab_id_type_default))) {
      cat(sprintf("Supported ID types: %s\n", paste(names(compound_map), collapse = ", ")))
      flush.console()
      id_type <- wb_smart_readline(
        "ID type: ",
        valid = function(ch) if (ch %in% names(compound_map)) TRUE else "Unsupported ID type, try again."
      )
      if (is.null(id_type)) {
        wb_msg("WARNING", "Skipped metabolite-ID setup. panoply_metaboanalyst requires this column.")
        return(invisible(gct_path))
      }
    }
    if (!any(!is.na(ids) & ids %in% compound_map[[id_type]])) { wb_msg("WARNING", "No valid IDs found, try again."); next }

    gct@rdesc[[metab_id_col_default]] <- if (identical(col, "0")) {
      gct@rid
    } else if (id_type != metab_id_type_default) {
      compound_map[[metab_id_type_default]][match(ids, compound_map[[id_type]])]
    } else {
      ids
    }
    cmapR::write_gct(gct, gct_path, appenddim = FALSE)
    wb_msg("INFO", sprintf("Using column '%s' (%s) as '%s'.", col_label, id_type, metab_id_col_default))
    break
  }
  invisible(gct_path)
}

wb_validate_inputs <- function(state) {
  if (is.null(state$typemap$annotation)) stop("No annotation file mapped -- run wb_load_and_map_inputs() first.")
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  wb_validate_annotation_table(annot)

  gct_categories <- intersect(names(state$typemap), c(PROTEOME_TYPES, "rna", "cna", "metabolome"))
  if (length(gct_categories) == 0) stop("No GCT files mapped -- run wb_load_and_map_inputs() first.")
  if (length(intersect(gct_categories, PROTEOME_TYPES)) == 0) stop("No proteomics dataset mapped.")

  params <- if (!is.null(state$typemap$parameters)) yaml::read_yaml(state$typemap$parameters) else wb_load_default_master_parameters(state$github_ref)

  wb_msg("INFO", "Validating sample IDs and gene-ID columns in GCT files...")
  for (cat_name in gct_categories) {
    gct_path <- state$typemap[[cat_name]]
    gct <- cmapR::parse_gctx(gct_path)
    wb_validate_sample_ids(cat_name, gct@cid, annot$Sample.ID)
    wb_validate_gene_id_column(gct, gct_path, cat_name, params)
  }

  wb_save_state(state)
}

wb_select_preprocessing_options <- function(state) {
  state$toggles$normalize_proteomics <- wb_confirm("Does proteomics data need normalization?")
  state$toggles$filter_proteomics    <- wb_confirm("Does proteomics data need filtering?")

  params <- if (!is.null(state$typemap$parameters)) yaml::read_yaml(state$typemap$parameters) else wb_load_default_master_parameters(state$github_ref)

  state$toggles$run_ptmsea <- FALSE
  if (!is.null(state$typemap$phosphoproteome) && !is.null(state$typemap$ptmseaDB)) {
    state$toggles$run_ptmsea <- wb_confirm("Phosphoproteome data detected. Should PTM-SEA be run?")
    if (state$toggles$run_ptmsea) wb_validate_flanking_sequence_column(state$typemap$phosphoproteome, params)
  }

  state$toggles$run_metab <- FALSE
  if (!is.null(state$typemap$metabolome)) {
    state$toggles$run_metab <- wb_confirm("Metabolomics data detected. Should MetaboAnalyst be run?")
    if (state$toggles$run_metab) wb_validate_metabolite_id_column(state$typemap$metabolome, params, state$github_ref)
  }

  wb_save_state(state)
}
