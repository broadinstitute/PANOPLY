# Input upload, category mapping, and validation.

# cmapR::write_gct() works poorly in server settings; write-times balloon out of control
# due to inefficient file-writing, and partial writes can result in corrupted files.
# To avoid this, the file is instead written to a temporary local directory, then copied
# to its final location.
wb_write_gct_atomic <- function(gct, path, quiet = FALSE) {
  # quiet = TRUE lets a caller that already printed its own more specific status line (e.g.
  # wb_write_subset(), writing several GCTs in a batch) skip this generic one.
  if (!quiet) wb_msg("INFO", "Writing changes -- this can take a while for large files, please wait...")
  local_tmp <- tempfile(fileext = ".gct")
  on.exit(unlink(local_tmp), add = TRUE) # clear temporary file on exit
  # invisible(capture.output()) used to silence cmapR::write_gct() printouts,
  # since the user doesn't need to see the temporary file-path.
  # note: genuine failure from write_gct() itself still surfaces normally.
  invisible(capture.output(cmapR::write_gct(gct, local_tmp, appenddim = FALSE)))
  tmp_path <- tempfile(tmpdir = dirname(path), fileext = ".gct")
  on.exit(unlink(tmp_path), add = TRUE) # clear temporary file on exit
  if (!file.copy(local_tmp, tmp_path)) { # attempt to copy file to server
    stop(sprintf("Failed to copy the written file into '%s'.", dirname(path)))
  }
  if (!file.rename(tmp_path, path)) { # rename file to permanent path (automatically overwrites)
    stop(sprintf("Failed to move the written file into place at '%s'.", path))
  }
  invisible(path)
}

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
  if (!file.exists(seeded_path) && !file.copy(source_path, seeded_path)) {
    stop(sprintf("Failed to seed default asset '%s' to '%s'.", source_path, seeded_path))
  }
  seeded_path
}

# Resolves a user-typed FASTA path to a local filesystem path. s3:// URIs are the only case
# that need special handling -- translated via wb_s3_to_local() (the inverse of
# wb_local_to_s3()); NA is returned if that can't resolve it (different project, malformed,
# S3_BUCKET/PROJECT_ID unset), so the caller can give a clear "use a local path instead"
# message rather than silently treating "s3://..." itself as a (nonexistent) local path.
# Anything else is already a normal local path -- path.expand() handles "~", and relative/
# absolute paths need no further wrangling (file.exists() resolves them against the cwd as-is).
wb_resolve_fasta_path <- function(p) {
  if (startsWith(p, "s3://")) {
    local <- wb_s3_to_local(p)
    return(if (is.null(local)) NA_character_ else local)
  }
  path.expand(p)
}

wb_validate_fasta_input <- function(raw) {
  resolved <- wb_resolve_fasta_path(raw)
  if (is.na(resolved)) {
    return(paste(
      "Cloud (s3://) paths aren't supported directly here -- please upload the FASTA under",
      "~/workbench/ and enter its local path instead, try again."
    ))
  }
  if (!file.exists(resolved)) return(sprintf("No file found at '%s', try again.", resolved))
  if (!grepl("\\.(fasta|fa)$", resolved, ignore.case = TRUE)) {
    return(sprintf("'%s' doesn't look like a .fasta/.fa file, try again.", resolved))
  }
  TRUE
}

# Suggests a CAT_MAP category for a filename by substring match (e.g. "ODG-v4-proteome-....gct"
# -> "proteome"), checking the LONGEST category names first -- otherwise a file matching
# "phosphoproteome" or "nglycoproteome" would incorrectly match the shorter "proteome" substring
# they both also contain. Just a suggestion (see wb_load_and_map_inputs()): the user can always
# accept or override it. Returns NA if nothing matches.
wb_smart_match_category <- function(filename) {
  ordered <- CAT_MAP[order(nchar(CAT_MAP), decreasing = TRUE)]
  filename_lower <- tolower(filename)
  # grepl() doesn't support ignore.case=TRUE together with fixed=TRUE (it's silently ignored,
  # with a warning) -- lowercase both sides instead for a case-insensitive literal match.
  hit <- ordered[vapply(ordered, function(cat) grepl(tolower(cat), filename_lower, fixed = TRUE), logical(1))]
  if (length(hit) == 0) return(NA_integer_)
  match(hit[1], CAT_MAP)
}

wb_load_and_map_inputs <- function(state, input_dir = file.path(wb_workbench_root(), "inputs"),
                                    zip_path = NULL) {
  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(zip_path)) {
    # Auto-detect rather than requiring zip_path= to be passed -- but still just a suggestion:
    # declining (or there being no zip at all) falls through to sorting input_dir's loose files
    # exactly as before.
    zips <- list.files(input_dir, pattern = "\\.zip$", ignore.case = TRUE, full.names = TRUE)
    if (length(zips) == 1) {
      if (wb_confirm(sprintf("Found '%s' -- unzip and use its contents?", basename(zips)))) zip_path <- zips
    } else if (length(zips) > 1) {
      cat(sprintf("Found multiple zip files in %s:\n", input_dir))
      for (i in seq_along(zips)) cat(sprintf("  %d: %s\n", i, basename(zips[i])))
      flush.console()
      idx <- wb_smart_readline(
        "Unzip and use one of these? Enter its number, or leave blank to skip all of them: ",
        allow_empty = TRUE,
        valid = function(ch) {
          if (!nzchar(ch)) return(TRUE)
          n <- suppressWarnings(as.integer(ch))
          if (is.na(n) || n < 1 || n > length(zips)) sprintf("Enter a number from 1 to %d.", length(zips)) else TRUE
        }
      )
      if (!is.null(idx) && nzchar(idx)) zip_path <- zips[as.integer(idx)]
    }
  } else {
    zip_path <- path.expand(zip_path)
    # unzip()'s own extraction backend doesn't distinguish "file not found" from "corrupt
    # archive" -- both surface as the same generic "error 1 in extracting from zip file"
    # warning, which makes a simple typo'd path look identical to a genuinely broken zip.
    # Checking existence up front makes that failure mode unambiguous.
    if (!file.exists(zip_path)) stop(sprintf("zip_path '%s' does not exist.", zip_path))
  }

  # Files that actually came from the zip, if one's being used -- so the "want to also map
  # other files sitting in this folder?" question below can tell them apart from whatever else
  # already happened to be in input_dir, rather than assuming both are wanted together.
  zip_files <- character(0)
  used_zip <- !is.null(zip_path)
  if (used_zip) zip_files <- basename(utils::unzip(zip_path, exdir = input_dir, junkpaths = TRUE))

  files <- list.files(input_dir, pattern = "\\.(gct|csv|ya?ml|gmt|fasta|fa)$", full.names = FALSE)
  if (length(files) == 0) {
    stop(sprintf(
      "No .gct/.csv/.yaml/.gmt/.fasta/.fa files found in %s (a .zip there is detected automatically). ",
      input_dir
    ), "Place your input files there and re-run.")
  }

  # Dedup against the ORIGINAL upload paths already consumed, not state$typemap -- typemap
  # holds the session-local COPY path (see below), which would never match input_dir's
  # paths and would cause every already-mapped file to be re-offered on every re-run.
  already_mapped <- unlist(state$typemap_originals, use.names = FALSE)
  files <- files[!file.path(input_dir, files) %in% already_mapped]

  if (used_zip) {
    other_files <- setdiff(files, zip_files)
    files <- intersect(files, zip_files)
    if (length(other_files) > 0 &&
        wb_confirm(sprintf("%d other file(s) in this folder weren't part of the zip -- map those too?", length(other_files)))) {
      files <- c(files, other_files)
    }
  }

  if (length(files) > 0) wb_list_data_categories()

  skipped_files <- character(0)
  cancelled <- FALSE

  for (f in files) {
    suggested <- wb_smart_match_category(f)
    prompt <- if (!is.na(suggested)) {
      sprintf("  %s -> category index [detected: %d) %s -- press Enter to accept, or type a different number]: ",
              f, suggested, CAT_MAP[suggested])
    } else {
      sprintf("  %s -> category index: ", f)
    }
    choice <- wb_smart_readline(
      prompt,
      allow_empty = !is.na(suggested),
      valid = function(ch) {
        if (!nzchar(ch)) return(TRUE)  # accepts the smart-detected suggestion, if any
        n <- suppressWarnings(as.integer(ch))
        if (is.na(n) || n < 0 || n > length(CAT_MAP)) {
          sprintf("Invalid category number (0-%d).", length(CAT_MAP))
        } else TRUE
      },
      cancel_msg = "Stopped mapping remaining files -- files already mapped are saved."
    )
    if (is.null(choice)) { cancelled <- TRUE; break }
    choice <- if (!nzchar(choice)) suggested else as.integer(choice)

    if (choice == 0) {
      if (!grepl("\\.gct$", f, ignore.case = TRUE)) {
        # Only GCTs get the "register as a new -ome" option below -- there's no equivalent
        # concept for an annotation/groups/parameters/database file that doesn't fit a category.
        skipped_files <- c(skipped_files, f)
        next
      }
      new_ome <- wb_smart_readline(
        paste(
          "If you would like to map this GCT to a new -ome, enter a name to register",
          "(it'll be subsetted and processed like other -omes, but NOT automatically wired into inputs.json),",
          "or leave blank to skip this file: "
        ),
        allow_empty = TRUE,
        valid = function(ch) {
          if (!nzchar(ch)) return(TRUE)
          if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", ch)) {
            return("Use a name starting with a letter (letters, numbers, underscore only), try again.")
          }
          if (ch %in% CAT_MAP) {
            return(sprintf("'%s' is already a standard category -- pick it by its index above instead, try again.", ch))
          }
          if (ch %in% names(state$typemap)) {
            return(sprintf("'%s' is already used for another mapped file, try again.", ch))
          }
          TRUE
        },
        cancel_msg = "Stopped mapping remaining files -- files already mapped are saved."
      )
      if (is.null(new_ome)) { cancelled <- TRUE; break }
      if (!nzchar(new_ome)) { skipped_files <- c(skipped_files, f); next }
      original_path <- file.path(input_dir, f)
      state$typemap[[new_ome]] <- wb_copy_into_session(original_path)
      state$typemap_originals[[new_ome]] <- original_path
      next
    }

    original_path <- file.path(input_dir, f)
    # Copy into the session rather than pointing typemap at the original upload directly --
    # downstream validators write in-place fixes (gene-ID column, etc.) to whatever path they're
    # given, so mapping the copy here is what keeps the original in ~/workbench/inputs/ untouched.
    state$typemap[[CAT_MAP[choice]]] <- wb_copy_into_session(original_path)
    state$typemap_originals[[CAT_MAP[choice]]] <- original_path
  }
  if (!cancelled && length(files) > 0) wb_msg("INFO", "All files have been sorted.")

  if (is.null(state$typemap$ptmseaDB)) {
    state$typemap$ptmseaDB <- wb_default_asset("^ptm\\.sig\\.db\\.all\\.flanking\\.human.*\\.gmt$")
  }
  if (is.null(state$typemap$gseaDB)) {
    state$typemap$gseaDB <- wb_default_asset("^h\\.all.*\\.symbols\\.gmt$")
  }

  cat("\nCurrent file mappings:\n")
  for (cat_name in intersect(CAT_MAP, names(state$typemap))) {
    cat(sprintf("  %-16s -> %s\n", cat_name, wb_display_path(state$typemap[[cat_name]])))
  }

  extra_omes <- setdiff(names(state$typemap), CAT_MAP)
  if (length(extra_omes) > 0) {
    cat("\nAdditional -omes (available for subsetting; not wired into inputs.json automatically):\n")
    for (cat_name in extra_omes) {
      cat(sprintf("  %-16s -> %s\n", cat_name, wb_display_path(state$typemap[[cat_name]])))
    }
  }

  if (length(skipped_files) > 0) {
    cat("\nSkipped (not mapped to anything):\n")
    for (f in skipped_files) cat(sprintf("  %s\n", f))
  }
  flush.console()

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

# A single lucky match is easy to hit by chance in a large dataset, so validity is judged by
# match RATE, not mere existence of one hit -- callers compare this against
# GENE_SYMBOL_MATCH_THRESHOLD. Any failure along the way (package missing, unexpected data,
# etc.) returns 0 rather than raising an error, which fails any positive threshold the same
# way a genuine zero match rate would.
GENE_SYMBOL_MATCH_THRESHOLD <- 0.10

wb_gene_symbol_match_rate <- function(values) {
  tryCatch({
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("org.Hs.eg.db not installed")
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    ids <- as.character(stats::na.omit(unique(values)))
    if (length(ids) == 0) return(0)
    # columns = "SYMBOL" here would be a degenerate self-lookup that just echoes every key
    # straight back regardless of whether it's a real symbol (AnnotationDbi::select() only
    # actually validates keys against a genuinely different target column) -- ENTREZID is
    # that real target; a row with a non-NA ENTREZID means the key resolved to an actual gene.
    result <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = ids, keytype = "SYMBOL", columns = "ENTREZID"))
    matched <- unique(result$SYMBOL[!is.na(result$ENTREZID)])
    length(matched) / length(ids)
  }, error = function(e) 0)
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
    match_rate <- wb_gene_symbol_match_rate(gct@rdesc[[gene_id_col_default]])
    valid <- match_rate >= GENE_SYMBOL_MATCH_THRESHOLD
    if (valid) {
      wb_msg("INFO", sprintf("Default gene-ID column '%s' detected and valid (%d%%) in %s data.",
                              gene_id_col_default, round(match_rate * 100), toupper(ome)))
    } else {
      wb_msg("WARNING", sprintf("Column '%s' in %s data does not contain valid HUGO gene symbols (%d%% valid).",
                                 gene_id_col_default, toupper(ome), round(match_rate * 100)))
    }
  } else {
    wb_msg("WARNING", sprintf("Default gene-ID column '%s' not found in %s data.", gene_id_col_default, toupper(ome)))
  }
  if (valid) return(invisible(gct_path))

  cat(sprintf("\n%s row-annotation columns: %s\n\n", toupper(ome), paste(rdesc_names, collapse = ", ")))
  flush.console()
  repeat {
    choice <- wb_smart_readline(
      paste0("To create a Gene ID column for ", toupper(ome), ", choose:\n",
             "  1) Select an existing column with HUGO gene symbols\n",
             "  2) Convert a protein-ID column ('", protein_id_col, "', type ", protein_id_type, ") to HUGO gene symbols\n",
             "  3) Skip (proceed without a validated gene-ID column)\n> "),
      valid = function(ch) if (ch %in% c("1", "2", "3")) TRUE else "Please enter 1, 2, or 3 (or 'quit' to cancel)."
    )
    if (is.null(choice) || choice == "3") {
      wb_msg("WARNING", sprintf("Skipping gene-ID validation for %s data. WARNING: Many PANOPLY modules require this column.", toupper(ome)))
      break
    }
    if (choice == "1") {
      col <- wb_smart_readline(
        "Column with HUGO gene symbols: ",
        valid = function(ch) {
          if (!(ch %in% rdesc_names)) return("Column not found, try again.")
          rate <- wb_gene_symbol_match_rate(gct@rdesc[[ch]])
          if (rate < GENE_SYMBOL_MATCH_THRESHOLD) {
            return(sprintf("'%s' does not appear to contain valid HUGO gene symbols (%d%% valid), try again.",
                           ch, round(rate * 100)))
          }
          TRUE
        }
      )
      if (is.null(col)) next
      # report successful mapping-- not just failure-- to confirm to the user that provided column is valid
      match_rate <- wb_gene_symbol_match_rate(gct@rdesc[[col]])
      wb_msg("INFO", sprintf("Column '%s' validated (%d%% matched HUGO gene symbols).", col, round(match_rate * 100)))
      gct@rdesc[[gene_id_col_default]] <- gct@rdesc[[col]]
      wb_msg("INFO", sprintf("Using column '%s' as '%s' for %s data.", col, gene_id_col_default, toupper(ome)))
      wb_write_gct_atomic(gct, gct_path)
      break
    } else {
      if (!(protein_id_col %in% rdesc_names)) {
        stop(sprintf("Configured protein-ID column '%s' not found in %s data.", protein_id_col, toupper(ome)))
      }
      if (!exists("map_id")) {
        stop("map_id() is unavailable (vendored proteomics-Rutil scripts failed to load) -- cannot convert protein IDs.")
      }
      converted <- map_id(gct@rdesc[[protein_id_col]], keytype_from = protein_id_type, keytype_to = "SYMBOL")
      match_rate <- wb_gene_symbol_match_rate(converted)
      if (match_rate < GENE_SYMBOL_MATCH_THRESHOLD) {
        wb_msg("WARNING", sprintf(
          "Converting '%s' (%s) produced mostly invalid gene symbols (%d%% matched) -- not using it.",
          protein_id_col, protein_id_type, round(match_rate * 100)
        ))
        next
      }
      wb_msg("INFO", sprintf("Conversion validated (%d%% matched HUGO gene symbols).", round(match_rate * 100)))
      gct@rdesc[[gene_id_col_default]] <- converted
      wb_msg("INFO", sprintf("Converted '%s' (%s) to gene symbols for %s data.", protein_id_col, protein_id_type, toupper(ome)))
      wb_write_gct_atomic(gct, gct_path)
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
  wb_msg("INFO", sprintf("Using column '%s' as '%s'.", col, seqwin_default))
  wb_write_gct_atomic(gct, gct_path)
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
    wb_msg("INFO", sprintf("Using column '%s' (%s) as '%s'.", col_label, id_type, metab_id_col_default))
    wb_write_gct_atomic(gct, gct_path)
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
