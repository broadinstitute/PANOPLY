# Local-folder sample subsetting -- replaces Terra sample sets.

# The set of mapped source files a subset's contents actually depend on -- every GCT/CSV
# category wb_write_subset() below reads from.
wb_subset_source_paths <- function(state) {
  gct_categories <- intersect(names(state$typemap), c(PROTEOME_TYPES, "rna", "cna", "metabolome"))
  csv_categories <- intersect(names(state$typemap), c("annotation", "groups", "groups_clumpsptm"))
  unlist(state$typemap[c(gct_categories, csv_categories)])
}

# A cheap staleness fingerprint: mtimes of every mapped source file, keyed by category. Cheap
# because it's metadata only -- no need to read a potentially huge GCT to know whether it
# changed. Comparing this (via identical()) against what was recorded the last time a subset
# was written detects both in-place edits (e.g. wb_validate_gene_id_column() rewriting a GCT --
# a fresh atomic write always bumps mtime, even if the write happens to reproduce the same
# bytes) and typemap changes (a category added/removed/remapped to a different file changes
# which paths are even in the fingerprint).
wb_subset_source_mtimes <- function(state) {
  paths <- wb_subset_source_paths(state)
  if (length(paths) == 0) return(list())
  setNames(as.list(as.character(file.mtime(paths))), names(paths))
}

# Does the actual file-writing for one subset -- no prompts, no wb_save_state()/wb_done()
# (the interactive wb_create_subset() below is the only public entry point, and saves once
# after everything it creates in one call, not once per subset). Skips the (potentially slow,
# for large GCTs) rewrite entirely if this exact subset was already built from these exact
# source files -- re-running wb_create_subset() is then cheap and idempotent, while an actual
# upstream change (a re-mapped input, an in-place GCT fix) still triggers a real rebuild.
wb_write_subset <- function(state, name, filter_col = NULL, filter_vals = NULL,
                             out_root = file.path(wb_session_dir(), "subsets"), force = FALSE) {
  source_mtimes <- wb_subset_source_mtimes(state)
  existing <- state$subsets[[name]]
  if (!force && !is.null(existing) &&
      identical(existing$filter_col, filter_col) &&
      identical(existing$filter_vals, filter_vals) &&
      identical(existing$source_mtimes, source_mtimes) &&
      dir.exists(existing$dir)) {
    wb_msg("INFO", sprintf("Subset '%s' is already up to date; skipping.", name))
    return(state)
  }

  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  sample_ids <- if (is.null(filter_col)) {
    annot$Sample.ID
  } else {
    if (!(filter_col %in% colnames(annot))) stop(sprintf("'%s' is not a column in the annotation table.", filter_col))
    annot$Sample.ID[annot[[filter_col]] %in% filter_vals]
  }
  if (length(sample_ids) == 0) stop("No samples matched the given filter.")

  subset_dir <- file.path(out_root, name)
  dir.create(subset_dir, showWarnings = FALSE, recursive = TRUE)

  gct_categories <- intersect(names(state$typemap), c(PROTEOME_TYPES, "rna", "cna", "metabolome"))
  for (cat_name in gct_categories) {
    # cmapR::parse_gctx()'s own "parsing as GCT v1.3" printout is silenced here -- the single
    # status line below is all that's needed per GCT, per subset.
    suppressMessages(invisible(capture.output(gct <- cmapR::parse_gctx(state$typemap[[cat_name]]))))
    keep <- intersect(sample_ids, gct@cid)
    if (length(keep) == 0) {
      wb_msg("WARNING", sprintf("No samples from this subset found in '%s'; skipping.", cat_name))
      next
    }
    sub <- cmapR::subset_gct(gct, cid = keep)
    wb_msg("INFO", sprintf("Writing %s GCT for '%s'.", cat_name, name))
    wb_write_gct_atomic(sub, file.path(subset_dir, paste0(cat_name, ".gct")), quiet = TRUE)
  }

  csv_categories <- intersect(names(state$typemap), c("annotation", "groups", "groups_clumpsptm"))
  for (cat_name in csv_categories) {
    csv <- read.csv(state$typemap[[cat_name]], stringsAsFactors = FALSE, quote = '"')
    write.csv(csv[csv$Sample.ID %in% sample_ids, , drop = FALSE],
              file.path(subset_dir, paste0(cat_name, ".csv")), row.names = FALSE, quote = TRUE)
  }

  state$subsets[[name]] <- list(filter_col = filter_col, filter_vals = filter_vals,
                                dir = subset_dir, n_samples = length(sample_ids),
                                source_mtimes = source_mtimes)
  wb_msg("INFO", sprintf("Subset '%s' created with %d sample(s) at %s", name, length(sample_ids), subset_dir))
  state
}

wb_create_subset <- function(state, out_root = file.path(wb_session_dir(), "subsets")) {
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')

  # Every subset (including 'all') is only actually written at the very end, in one batch --
  # writing GCTs is slow, so collecting every name/filter first means the user answers all the
  # prompts up front instead of waiting between each one.
  wb_msg("INFO", "An 'all' subset (every sample) will always be created.")
  requests <- list(list(name = "all", filter_col = NULL, filter_vals = NULL))
  pending_names <- function() vapply(requests, function(r) r$name, character(1))

  filter_cols <- setdiff(colnames(annot), "Sample.ID")
  printed_cols <- FALSE
  repeat {
    if (!wb_confirm("Create an additional subset?")) break

    if (!printed_cols) {
      cat("Annotation columns:\n")
      for (i in seq_along(filter_cols)) cat(sprintf("  %2d: %s\n", i, filter_cols[i]))
      flush.console()
      printed_cols <- TRUE
    }
    col_idx <- wb_smart_readline(
      "Select an annotation column to filter on (or 'quit' to cancel): ",
      valid = function(ch) {
        n <- suppressWarnings(as.integer(ch))
        if (is.na(n) || n < 1 || n > length(filter_cols)) {
          sprintf("Please enter a number from 1 to %d.", length(filter_cols))
        } else TRUE
      }
    )
    if (is.null(col_idx)) { wb_msg("CANCELLED", "No additional subset created."); next }
    filter_col <- filter_cols[as.integer(col_idx)]

    # NA values aren't independently selectable here (same as passing filter_vals = NA to
    # wb_write_subset() directly never matched anything, via %in%'s own NA handling).
    values <- sort(unique(as.character(annot[[filter_col]])))
    values <- values[!is.na(values)]

    cat(sprintf("\n'%s' values:\n", filter_col))
    for (i in seq_along(values)) cat(sprintf("  %2d: %s\n", i, values[i]))
    flush.console()

    sel <- wb_smart_readline(
      "Select value(s) to include -- comma-separated indexes or ranges (e.g. 1,3:5) (or 'quit' to cancel): ",
      valid = function(ch) {
        tokens <- wb_trim(strsplit(ch, ",")[[1]])
        if (!all(grepl("^[0-9]+(:[0-9]+)?$", tokens))) {
          return("Use indexes or ranges only (e.g. 1,3:5), try again.")
        }
        idx <- wb_parse_index_ranges(tokens)
        if (any(idx < 1 | idx > length(values))) {
          return(sprintf("Index out of range (1-%d), try again.", length(values)))
        }
        vals <- values[idx]
        if (sum(annot[[filter_col]] %in% vals) == 0) {
          return(sprintf("No samples match '%s' in {%s}, try again.", filter_col, paste(vals, collapse = ", ")))
        }
        TRUE
      }
    )
    if (is.null(sel)) { wb_msg("CANCELLED", "No additional subset created."); next }
    filter_vals <- values[wb_parse_index_ranges(wb_trim(strsplit(sel, ",")[[1]]))]
    n_matched <- sum(annot[[filter_col]] %in% filter_vals)
    wb_msg("INFO", sprintf("%d sample(s) match '%s' in {%s}.", n_matched, filter_col, paste(filter_vals, collapse = ", ")))

    name <- NULL
    repeat {
      candidate <- wb_smart_readline(
        "Name for this subset: ",
        valid = function(ch) if (grepl("^[A-Za-z0-9_.-]+$", ch)) TRUE else "Use only letters, numbers, '-', '_', and '.', try again."
      )
      if (is.null(candidate)) break
      if (candidate %in% union(names(state$subsets), pending_names())) {
        if (wb_confirm(sprintf("A subset named '%s' already exists (or is already queued). Overwrite/replace it?", candidate))) {
          name <- candidate
          break
        }
        # else: loop back and ask for a different name
      } else {
        name <- candidate
        break
      }
    }
    if (is.null(name)) { wb_msg("CANCELLED", "No additional subset created."); next }

    # Replace any earlier queued request with the same name -- only the latest spec is built.
    requests <- Filter(function(r) r$name != name, requests)
    requests[[length(requests) + 1]] <- list(name = name, filter_col = filter_col, filter_vals = filter_vals)
  }

  wb_msg("INFO", sprintf("Creating %d subset(s): %s", length(requests), paste(pending_names(), collapse = ", ")))
  for (req in requests) {
    state <- wb_write_subset(state, req$name, filter_col = req$filter_col, filter_vals = req$filter_vals, out_root = out_root)
  }

  wb_save_state(state)
}

wb_list_subsets <- function(state) names(state$subsets)

wb_subset_files <- function(state, name) {
  if (is.null(state$subsets[[name]])) stop(sprintf("No subset named '%s' -- run wb_create_subset() first.", name))
  files <- list.files(state$subsets[[name]]$dir, full.names = TRUE)
  setNames(as.list(files), sub("\\.(gct|csv)$", "", basename(files)))
}
