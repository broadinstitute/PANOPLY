# Local-folder sample subsetting -- replaces Terra sample sets.

# Every typemap category whose mapped file is a GCT -- covers both the standard proteomics/
# genomics categories AND any custom "extra -ome" registered via wb_load_and_map_inputs()'s
# "register as a new -ome" option (see CAT_MAP there), so custom data gets subsetted the same
# way without needing a fixed category-name whitelist here.
wb_gct_typemap_categories <- function(state) {
  if (length(state$typemap) == 0) return(character(0))
  names(state$typemap)[vapply(state$typemap, function(p) grepl("\\.gct$", p, ignore.case = TRUE), logical(1))]
}

# The set of mapped source files a subset's contents actually depend on -- every GCT/CSV
# category wb_write_subset() below reads from. Deliberately excludes "groups" -- unlike
# annotation/groups_clumpsptm, the per-subset groups.csv is generated fresh from
# annotation+groups_cols (see wb_write_subset()), not read from state$typemap$groups at all, so
# that upload's own mtime isn't a real dependency of the output.
wb_subset_source_paths <- function(state) {
  gct_categories <- wb_gct_typemap_categories(state)
  csv_categories <- intersect(names(state$typemap), c("annotation", "groups_clumpsptm"))
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

# Is an already-built subset stale relative to the CURRENT session data? Same underlying
# comparison wb_write_subset() uses to decide whether to skip a rewrite, minus the
# filter_col/filter_vals check (irrelevant here -- we're asking whether name's own recorded
# filter is still current, not proposing a different one). Used by wb_create_subset() to detect
# and offer to regenerate stale subsets up front, without performing any rebuild itself.
wb_subset_is_stale <- function(state, name) {
  existing <- state$subsets[[name]]
  if (is.null(existing)) return(FALSE)
  !identical(existing$source_mtimes, wb_subset_source_mtimes(state)) ||
    !identical(existing$groups_cols, state$groups_cols) ||
    !dir.exists(existing$dir)
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
      identical(existing$groups_cols, state$groups_cols) &&
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

  gct_categories <- wb_gct_typemap_categories(state)
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

  csv_categories <- intersect(names(state$typemap), c("annotation", "groups_clumpsptm"))
  for (cat_name in csv_categories) {
    csv <- read.csv(state$typemap[[cat_name]], stringsAsFactors = FALSE, quote = '"')
    write.csv(csv[csv$Sample.ID %in% sample_ids, , drop = FALSE],
              file.path(subset_dir, paste0(cat_name, ".csv")), row.names = FALSE, quote = TRUE)
  }

  # groups.csv (the groups_file WDL input) is a per-sample table -- Sample.ID plus the
  # currently-selected group columns' values -- generated fresh from the annotation table each
  # time, mirroring panda-src/bin/r-source/create-groups.r's write_groups(). It is NOT a copy or
  # filter of state$typemap$groups: that upload is just an optional list of column NAMES used to
  # help wb_select_groups() pick state$groups_cols in the first place, and has no Sample.ID
  # column at all -- filtering it as if it were already a per-sample table (the previous
  # behavior) always produced a header with zero data rows.
  wb_write_groups_file(annot[annot$Sample.ID %in% sample_ids, , drop = FALSE],
                       state$groups_cols, file.path(subset_dir, "groups.csv"))

  state$subsets[[name]] <- list(filter_col = filter_col, filter_vals = filter_vals,
                                dir = subset_dir, n_samples = length(sample_ids),
                                source_mtimes = source_mtimes, groups_cols = state$groups_cols)
  wb_msg("INFO", sprintf("Subset '%s' created with %d sample(s) at %s", name, length(sample_ids), subset_dir))
  state
}

wb_describe_subset <- function(name, s) {
  desc <- if (is.null(s$filter_col)) "all samples" else sprintf("%s in {%s}", s$filter_col, paste(s$filter_vals, collapse = ", "))
  sprintf("  %-16s -> %s (%d sample(s))", name, desc, s$n_samples)
}

# Column-index -> value-index(es)-> name prompt sequence for one subset, reused by the "add a
# new subset" menu action below. Returns NULL if cancelled at any point, otherwise
# list(name=, filter_col=, filter_vals=).
wb_prompt_new_subset <- function(state, annot, filter_cols) {
  cat("Annotation columns:\n")
  for (i in seq_along(filter_cols)) cat(sprintf("  %2d: %s\n", i, filter_cols[i]))
  flush.console()
  col_idx <- wb_smart_readline(
    "Select an annotation column to filter on (or 'quit' to cancel): ",
    valid = function(ch) {
      n <- suppressWarnings(as.integer(ch))
      if (is.na(n) || n < 1 || n > length(filter_cols)) sprintf("Please enter a number from 1 to %d.", length(filter_cols)) else TRUE
    },
    cancel_msg = "Previous subset changes saved."
  )
  if (is.null(col_idx)) return(NULL)
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
      if (!all(grepl("^[0-9]+(:[0-9]+)?$", tokens))) return("Use indexes or ranges only (e.g. 1,3:5), try again.")
      idx <- wb_parse_index_ranges(tokens)
      if (any(idx < 1 | idx > length(values))) return(sprintf("Index out of range (1-%d), try again.", length(values)))
      vals <- values[idx]
      if (sum(annot[[filter_col]] %in% vals) == 0) return(sprintf("No samples match '%s' in {%s}, try again.", filter_col, paste(vals, collapse = ", ")))
      TRUE
    },
    cancel_msg = "Previous subset changes saved."
  )
  if (is.null(sel)) return(NULL)
  filter_vals <- values[wb_parse_index_ranges(wb_trim(strsplit(sel, ",")[[1]]))]
  n_matched <- sum(annot[[filter_col]] %in% filter_vals)
  wb_msg("INFO", sprintf("%d sample(s) match '%s' in {%s}.", n_matched, filter_col, paste(filter_vals, collapse = ", ")))

  name <- NULL
  repeat {
    candidate <- wb_smart_readline(
      "Name for this subset: ",
      valid = function(ch) {
        if (!grepl("^[A-Za-z0-9_.-]+$", ch)) return("Use only letters, numbers, '-', '_', and '.', try again.")
        # 'all' always means every sample -- never offer to overwrite it with a custom filter.
        if (identical(ch, "all")) return("'all' is reserved for every sample -- choose a different name, try again.")
        TRUE
      },
      cancel_msg = "Previous subset changes saved."
    )
    if (is.null(candidate)) break
    if (candidate %in% names(state$subsets)) {
      if (wb_confirm(sprintf("A subset named '%s' already exists. Overwrite it?", candidate), cancel_msg = "Previous subset changes saved.")) { name <- candidate; break }
      # else: loop back and ask for a different name
    } else {
      name <- candidate
      break
    }
  }
  if (is.null(name)) return(NULL)
  list(name = name, filter_col = filter_col, filter_vals = filter_vals)
}

# 'all' always exists (created once, the first time there are no subsets yet at all) and can't
# be removed -- beyond that it's treated like any other subset: no automatic re-checking or
# rebuilding on every call, since the staleness-detection prompt below and the explicit
# add/remove/regenerate menu already give full control over when anything actually gets rebuilt.
wb_create_subset <- function(state, out_root = file.path(wb_session_dir(), "subsets")) {
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  filter_cols <- setdiff(colnames(annot), "Sample.ID")

  if (length(state$subsets) == 0) {
    wb_msg("INFO", "An 'all' subset (every sample) will always be created.")
    state <- wb_write_subset(state, "all", out_root = out_root)
  }

  stale <- Filter(function(n) wb_subset_is_stale(state, n), names(state$subsets))
  if (length(stale) > 0 &&
      wb_confirm(sprintf(
        "%d subset(s) look out of date (source data or group columns changed since they were built: %s) -- regenerate them now?",
        length(stale), paste(stale, collapse = ", ")
      ), cancel_msg = "Previous subset changes saved.")) {
    for (n in stale) {
      s <- state$subsets[[n]]
      state <- wb_write_subset(state, n, filter_col = s$filter_col, filter_vals = s$filter_vals, out_root = out_root, force = TRUE)
    }
  }

  repeat {
    if (length(state$subsets) > 0) {
      cat("\nExisting subsets:\n")
      for (n in names(state$subsets)) cat(wb_describe_subset(n, state$subsets[[n]]), "\n")
      flush.console()
    }

    action <- wb_smart_readline(
      paste0(
        "What would you like to do?\n",
        "  1) Add a new subset (or overwrite an existing one by reusing its name)\n",
        "  2) Remove an existing subset (not 'all')\n",
        "  3) Refresh all existing subsets from current session data\n",
        "(or 'quit' to finish): "
      ),
      valid = function(ch) if (ch %in% c("1", "2", "3")) TRUE else "Please enter 1, 2, or 3 (or 'quit' to finish).",
      cancel_msg = "Previous subset changes saved."
    )
    if (is.null(action)) break

    if (action == "1") {
      spec <- wb_prompt_new_subset(state, annot, filter_cols)
      if (is.null(spec)) { wb_msg("CANCELLED", "No subset created."); next }
      state <- wb_write_subset(state, spec$name, filter_col = spec$filter_col, filter_vals = spec$filter_vals, out_root = out_root)
    } else if (action == "2") {
      removable <- setdiff(names(state$subsets), "all")
      if (length(removable) == 0) { wb_msg("WARNING", "Nothing to remove ('all' can't be removed)."); next }
      cat("Removable subsets:\n")
      for (i in seq_along(removable)) cat(sprintf("  %d: %s\n", i, removable[i]))
      flush.console()
      idx <- wb_smart_readline(
        "Remove which subset? Enter its number (or 'quit' to cancel): ",
        valid = function(ch) {
          n <- suppressWarnings(as.integer(ch))
          if (is.na(n) || n < 1 || n > length(removable)) sprintf("Enter a number from 1 to %d.", length(removable)) else TRUE
        },
        cancel_msg = "Previous subset changes saved."
      )
      if (is.null(idx)) next
      target <- removable[as.integer(idx)]
      if (wb_confirm(sprintf("Remove subset '%s'? This deletes its folder (%s) and cannot be undone.",
                             target, state$subsets[[target]]$dir), cancel_msg = "Previous subset changes saved.")) {
        unlink(state$subsets[[target]]$dir, recursive = TRUE)
        state$subsets[[target]] <- NULL
        wb_msg("INFO", sprintf("Removed subset '%s'.", target))
      }
    } else {
      if (length(state$subsets) == 0) { wb_msg("INFO", "No subsets to regenerate yet."); next }
      wb_msg("INFO", sprintf("Regenerating %d subset(s): %s", length(state$subsets), paste(names(state$subsets), collapse = ", ")))
      for (n in names(state$subsets)) {
        s <- state$subsets[[n]]
        state <- wb_write_subset(state, n, filter_col = s$filter_col, filter_vals = s$filter_vals, out_root = out_root, force = TRUE)
      }
    }
  }

  wb_save_state(state)
}

wb_list_subsets <- function(state) names(state$subsets)

# Lists a subset's files -- resolved against the NAMED session (not state$subsets[[name]]$dir's
# raw current-session path) whenever one's active, since this is only ever used to build
# inputs.json, and current-session/'s copy may be stale, empty, or (if a saved session's state
# was opened directly via wb_open_saved_session() rather than fully copied in) never populated
# at all. See wb_in_named_session() in wdl.r.
wb_subset_files <- function(state, name) {
  if (is.null(state$subsets[[name]])) stop(sprintf("No subset named '%s' -- run wb_create_subset() first.", name))
  dir <- state$subsets[[name]]$dir
  if (!is.null(state$active_named_session)) dir <- wb_in_named_session(state, dir)
  files <- list.files(dir, full.names = TRUE)
  setNames(as.list(files), sub("\\.(gct|csv)$", "", basename(files)))
}
