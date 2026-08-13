# Group selection, colors, COSMO attributes, ClumpsPTM groups.

wb_list_annotation_columns <- function(state, done = TRUE) {
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  cat("Annotation columns:\n")
  for (col in colnames(annot)) cat(" -", col, "\n")
  flush.console()
  # done = FALSE for internal callers (e.g. wb_select_clumpsptm_groups()) that use this as a
  # mid-function listing, not their final action -- printing "done" here would be premature.
  if (done) wb_done()
  invisible(colnames(annot))
}

wb_verify_group_validity <- function(annot, columns, max_categories) {
  if (length(columns) == 0) return(list(groups_cols = character(0), groups_cols_continuous = character(0)))
  keep <- character(0)
  continuous <- character(0)
  for (col in columns) {
    vals <- unique(annot[[col]])
    n <- length(vals)
    if (n <= 1) {
      wb_msg("WARNING", sprintf("'%s' has %d unique value(s); excluded.", col, n))
    } else if (n > max_categories) {
      if (is.numeric(vals)) {
        continuous <- c(continuous, col)
        wb_msg("INFO", sprintf("'%s' has more than %d unique values; treated as continuous.", col, max_categories))
      } else {
        wb_msg("WARNING", sprintf("'%s' has more than %d unique values and is not numeric; excluded.", col, max_categories))
      }
    } else {
      keep <- c(keep, col)
    }
  }
  list(groups_cols = keep, groups_cols_continuous = continuous)
}

wb_select_groups <- function(state, columns = NULL, max_categories = 10) {
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  all_cols <- colnames(annot)

  if (is.null(columns)) {
    if (!is.null(state$typemap$groups)) {
      columns <- as.vector(unlist(read.csv(state$typemap$groups, stringsAsFactors = FALSE, quote = '"')))
      wb_msg("INFO", "Using columns from the provided groups file.")
    } else {
      columns <- setdiff(all_cols, c(REQUIRED_COLS, IGNORE_COLS))
      wb_msg("INFO", "No groups specified; using all valid annotation columns.")
    }
  }
  columns <- intersect(columns, all_cols)

  validated <- wb_verify_group_validity(annot, columns, max_categories)
  state$groups_cols <- validated$groups_cols
  state$groups_cols_continuous <- validated$groups_cols_continuous

  cat("Selected groups:\n"); for (col in state$groups_cols) cat(" -", col, "\n")
  if (length(state$groups_cols_continuous) > 0) {
    cat("Continuous groups:\n"); for (col in state$groups_cols_continuous) cat(" -", col, "\n")
  }
  flush.console()

  wb_reset_colors(state, annot = annot)
}

wb_assign_colors <- function(annot, groups_cols) {
  if (length(groups_cols) == 0) return(list())
  if (!exists("set_annot_colors")) {
    stop("set_annot_colors() is unavailable (vendored proteomics-Rutil scripts failed to load).")
  }
  raw <- set_annot_colors(annot[, groups_cols, drop = FALSE])
  colors <- list()
  for (group in names(raw)) {
    vals <- ifelse(is.na(raw[[group]]$vals), "NA", as.character(raw[[group]]$vals))
    colors[[group]] <- setNames(as.list(raw[[group]]$colors), vals)
  }
  colors
}

wb_show_colors <- function(state) {
  if (length(state$groups_colors) == 0) {
    wb_msg("WARNING", "No colors assigned yet. Run wb_reset_colors() first.")
    wb_done()
    return(invisible(NULL))
  }
  for (group in names(state$groups_colors)) {
    cat(sprintf("\n%s:\n", group))
    vals <- state$groups_colors[[group]]
    for (v in names(vals)) cat(sprintf("  %-20s %s\n", v, vals[[v]]))
  }
  flush.console()
  wb_done()
  invisible(state$groups_colors)
}

wb_reset_colors <- function(state, annot = NULL) {
  if (is.null(annot)) annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  state$groups_colors <- wb_assign_colors(annot, state$groups_cols)
  wb_save_state(state)
}

wb_edit_color <- function(state, group, value, hex_color) {
  if (!grepl("^#[A-Fa-f0-9]{6}$", hex_color)) stop("hex_color must look like '#RRGGBB'.")
  if (!(group %in% names(state$groups_colors))) stop(sprintf("'%s' is not a known group.", group))
  key <- as.character(value)
  if (!(key %in% names(state$groups_colors[[group]]))) {
    stop(sprintf("'%s' is not a known value for group '%s'.", value, group))
  }
  state$groups_colors[[group]][[key]] <- hex_color
  wb_save_state(state)
}

wb_write_groups_file <- function(annot, groups_cols, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(annot[, c("Sample.ID", groups_cols), drop = FALSE], out_path, row.names = FALSE, quote = TRUE)
  out_path
}

wb_select_cosmo_attributes <- function(state, columns = NULL) {
  if (!wb_confirm("Run COSMO?")) {
    if (!wb_confirm("Select COSMO attributes anyway, to run COSMO later?")) {
      state$cosmo_params <- list(run_cosmo = FALSE, sample_label = "")
      return(wb_save_state(state))
    }
  }

  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  candidates <- setdiff(colnames(annot), IGNORE_COLS)
  valid_attrs <- Filter(function(col) {
    tab <- table(annot[[col]])
    length(tab) == 2 && min(tab) >= min(10, nrow(annot) / 5) && !any(is.na(annot[[col]]))
  }, candidates)

  if (length(valid_attrs) == 0) {
    wb_msg("WARNING", "No valid (binary, balanced, NA-free) attributes found. COSMO will not be run.")
    state$cosmo_params <- list(run_cosmo = FALSE, sample_label = "")
    return(wb_save_state(state))
  }

  cat("Valid COSMO attributes:\n"); for (col in valid_attrs) cat(" -", col, "\n")
  flush.console()
  if (is.null(columns)) {
    columns <- wb_smart_readline(
      "Select attribute(s), comma-separated: ",
      valid = function(ch) {
        picked <- intersect(strsplit(ch, "\\s*,\\s*")[[1]], valid_attrs)
        if (length(picked) == 0) "None of those match a valid attribute above, try again." else TRUE
      }
    )
    columns <- if (is.null(columns)) character(0) else strsplit(columns, "\\s*,\\s*")[[1]]
  }
  columns <- intersect(columns, valid_attrs)

  if (length(columns) == 0) {
    wb_msg("WARNING", "No valid attributes selected. COSMO will not be run.")
    state$cosmo_params <- list(run_cosmo = FALSE, sample_label = "")
  } else {
    state$cosmo_params <- list(run_cosmo = TRUE, sample_label = paste(columns, collapse = ","))
    wb_msg("INFO", sprintf("COSMO will run with attribute(s): %s", paste(columns, collapse = ", ")))
  }
  wb_save_state(state)
}

wb_select_clumpsptm_groups <- function(state, columns = NULL, fasta_path = NULL) {
  ptm_types <- c("phosphoproteome", "acetylome", "ubiquitylome")
  if (length(intersect(ptm_types, names(state$typemap))) < 2) {
    wb_msg("INFO", "Fewer than 2 PTM datasets detected; Clumps-PTM will not be run.")
    state$toggles$run_clumpsptm <- FALSE
    return(wb_save_state(state))
  }

  state$toggles$run_clumpsptm <- wb_confirm("PTM data detected. Should Clumps-PTM be run?")
  if (!state$toggles$run_clumpsptm) return(wb_save_state(state))

  if (is.null(fasta_path) || !file.exists(fasta_path)) {
    stop("Clumps-PTM requires a reference FASTA file -- pass its local path as fasta_path=.")
  }
  state$typemap$clumpsFASTA <- wb_copy_into_session(fasta_path)

  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  if (is.null(columns)) {
    wb_list_annotation_columns(state, done = FALSE)
    columns <- wb_smart_readline(
      "Select up to 3 categorical annotations for Clumps-PTM, comma-separated: ",
      valid = function(ch) {
        picked <- intersect(strsplit(ch, "\\s*,\\s*")[[1]], colnames(annot))
        if (length(picked) == 0) "None of those match a known annotation column, try again." else TRUE
      }
    )
    if (is.null(columns)) {
      wb_msg("WARNING", "Skipped Clumps-PTM annotation selection. Clumps-PTM will not be run.")
      state$toggles$run_clumpsptm <- FALSE
      return(wb_save_state(state))
    }
    columns <- strsplit(columns, "\\s*,\\s*")[[1]]
  }
  columns <- intersect(columns, colnames(annot))
  if (length(columns) > 3) {
    wb_msg("WARNING", "More than 3 annotations selected; Clumps-PTM runs can become long and expensive.")
  }

  out_path <- file.path(dirname(state$typemap$groups %||% state$typemap$annotation), "groups-clumpsptm.csv")
  state$typemap$groups_clumpsptm <- wb_write_groups_file(annot, columns, out_path)
  wb_save_state(state)
}
