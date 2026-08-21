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

# Clumps-PTM has no "treat as continuous" path -- clumps_diffexp.r just skips any annotation
# whose unique-value count exceeds max_annot_levels outright, regardless of numeric-ness (see
# `if (length(unique(annots[[annot_of_interest]])) > opt$max_annot_levels) next`). Filtering
# here mirrors that exclusion locally, so a column that's certain to be skipped downstream
# doesn't waste a same-cost Cromwell run finding that out. Each retained annotation multiplies
# how many sub-value comparisons Clumps-PTM has to run, hence the cost framing below.
wb_verify_clumpsptm_group_validity <- function(annot, columns, max_categories) {
  Filter(function(col) {
    n <- length(unique(annot[[col]]))
    if (n > max_categories) {
      wb_msg("WARNING", sprintf(
        "'%s' has %d unique values (> %d) -- excluded from Clumps-PTM annotations. Each additional subvalue increases Clumps-PTM's analysis cost.",
        col, n, max_categories
      ))
      FALSE
    } else TRUE
  }, columns)
}

# Parses a comma-separated list of 1-based indexes and/or "start:end" ranges (e.g.
# "1,3:5,8") into a plain integer vector -- mirrors the index/range selection from
# panda-src/build-config.r's select_groups_case(), which offered the same shorthand so users
# didn't have to type out every annotation name by hand.
wb_parse_index_ranges <- function(tokens) {
  unlist(lapply(tokens, function(t) {
    bounds <- as.integer(strsplit(t, ":")[[1]])
    if (length(bounds) == 1) bounds else bounds[1]:bounds[2]
  }))
}

wb_prompt_indexed_columns <- function(valid_cols) {
  cat("Valid annotation columns:\n")
  for (i in seq_along(valid_cols)) cat(sprintf("  %2d: %s\n", i, valid_cols[i]))
  flush.console()
  picked <- wb_smart_readline(
    "Select column(s) to use as groups -- comma-separated indexes or ranges (e.g. 1,3:5) (or 'quit' for all valid columns): ",
    valid = function(ch) {
      tokens <- wb_trim(strsplit(ch, ",")[[1]])
      if (!all(grepl("^[0-9]+(:[0-9]+)?$", tokens))) {
        return("Use indexes or ranges only (e.g. 1,3:5), try again.")
      }
      idx <- wb_parse_index_ranges(tokens)
      if (any(idx < 1 | idx > length(valid_cols))) {
        return(sprintf("Index out of range (1-%d), try again.", length(valid_cols)))
      }
      TRUE
    }
  )
  if (is.null(picked)) return(valid_cols)
  valid_cols[wb_parse_index_ranges(wb_trim(strsplit(picked, ",")[[1]]))]
}

wb_select_groups <- function(state, columns = NULL, max_categories = 10) {
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  all_cols <- colnames(annot)
  valid_cols <- setdiff(all_cols, c(REQUIRED_COLS, IGNORE_COLS))

  if (is.null(columns)) {
    if (!is.null(state$typemap$groups) &&
        wb_confirm("A groups file was provided. Use its columns as groups?")) {
      # header = FALSE -- this file is documented as a plain list, one annotation-column name
      # per line, with no header row. The default header = TRUE would otherwise silently
      # consume the first name as a column header and drop it from the result.
      columns <- as.vector(unlist(read.csv(state$typemap$groups, header = FALSE, stringsAsFactors = FALSE, quote = '"')))
      wb_msg("INFO", "Using columns from the provided groups file.")
    } else if (!is.null(state$typemap$groups)) {
      # Declined the groups file above -- let the user specify columns manually instead of
      # falling all the way back to "every valid column" (that's still available by quitting).
      wb_msg("INFO", "Discarding the provided groups file for this selection.")
      columns <- wb_prompt_indexed_columns(valid_cols)
    } else {
      columns <- valid_cols
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

wb_html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

wb_print_colors_plain <- function(state) {
  for (group in names(state$groups_colors)) {
    cat(sprintf("\n%s:\n", group))
    vals <- state$groups_colors[[group]]
    for (v in names(vals)) cat(sprintf("  %-20s %s\n", v, vals[[v]]))
  }
  flush.console()
}

# A bare hex code in cat() output is just text -- IRkernel's cell output is plain text, not a
# terminal that renders ANSI color codes as swatches, so it never shows an actual color preview
# no matter the Jupyter version. Renders `build_html()` via IRdisplay::display_html() instead
# when running live inside a Jupyter/IRkernel session, giving each value a colored box next to
# its hex code; falls back to `plain_fn()` outside a live kernel session (IRdisplay missing,
# running under plain Rscript, etc.), or if the rich-display call fails for any reason -- both
# the HTML build and the display call are covered by the same fallback.
wb_display_or_plain <- function(build_html, plain_fn) {
  in_jupyter <- requireNamespace("IRdisplay", quietly = TRUE) && isTRUE(getOption("jupyter.in_kernel"))
  if (in_jupyter) {
    ok <- tryCatch({ IRdisplay::display_html(build_html()); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(NULL))
  }
  plain_fn()
}

wb_print_colors <- function(state) {
  wb_display_or_plain(
    build_html = function() {
      html <- "<div style='font-family:monospace'>"
      for (group in names(state$groups_colors)) {
        html <- paste0(html, sprintf("<b>%s</b><br>", wb_html_escape(group)))
        vals <- state$groups_colors[[group]]
        for (v in names(vals)) {
          hex <- vals[[v]]
          html <- paste0(html, sprintf(
            "<span style='display:inline-block;width:14px;height:14px;background:%s;border:1px solid #888;margin-right:6px;vertical-align:middle'></span>%s: %s<br>",
            hex, wb_html_escape(v), hex
          ))
        }
      }
      paste0(html, "</div>")
    },
    plain_fn = function() wb_print_colors_plain(state)
  )
}

# Same swatch preview as wb_print_colors(), but for the indexed "value: hex" listing shown
# mid-edit in wb_edit_color() -- shows the current color next to each value so you can see
# what you're about to change, not just its hex code.
wb_print_indexed_values <- function(group, values, colors) {
  wb_display_or_plain(
    build_html = function() {
      html <- sprintf("<div style='font-family:monospace'><b>%s</b><br>", wb_html_escape(group))
      for (i in seq_along(values)) {
        hex <- colors[[i]]
        html <- paste0(html, sprintf(
          "<span style='display:inline-block;width:14px;height:14px;background:%s;border:1px solid #888;margin-right:6px;vertical-align:middle'></span>%2d: %-20s %s<br>",
          hex, i, wb_html_escape(values[i]), hex
        ))
      }
      paste0(html, "</div>")
    },
    plain_fn = function() {
      cat(sprintf("\n%s values:\n", group))
      for (i in seq_along(values)) cat(sprintf("  %2d: %-20s %s\n", i, values[i], colors[[i]]))
      flush.console()
    }
  )
}

wb_show_colors <- function(state) {
  if (length(state$groups_colors) == 0) {
    wb_msg("WARNING", "No colors assigned yet. Run wb_reset_colors() first.")
    wb_done()
    return(invisible(NULL))
  }
  wb_print_colors(state)
  wb_done()
  invisible(state$groups_colors)
}

wb_reset_colors <- function(state, annot = NULL) {
  if (is.null(annot)) annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  state$groups_colors <- wb_assign_colors(annot, state$groups_cols)
  wb_save_state(state)
}

# Accepts a hex color with or without a leading '#' (normalizing to include it) -- matches
# either mode below.
wb_normalize_hex <- function(hex) if (startsWith(hex, "#")) hex else paste0("#", hex)

wb_edit_color <- function(state, group = NULL, value = NULL, hex_color = NULL) {
  if (!is.null(group) && !is.null(value) && !is.null(hex_color)) {
    # Direct, non-interactive single edit -- e.g. for scripting.
    if (!grepl("^#[A-Fa-f0-9]{6}$", hex_color)) stop("hex_color must look like '#RRGGBB'.")
    if (!(group %in% names(state$groups_colors))) stop(sprintf("'%s' is not a known group.", group))
    key <- as.character(value)
    if (!(key %in% names(state$groups_colors[[group]]))) {
      stop(sprintf("'%s' is not a known value for group '%s'.", value, group))
    }
    state$groups_colors[[group]][[key]] <- hex_color
    return(wb_save_state(state))
  }

  # Interactive session -- ported from panda-src/build-config.r's change_current_colors():
  # pick a group, then either overwrite one value's color at a time, or replace every color
  # for that group at once from a comma-separated hex list. Loops across as many groups (and,
  # within a group, as many values) as wanted, until 'quit'.
  if (length(state$groups_colors) == 0) {
    wb_msg("WARNING", "No colors assigned yet. Run wb_reset_colors() first.")
    wb_done()
    return(invisible(state))
  }

  groups <- names(state$groups_colors)
  repeat {
    cat("\nGroups:\n")
    for (i in seq_along(groups)) cat(sprintf("  %2d: %s\n", i, groups[i]))
    flush.console()

    g_idx <- wb_smart_readline(
      "Enter group index to edit (or 'quit' to finish): ",
      valid = function(ch) {
        n <- suppressWarnings(as.integer(ch))
        if (is.na(n) || n < 1 || n > length(groups)) sprintf("Please enter a number from 1 to %d.", length(groups)) else TRUE
      },
      cancel_msg="Previous color-edits saved."
    )
    if (is.null(g_idx)) break
    this_group <- groups[as.integer(g_idx)]
    values <- names(state$groups_colors[[this_group]])

    wb_print_indexed_values(this_group, values, state$groups_colors[[this_group]])

    bulk <- wb_confirm(sprintf(
      "Set every color for '%s' at once (comma-separated hex list), instead of one at a time?", this_group),
      cancel_msg="Previous color-edits saved."
    )

    if (bulk) {
      nvals <- length(values)
      hexes <- wb_smart_readline(
        sprintf("Enter %d hex colors for '%s' (%s), comma-separated: ",
                nvals, this_group, paste(values, collapse = ", ")),
        valid = function(ch) {
          tokens <- wb_trim(strsplit(ch, ",")[[1]])
          if (length(tokens) != nvals) return(sprintf("Please enter exactly %d hex colors, comma-separated.", nvals))
          if (!all(grepl("^#?[A-Fa-f0-9]{6}$", tokens))) return("Each entry must look like '#RRGGBB', try again.")
          TRUE
        },
        cancel_msg="Previous color-edits saved."
      )
      if (!is.null(hexes)) {
        tokens <- vapply(wb_trim(strsplit(hexes, ",")[[1]]), wb_normalize_hex, character(1), USE.NAMES = FALSE)
        state$groups_colors[[this_group]] <- setNames(as.list(tokens), values)
        wb_msg("INFO", sprintf("Updated all %d colors for '%s'.", nvals, this_group))
      }
    } else {
      repeat {
        v_idx <- wb_smart_readline(
          sprintf("Enter value index within '%s' to edit (or 'quit' to stop editing this group): ", this_group),
          valid = function(ch) {
            n <- suppressWarnings(as.integer(ch))
            if (is.na(n) || n < 1 || n > length(values)) sprintf("Please enter a number from 1 to %d.", length(values)) else TRUE
          },
          cancel_msg="Previous color-edits saved."
        )
        if (is.null(v_idx)) break
        val_name <- values[as.integer(v_idx)]

        hex <- wb_smart_readline(
          sprintf("Enter hex color for '%s' (e.g. #RRGGBB): ", val_name),
          valid = function(ch) if (grepl("^#?[A-Fa-f0-9]{6}$", ch)) TRUE else "Invalid hex color, try again.",
          cancel_msg="Previous color-edits saved."
        )
        if (is.null(hex)) break
        state$groups_colors[[this_group]][[val_name]] <- wb_normalize_hex(hex)
        wb_msg("INFO", sprintf("Set '%s' -> '%s' to %s.", this_group, val_name, wb_normalize_hex(hex)))
      }
    }
  }

  wb_print_colors(state)
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

# Lightweight local approximation of what panoply_clumps_ptm_mapping (clumps_ptm_mapping.py)
# does at runtime: it splits the FASTA into per-sequence files (named per header, parsed using
# the FASTA_sep_type delimiter), then keeps only the ones whose name contains an accession_col
# value as a substring, hard-failing the whole Cromwell job if none match. That exact per-ID
# split isn't reproducible here (the parsing package isn't available in this repo), so this
# just checks for ANY substring overlap between accession values and raw FASTA header text --
# looser, but enough to flag a wholly mismatched ID scheme before a job gets submitted, without
# adding a FASTA-parsing dependency or reimplementing FASTA_sep_type's own delimiter logic.
wb_fasta_header_sample <- function(fasta_path, max_lines = 200000) {
  lines <- readLines(fasta_path, n = max_lines, warn = FALSE)
  sub("^>", "", grep("^>", lines, value = TRUE))
}

wb_accession_overlaps_fasta <- function(values, fasta_headers, max_values = 2000) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0 || length(fasta_headers) == 0) return(FALSE)
  values <- values[seq_len(min(length(values), max_values))]
  escaped <- values
  # Escape regex metacharacters one at a time (backslash first) rather than via a single
  # character-class pattern -- a class mixing '[', ']', and '{', '}' is awkward to get right
  # and portable across regex engines; this is simpler to verify correct.
  for (ch in c("\\", ".", "|", "(", ")", "[", "]", "{", "}", "^", "$", "*", "+", "?")) {
    escaped <- gsub(ch, paste0("\\", ch), escaped, fixed = TRUE)
  }
  pattern <- paste(escaped, collapse = "|")
  any(grepl(pattern, fasta_headers, perl = TRUE))
}

# Mirrors wb_validate_gene_id_column()'s shape, with one deliberate difference: if the
# configured accession column (default 'id.description', from
# panoply_ptm_normalization.accession_number_colname) is MISSING outright, this interactively
# fixes it (pick + validate + write an existing column in its place, like the gene-ID/
# metabolite-ID fixups elsewhere). But if the column EXISTS and just looks malformed or doesn't
# overlap the FASTA at all, that's only a warning -- unlike a missing column, there's no
# confident way to guess a better one, and Clumps-PTM's own mapping step remains the authority
# on whether it actually works.
wb_validate_clumpsptm_accession_column <- function(gct, gct_path, ome, accession_col, fasta_headers) {
  rdesc_names <- colnames(gct@rdesc)

  if (accession_col %in% rdesc_names) {
    vals <- gct@rdesc[[accession_col]]
    non_blank <- as.character(vals)[!is.na(vals) & nzchar(as.character(vals))]
    if (length(non_blank) == 0) {
      wb_msg("WARNING", sprintf(
        "Accession column '%s' in %s data is empty/all-NA -- Clumps-PTM's mapping step requires valid accession IDs here.",
        accession_col, toupper(ome)
      ))
    } else if (!wb_accession_overlaps_fasta(vals, fasta_headers)) {
      wb_msg("WARNING", sprintf(paste(
        "No overlap detected between '%s' values in %s data and identifiers in the provided",
        "FASTA -- Clumps-PTM's mapping step will likely fail. Double-check that the FASTA and",
        "this accession column use matching ID types (and FASTA_sep_type in",
        "master-parameters.yaml, if relevant)."
      ), accession_col, toupper(ome)))
    } else {
      wb_msg("INFO", sprintf("Accession column '%s' detected and appears to match the provided FASTA in %s data.", accession_col, toupper(ome)))
    }
    return(TRUE)
  }

  wb_msg("WARNING", sprintf("Accession column '%s' not found in %s data.", accession_col, toupper(ome)))
  cat(sprintf("\n%s row-annotation columns: %s\n\n", toupper(ome), paste(rdesc_names, collapse = ", ")))
  flush.console()
  col <- wb_smart_readline(
    sprintf("Column to use as '%s' for %s data (or 'quit' to skip): ", accession_col, toupper(ome)),
    valid = function(ch) {
      if (!(ch %in% rdesc_names)) return("Column not found, try again.")
      if (!wb_accession_overlaps_fasta(gct@rdesc[[ch]], fasta_headers)) {
        return(sprintf("'%s' does not appear to overlap the provided FASTA, try again.", ch))
      }
      TRUE
    }
  )
  if (is.null(col)) {
    wb_msg("WARNING", sprintf("Skipped accession-column setup for %s data. Clumps-PTM will not be run.", toupper(ome)))
    return(FALSE)
  }
  gct@rdesc[[accession_col]] <- gct@rdesc[[col]]
  wb_msg("INFO", sprintf("Using column '%s' as '%s' for %s data.", col, accession_col, toupper(ome)))
  wb_write_gct_atomic(gct, gct_path)
  TRUE
}

wb_select_clumpsptm_groups <- function(state, columns = NULL, fasta_path = NULL) {
  # panoply_clumps_ptm_workflow.wdl declares pSTY_gct/acK_gct/ubK_gct all optional ("must
  # include at least one") -- so any single PTM dataset is enough to offer Clumps-PTM.
  ptm_types <- c("phosphoproteome", "acetylome", "ubiquitylome")
  if (length(intersect(ptm_types, names(state$typemap))) < 1) {
    wb_msg("INFO", "No PTM data detected; Clumps-PTM will not be run.")
    state$toggles$run_clumpsptm <- FALSE
    return(wb_save_state(state))
  }

  state$toggles$run_clumpsptm <- wb_confirm("PTM data detected. Should Clumps-PTM be run?")
  if (!state$toggles$run_clumpsptm) return(wb_save_state(state))

  if (!is.null(fasta_path)) {
    # Explicit override -- validate and copy in like any other manually-specified path.
    if (!file.exists(fasta_path) || !grepl("\\.(fasta|fa)$", fasta_path, ignore.case = TRUE)) {
      stop(sprintf("'%s' is not an existing .fasta/.fa file.", fasta_path))
    }
    state$typemap$clumpsFASTA <- wb_copy_into_session(fasta_path)
  } else if (is.null(state$typemap$clumpsFASTA)) {
    # Not already mapped via wb_load_and_map_inputs() either -- prompt for it.
    path <- wb_smart_readline(
      "Path to the reference FASTA file (.fasta/.fa): ",
      valid = wb_validate_fasta_input
    )
    if (is.null(path)) {
      wb_msg("WARNING", "No reference FASTA provided. Clumps-PTM will not be run.")
      state$toggles$run_clumpsptm <- FALSE
      return(wb_save_state(state))
    }
    state$typemap$clumpsFASTA <- wb_copy_into_session(wb_resolve_fasta_path(path))
  }
  # else: state$typemap$clumpsFASTA was already mapped via wb_load_and_map_inputs() -- use it as-is.

  params <- if (!is.null(state$typemap$parameters)) yaml::read_yaml(state$typemap$parameters) else wb_load_default_master_parameters(state$github_ref)

  accession_col <- params$panoply_ptm_normalization$accession_number_colname %||% "id.description"
  fasta_headers <- wb_fasta_header_sample(state$typemap$clumpsFASTA)
  for (ptm_type in intersect(ptm_types, names(state$typemap))) {
    gct_path <- state$typemap[[ptm_type]]
    gct <- cmapR::parse_gctx(gct_path)
    if (!wb_validate_clumpsptm_accession_column(gct, gct_path, ptm_type, accession_col, fasta_headers)) {
      state$toggles$run_clumpsptm <- FALSE
      return(wb_save_state(state))
    }
  }

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

  max_categories <- params$panoply_clumps_ptm$diff_exp$max_annot_levels %||% 10
  columns <- wb_verify_clumpsptm_group_validity(annot, columns, max_categories)
  if (length(columns) == 0) {
    wb_msg("WARNING", "No annotation columns remain after filtering by subvalue count. Clumps-PTM will not be run.")
    state$toggles$run_clumpsptm <- FALSE
    return(wb_save_state(state))
  }
  # Not a hard cap -- just a strong nudge. Each additional annotation column (independent of
  # per-column subvalue count above) multiplies the number of Clumps-PTM sub-runs, so the user
  # is warned but free to proceed with as many as they like.
  if (length(columns) > 3) {
    wb_msg("WARNING", "More than 3 annotations selected; Clumps-PTM runs can become long and expensive.")
  }

  out_path <- file.path(dirname(state$typemap$groups %||% state$typemap$annotation), "groups-clumpsptm.csv")
  state$typemap$groups_clumpsptm <- wb_write_groups_file(annot, columns, out_path)
  wb_save_state(state)
}
