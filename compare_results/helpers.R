#
# compare_results/helpers.R
#
# Shared primitives for the unified-run comparator: filename normalization &
# pairing, generic file comparators (gct/csv/tsv/text/binary/parquet),
# parameter-file comparison, the module registry lookup, and the recursive
# tree walker. `compare.R` sources this file; `module_checks.R` (custom
# per-module comparators) relies on the generic comparators defined here.
#
# See rules.yaml for the config schema these functions implement:
#   strict:    exact/near-exact match (extends test-panoply.r's approach)
#   threshold: numeric-tolerance match (for stochastic/float-drift-prone files)
#   custom:    dispatch to a named function in module_checks.R, with args
#   ignore:    skip entirely, no report line
#   params:    compare a parameters file BEFORE results; annotates downstream
#              results if parameters differ meaningfully
#

suppressMessages(library(yaml))
suppressMessages(library(cmapR))

# Fallback defaults ONLY - for manual dev-testing with no YAML loaded. A real
# `compare.R` run overrides these via task_rules.yaml's `global:` block.
#
# DEFAULT_TOLERANCE:  GCTs are rounded to ndigits=5 pipeline-wide,
# so we can expect rounding error in the 5th decimal place.
# Most CSV/TSV data is GCT-derived and inherits the same rounding noise,
# even though R's write.csv/write.table itself doesn't  round.
# `strict:`-classified files resolve to this floor
DEFAULT_TOLERANCE <- 2.0e-05
DEFAULT_SIZE_DIFF_KB <- 5   # KB of absolute file-size difference considered "allowable" for binary_presence
DEFAULT_QUICK_SAMPLE_ROWS <- 2000
DEFAULT_QUICK_SAMPLE_COLS <- 200

# Whether a given module should run in "quick" (partial-sample) mode: an
# explicit --quick-only list takes precedence entirely; otherwise --quick
# applies to everything except modules named in --quick-except.
is_quick_module <- function(module_key, cfg) {
  g <- cfg$global
  if (!is.null(g$quick_include) && length(g$quick_include) > 0) return(module_key %in% g$quick_include)
  if (isTRUE(g$quick_default)) return(!(module_key %in% (g$quick_exclude %||% character(0))))
  FALSE
}

# all.equal() on a wide data.frame/matrix can emit one line per mismatched
# column/component - with thousands of columns that becomes megabytes of
# text per file. Cap it so the report stays readable; the important signal
# (how many components differ, the first few, the largest deltas) survives.
truncate_detail <- function(x, n = 600) {
  x <- paste(x, collapse = " | ")
  if (nchar(x) > n) paste0(substr(x, 1, n), sprintf(" ... [truncated, %d chars total]", nchar(x))) else x
}

# ---------------------------------------------------------------------------
# Result rows
# ---------------------------------------------------------------------------

# Status vocabulary:
#   PASS             matched (within tolerance where applicable)
#   FAIL             did not match, no legitimate explanation found
#   WARN             a custom comparator judged this "similar enough" but not
#                     identical (e.g. stochastic modules), or params differ
#   ABSENT           present in exactly one run
#   IGNORE           deliberately skipped per YAML `ignore` rule
#   INCOMPLETE_RULE  no rule (strict/threshold/custom/ignore) covers this file
STATUS_ORDER <- c("FAIL", "INCOMPLETE_RULE", "WARN", "ABSENT", "PASS", "IGNORE", "SKIPPED")

new_result <- function(module, relpath, status, check, detail = "", note = "") {
  data.frame(
    module = module, path = relpath, status = status, check = check,
    detail = detail, note = note, stringsAsFactors = FALSE
  )
}

empty_results <- function() {
  new_result(character(0), character(0), character(0), character(0), character(0), character(0))
}

# ---------------------------------------------------------------------------
# Filename normalization & pairing
# ---------------------------------------------------------------------------

apply_patterns <- function(x, patterns) {
  for (p in patterns) x <- gsub(p, "", x, perl = TRUE)
  x
}

# Strip a shared leading run/dataset prefix (e.g. "Tutorial-BRCA-all-" vs
# "Tutorial_v1_6-") from every name in a single directory listing, cutting
# only at a "-"/"_" boundary so tokens aren't split mid-word. No-op when
# fewer than 2 names are given (nothing to compare against).
strip_common_prefix <- function(x) {
  if (length(x) < 2) return(x)
  common <- Reduce(function(a, b) {
    n <- min(nchar(a), nchar(b))
    i <- 0
    while (i < n && substr(a, i + 1, i + 1) == substr(b, i + 1, i + 1)) i <- i + 1
    substr(a, 1, i)
  }, x)
  seps <- gregexpr("[-_/]", common)[[1]]
  cut <- if (seps[1] == -1) 0 else max(seps)
  if (cut > 0) substring(x, cut + 1) else x
}

# Pair up two vectors of names (files or directories) using progressively
# looser strategies: exact match, match after stripping configured "noise"
# regexes (job ids / cache digests / timestamps / uuids), match after
# stripping each side's own common dataset-prefix, and finally a positional
# fallback when the same (small) count of names remain on both sides. Always
# called on a single level of SIBLING names (see pair_tree below for how a
# whole module's recursive contents get paired level-by-level, which is what
# makes a flat common-prefix strip valid here - siblings share a directory,
# so nothing but the basename needs to vary).
# Returns list(pairs = list(list(a=,b=,note=)), only_a = chr, only_b = chr).
pair_names <- function(names1, names2, normalize_patterns = character(0), positional_fallback_max = 8,
                        prefix_strip_fn = strip_common_prefix) {
  pairs <- list()
  left <- unique(names1); right <- unique(names2)

  take_matches <- function(key_fn, note) {
    if (length(left) == 0 || length(right) == 0) return(invisible())
    # Key BOTH sides as a single batch call, not left one-at-a-time - a
    # per-item call loses the sibling context strip_common_prefix(_by_dir)
    # needs (it can only find a shared prefix across >= 2 names at once).
    left_keys <- key_fn(left)
    right_keys <- key_fn(right)
    matched_left <- character(0)
    for (i in seq_along(left)) {
      hit <- which(right_keys == left_keys[i])
      if (length(hit) >= 1) {
        hit <- hit[1]
        pairs[[length(pairs) + 1]] <<- list(a = left[i], b = right[hit], note = note)
        matched_left <- c(matched_left, left[i])
        right <<- right[-hit]
        right_keys <- right_keys[-hit]
      }
    }
    left <<- setdiff(left, matched_left)
  }

  take_matches(function(x) x, "exact")
  if (length(normalize_patterns) > 0) {
    take_matches(function(x) apply_patterns(x, normalize_patterns), "normalized")
  }
  take_matches(prefix_strip_fn, "prefix-stripped")
  if (length(normalize_patterns) > 0) {
    # Combined pass: some files carry BOTH a dataset prefix AND a
    # run-specific embedded token (e.g. a CACHE_DIGEST hash) - neither
    # strip alone pairs them, but stripping noise first and then the
    # dataset prefix does.
    take_matches(function(x) prefix_strip_fn(apply_patterns(x, normalize_patterns)), "prefix+normalized")
  }

  # Pairwise fallback for a dataset token embedded mid-filename rather than
  # as a clean leading prefix (e.g. "cna-analysis_Tutorial_v1_6-proteome.html"
  # vs "cna-analysis_Tutorial-BRCA-all-proteome.html" - shares a prefix AND a
  # suffix, just not a prefix alone). Greedily pairs whichever remaining
  # left/right names share the most combined prefix+suffix length, requiring
  # that shared span to cover at least half of the shorter name so two
  # genuinely different files don't get matched on a short coincidental
  # overlap. O(n*m) pairwise, but only ever runs on the leftovers after the
  # cheap stages above, which is small in practice.
  if (length(left) > 0 && length(right) > 0) {
    prefix_suffix_len <- function(a, b) {
      na <- nchar(a); nb <- nchar(b); n <- min(na, nb)
      i <- 0
      while (i < n && substr(a, i + 1, i + 1) == substr(b, i + 1, i + 1)) i <- i + 1
      j <- 0; max_j <- n - i
      while (j < max_j && substr(a, na - j, na - j) == substr(b, nb - j, nb - j)) j <- j + 1
      i + j
    }
    repeat {
      if (length(left) == 0 || length(right) == 0) break
      scores <- outer(left, right, Vectorize(prefix_suffix_len))
      best <- which(scores == max(scores), arr.ind = TRUE)[1, ]
      a <- left[best[1]]; b <- right[best[2]]
      if (scores[best[1], best[2]] < 0.5 * min(nchar(a), nchar(b))) break
      pairs[[length(pairs) + 1]] <- list(a = a, b = b, note = "common-prefix-suffix")
      left <- setdiff(left, a); right <- setdiff(right, b)
    }
  }

  if (length(left) > 0 && length(left) == length(right) && length(left) <= positional_fallback_max) {
    ord1 <- sort(left); ord2 <- sort(right)
    for (i in seq_along(ord1)) {
      pairs[[length(pairs) + 1]] <- list(a = ord1[i], b = ord2[i], note = "positional (names differ)")
    }
    left <- character(0); right <- character(0)
  }

  list(pairs = pairs, only_a = left, only_b = right)
}

# ---------------------------------------------------------------------------
# Glob matching (base R utils::glob2rx does the heavy lifting)
# ---------------------------------------------------------------------------

# A pattern is tested against the full relative path AND against the bare
# basename - so a simple pattern like "config.r" intuitively matches that
# file at any depth (mirroring how list.files(pattern=) behaves, which is
# basename-only even when recursive=TRUE), while a pattern that includes a
# "/" (e.g. "*/signature_gct/*") still gets genuine path-qualified matching.
glob_match_one <- function(relpath, pattern) {
  rx <- utils::glob2rx(pattern)
  grepl(rx, relpath) || grepl(rx, basename(relpath))
}

glob_hit <- function(relpath, patterns) {
  if (is.null(patterns) || length(patterns) == 0) return(FALSE)
  patterns <- unlist(patterns, use.names = FALSE)
  any(vapply(patterns, function(p) glob_match_one(relpath, p), logical(1)))
}

# Like glob_hit(), but returns the specific pattern that matched (first hit,
# NA if none) - used to report WHY a file was ignored, not just that it was.
glob_which <- function(relpath, patterns) {
  if (is.null(patterns) || length(patterns) == 0) return(NA_character_)
  patterns <- unlist(patterns, use.names = FALSE)
  for (p in patterns) if (glob_match_one(relpath, p)) return(p)
  NA_character_
}

custom_hit <- function(relpath, custom) {
  if (is.null(custom) || length(custom) == 0) return(NULL)
  for (pat in names(custom)) {
    if (glob_match_one(relpath, pat)) {
      spec <- custom[[pat]]
      if (is.character(spec)) return(list(fn = spec, args = list(), tolerance = NULL))
      return(list(fn = spec$fn, args = if (is.null(spec$args)) list() else spec$args, tolerance = spec$tolerance))
    }
  }
  NULL
}

# threshold entries are a list whose elements are EITHER a bare glob string
# (use the default tolerance) OR a single-key mapping {glob: tolerance}.
# Returns NA (matched, use default), a numeric tolerance (matched, explicit),
# or NULL (no match).
threshold_hit <- function(relpath, thr) {
  if (is.null(thr) || length(thr) == 0) return(NULL)
  for (item in thr) {
    if (is.character(item) && length(item) == 1) {
      if (glob_match_one(relpath, item)) return(NA)
    } else if (is.list(item) && length(item) >= 1 && !is.null(names(item))) {
      pat <- names(item)[1]
      if (nzchar(pat) && glob_match_one(relpath, pat)) return(item[[1]])
    }
  }
  NULL
}

# Classify a single relative file path against one module's rule set.
# Precedence: reproducible (only when `reproducible=TRUE`, i.e. the seed
# comparison in this module's `params:` block matched) > ignore > strict >
# custom > threshold > UNCOVERED. Strict is checked before threshold so an
# exact filename listed under `strict` wins over a broad wildcard the same
# file also happens to match under `threshold`.
classify_file <- function(relpath, rules, reproducible = FALSE) {
  if (reproducible && !is.null(rules$reproducible) && glob_hit(relpath, rules$reproducible$files)) {
    return(list(kind = rules$reproducible$mode %||% "threshold"))
  }
  if (glob_hit(relpath, rules$ignore)) return(list(kind = "ignore", pattern = glob_which(relpath, rules$ignore)))
  if (glob_hit(relpath, rules$strict)) return(list(kind = "strict"))
  spec <- custom_hit(relpath, rules$custom)
  if (!is.null(spec)) return(list(kind = "custom", spec = spec))
  thr <- threshold_hit(relpath, rules$threshold)
  if (!is.null(thr)) return(list(kind = "threshold", tolerance = thr))
  list(kind = "UNCOVERED")
}

# ---------------------------------------------------------------------------
# Generic file comparators
# ---------------------------------------------------------------------------

# cmapR's parser is chatty by default (one message per file: "parsing as GCT
# vX.X", then a rows/cols/descriptors summary) - not useful to a user running
# a comparison across hundreds of GCTs, so it's suppressed here rather than
# at every call site. Prefers parse_gctx (parse.gctx is deprecated and
# additionally emits a per-call deprecation warning).
read_gct <- function(path) {
  suppressMessages(
    if (exists("parse_gctx", where = asNamespace("cmapR"), inherits = FALSE)) {
      cmapR::parse_gctx(path)
    } else {
      cmapR::parse.gctx(path)
    }
  )
}

# Parameterized GCT comparator - the "some GCTs need matrix-only comparison,
# others need full rdesc/cdesc, others need a specific column subset, others
# need a looser tolerance for a specific computed column" knob panel lives
# here, driven by `args` from a `custom:` YAML entry. rdesc_tolerance/
# cdesc_tolerance each default to the matrix `tolerance` when unset.
# rdesc_cols/cdesc_cols are a KEEP-list (only these columns are compared at
# all); rdesc_ignore_cols/cdesc_ignore_cols are a DROP-list (these columns
# are excluded, everything else still compared) - applied AFTER the
# keep-list, so a column named in both is ignored (ignore always wins, same
# precedence convention as task_rules.yaml's own ignore: vs strict/
# threshold). Both act before rdesc_set_cols/rdesc_numeric_delim_cols get a
# chance to see the columns at all. Used directly (with no args) as the
# default comparator for *.gct under
# strict/threshold too.
#
# rdesc_set_cols/cdesc_set_cols: for rdesc/cdesc columns holding a
# delimiter-joined list (e.g. ssGSEA's "Signature.set.overlap.<sample>"
# gene lists) where element ORDER isn't stable across platforms even when
# the content is identical - compared per-cell as a SET (sorted,
# order-ignored) instead of exact string equality. Confirmed empirically on
# real same-config runs: the delimited text differs in ~84% of such cells,
# purely from ordering, while the set itself never differs. These are
# REGEXes (not globs, unlike the rest of this codebase), specifically so
# one can exclude a sibling column via a negative lookahead - e.g.
# "^Signature\\.set\\.overlap\\.(?!percent)" to match the gene-list columns
# but not their "...percent.<sample>" numeric siblings. Both share one
# `delimiter` (default "|") since nothing so far needs them to differ.
compare_delimited_set_cols <- function(d1, d2, pattern, delimiter) {
  cols <- if (!is.null(pattern)) intersect(grep(pattern, names(d1), value = TRUE, perl = TRUE), names(d2)) else character(0)
  notes <- character(0)
  for (col in cols) {
    v1 <- d1[[col]]; v2 <- d2[[col]]
    for (i in seq_len(min(length(v1), length(v2)))) {
      s1 <- v1[i]; s2 <- v2[i]
      if (identical(s1, s2)) next
      if (is.na(s1) || is.na(s2)) {
        notes <- c(notes, sprintf("%s[%s]: NA mismatch", col, rownames(d1)[i]))
        next
      }
      set1 <- sort(strsplit(s1, delimiter, fixed = TRUE)[[1]])
      set2 <- sort(strsplit(s2, delimiter, fixed = TRUE)[[1]])
      if (!identical(set1, set2)) {
        notes <- c(notes, sprintf("%s[%s]: only in run1: %s; only in run2: %s",
                                   col, rownames(d1)[i],
                                   paste(setdiff(set1, set2), collapse = ","),
                                   paste(setdiff(set2, set1), collapse = ",")))
      }
    }
  }
  list(notes = notes, cols = cols)
}

# General-purpose (not GCT-specific - usable for CSV/TSV data.frames too)
# per-cell NUMERIC-SEQUENCE comparison, order-preserving and tolerance-based,
# for delimiter-joined columns matched by a regex. Written for cmapR's
# all-or-nothing numeric coercion: fix_datatypes() reads a whole rdesc
# column back as character if even ONE cell fails to parse as a single
# plain number - which happens whenever panoply_association's
# duplicate.gene.policy concatenates multiple collapsed-duplicate-gene rows
# into one cell (e.g. AveExpr: "-1.98e-17|-1.15e-17"). That reverts every
# OTHER (perfectly ordinary, single-value) cell in the column to exact
# string comparison too - this recovers numeric-tolerance comparison for
# all of them uniformly, since a single-value cell is just a length-1
# sequence. Unlike compare_delimited_set_cols(), order is preserved (these
# are positionally meaningful duplicate-row values, not an arbitrary-order
# set) and values are compared with `tolerance`, not exact equality.
# Returns list(notes=<mismatch detail strings>, cols=<matched column names,
# to exclude from an ordinary all.equal() pass over the rest>).
compare_delimited_numeric_cols <- function(d1, d2, pattern, delimiter, tolerance = DEFAULT_TOLERANCE) {
  cols <- if (!is.null(pattern)) intersect(grep(pattern, names(d1), value = TRUE, perl = TRUE), names(d2)) else character(0)
  notes <- character(0)
  for (col in cols) {
    v1 <- d1[[col]]; v2 <- d2[[col]]
    for (i in seq_len(min(length(v1), length(v2)))) {
      s1 <- v1[i]; s2 <- v2[i]
      if (identical(s1, s2)) next
      if (is.na(s1) || is.na(s2)) {
        notes <- c(notes, sprintf("%s[%s]: NA mismatch", col, rownames(d1)[i]))
        next
      }
      n1 <- suppressWarnings(as.numeric(strsplit(as.character(s1), delimiter, fixed = TRUE)[[1]]))
      n2 <- suppressWarnings(as.numeric(strsplit(as.character(s2), delimiter, fixed = TRUE)[[1]]))
      if (length(n1) != length(n2) || anyNA(n1) || anyNA(n2)) {
        notes <- c(notes, sprintf("%s[%s]: not a comparable numeric sequence (run1=%s run2=%s)",
                                   col, rownames(d1)[i], s1, s2))
        next
      }
      check <- all.equal(n1, n2, tolerance = tolerance)
      if (!isTRUE(check)) {
        notes <- c(notes, sprintf("%s[%s]: %s", col, rownames(d1)[i], paste(check, collapse = "; ")))
      }
    }
  }
  list(notes = notes, cols = cols)
}

compare_gct <- function(f1, f2, tolerance = DEFAULT_TOLERANCE,
                         quick = FALSE, quick_sample_rows = DEFAULT_QUICK_SAMPLE_ROWS,
                         ignore_rdesc = FALSE, rdesc_cols = NULL, rdesc_ignore_cols = NULL, rdesc_tolerance = NULL,
                         rdesc_set_cols = NULL, rdesc_numeric_delim_cols = NULL,
                         ignore_cdesc = FALSE, cdesc_cols = NULL, cdesc_ignore_cols = NULL, cdesc_tolerance = NULL,
                         cdesc_set_cols = NULL, cdesc_numeric_delim_cols = NULL,
                         delimiter = "|", ...) {
  g1 <- read_gct(f1); g2 <- read_gct(f2)
  mat_tol <- tolerance
  notes <- character(0)

  # Positional comparison (matrix arithmetic, rdesc/cdesc row-by-row) only
  # makes sense if g1/g2 are actually in the same row/column order -
  # matching dimensions alone doesn't guarantee that, and nothing upstream
  # of this ever checked it. GCTs derived from topTable()-style output (e.g.
  # panoply_association's contrast GCTs) are typically sorted by
  # significance, which is sensitive to tiny platform-dependent noise in
  # the test statistic - two near-tied rows can swap rank between
  # platforms, silently reordering every GCT built from that result, which
  # would otherwise look like total nonsense (wrong genes compared against
  # each other) rather than the real, small difference it actually is.
  # Same ID SET in a different order -> silently realign run2 to run1's
  # order before comparing (a reordering is not a real difference). Actual
  # different ID sets (or duplicate-count differences, via sort()) -> FAIL
  # immediately and precisely, replacing the old dims-only check (which
  # couldn't tell "same count, different genes" from "same genes" at all).
  resolve_gct_alignment <- function(ids1, ids2, label) {
    if (identical(ids1, ids2)) return(list(kind = "ok"))
    if (!identical(sort(ids1), sort(ids2))) {
      return(list(kind = "mismatch", detail = sprintf(
        "%s IDs differ: %d only in run1, %d only in run2 (run1 n=%d, run2 n=%d)", label,
        length(setdiff(ids1, ids2)), length(setdiff(ids2, ids1)), length(ids1), length(ids2))))
    }
    list(kind = "reorder", order = match(ids1, ids2))
  }

  rows <- resolve_gct_alignment(g1@rid, g2@rid, "row")
  if (rows$kind == "mismatch") return(list(status = "FAIL", detail = rows$detail))
  if (rows$kind == "reorder") {
    g2@mat <- g2@mat[rows$order, , drop = FALSE]
    g2@rdesc <- g2@rdesc[rows$order, , drop = FALSE]
    g2@rid <- g2@rid[rows$order]
  }
  cols <- resolve_gct_alignment(g1@cid, g2@cid, "column")
  if (cols$kind == "mismatch") return(list(status = "FAIL", detail = cols$detail))
  if (cols$kind == "reorder") {
    g2@mat <- g2@mat[, cols$order, drop = FALSE]
    g2@cdesc <- g2@cdesc[cols$order, , drop = FALSE]
    g2@cid <- g2@cid[cols$order]
  }

  # Quick mode: cmapR has no partial-read shortcut used here, so this
  # doesn't save parse time, but it does cut the (often dominant) cost of
  # all.equal() over the full matrix/rdesc/cdesc - meaningful for the
  # larger cmap/NMF GCTs. Sample rows once so mat/rdesc stay aligned.
  quick_note <- ""
  if (isTRUE(quick)) {
    n <- quick_sample_rows
    orig_n <- nrow(g1@mat)
    if (orig_n > n) {
      idx <- sort(sample.int(orig_n, n))
      g1@mat <- g1@mat[idx, , drop = FALSE]; g2@mat <- g2@mat[idx, , drop = FALSE]
      g1@rdesc <- g1@rdesc[idx, , drop = FALSE]; g2@rdesc <- g2@rdesc[idx, , drop = FALSE]
      quick_note <- sprintf(" [quick mode: %d/%d rows sampled, not exhaustive]", n, orig_n)
    }
  }

  mat_check <- all.equal(g1@mat, g2@mat, tolerance = mat_tol)
  if (!isTRUE(mat_check)) notes <- c(notes, paste("matrix:", paste(mat_check, collapse = "; ")))

  if (!isTRUE(ignore_rdesc)) {
    r1 <- g1@rdesc; r2 <- g2@rdesc
    if (!is.null(rdesc_cols)) {
      r1 <- r1[, intersect(rdesc_cols, names(r1)), drop = FALSE]
      r2 <- r2[, intersect(rdesc_cols, names(r2)), drop = FALSE]
    }
    if (!is.null(rdesc_ignore_cols)) {
      r1 <- r1[, setdiff(names(r1), rdesc_ignore_cols), drop = FALSE]
      r2 <- r2[, setdiff(names(r2), rdesc_ignore_cols), drop = FALSE]
    }
    rset <- compare_delimited_set_cols(r1, r2, rdesc_set_cols, delimiter)
    if (length(rset$notes) > 0) notes <- c(notes, paste("rdesc set-columns:", paste(rset$notes, collapse = "; ")))
    rnum <- compare_delimited_numeric_cols(r1, r2, rdesc_numeric_delim_cols, delimiter, rdesc_tolerance %||% mat_tol)
    if (length(rnum$notes) > 0) notes <- c(notes, paste("rdesc numeric-delim columns:", paste(rnum$notes, collapse = "; ")))
    r1 <- r1[, setdiff(names(r1), c(rset$cols, rnum$cols)), drop = FALSE]
    r2 <- r2[, setdiff(names(r2), c(rset$cols, rnum$cols)), drop = FALSE]
    rcheck <- all.equal(r1, r2, tolerance = rdesc_tolerance %||% mat_tol, check.attributes = FALSE)
    if (!isTRUE(rcheck)) notes <- c(notes, paste("rdesc:", paste(rcheck, collapse = "; ")))
  }
  if (!isTRUE(ignore_cdesc)) {
    c1 <- g1@cdesc; c2 <- g2@cdesc
    if (!is.null(cdesc_cols)) {
      c1 <- c1[, intersect(cdesc_cols, names(c1)), drop = FALSE]
      c2 <- c2[, intersect(cdesc_cols, names(c2)), drop = FALSE]
    }
    if (!is.null(cdesc_ignore_cols)) {
      c1 <- c1[, setdiff(names(c1), cdesc_ignore_cols), drop = FALSE]
      c2 <- c2[, setdiff(names(c2), cdesc_ignore_cols), drop = FALSE]
    }
    cset <- compare_delimited_set_cols(c1, c2, cdesc_set_cols, delimiter)
    if (length(cset$notes) > 0) notes <- c(notes, paste("cdesc set-columns:", paste(cset$notes, collapse = "; ")))
    cnum <- compare_delimited_numeric_cols(c1, c2, cdesc_numeric_delim_cols, delimiter, cdesc_tolerance %||% mat_tol)
    if (length(cnum$notes) > 0) notes <- c(notes, paste("cdesc numeric-delim columns:", paste(cnum$notes, collapse = "; ")))
    c1 <- c1[, setdiff(names(c1), c(cset$cols, cnum$cols)), drop = FALSE]
    c2 <- c2[, setdiff(names(c2), c(cset$cols, cnum$cols)), drop = FALSE]
    ccheck <- all.equal(c1, c2, tolerance = cdesc_tolerance %||% mat_tol, check.attributes = FALSE)
    if (!isTRUE(ccheck)) notes <- c(notes, paste("cdesc:", paste(ccheck, collapse = "; ")))
  }

  if (length(notes) == 0) list(status = "PASS", detail = paste0("matched", quick_note))
  else list(status = "FAIL", detail = paste0(truncate_detail(notes), quick_note))
}

compare_delimited <- function(f1, f2, tolerance = DEFAULT_TOLERANCE, sep = ",",
                               quick = FALSE, quick_sample_rows = DEFAULT_QUICK_SAMPLE_ROWS,
                               quick_sample_cols = DEFAULT_QUICK_SAMPLE_COLS,
                               ignore_cols = NULL, align_by = NULL, ...) {
  # Quick mode: cap ROWS via nrows= (read.delim stops scanning the file
  # after N lines) and cap COLUMNS via colClasses="NULL" (read.table skips
  # parsing those fields entirely rather than reading-then-discarding them)
  # - a wide genome-wide correlation matrix (many columns, moderate rows)
  # needs the column cap to actually get faster; a tall one needs the row
  # cap. Column sampling is done by NAME so both sides read the same
  # logical columns even if their physical order/count differs.
  quick <- isTRUE(quick)
  nrows <- if (quick) quick_sample_rows else -1L
  max_cols <- if (quick) quick_sample_cols else Inf

  # read.delim(..., nrows=0) is NOT a cheap header peek in practice (R's
  # read.table internals do a full pre-scan regardless of nrows for some
  # defaults) - measured at ~1 minute on a 485MB file, versus <1 second via
  # readLines(n=1). Get the header this way instead whenever we might skip
  # reading the rest of the file anyway.
  peek_header <- function(f) strsplit(readLines(f, n = 1), sep, fixed = TRUE)[[1]]

  read_one <- function(f, header, keep_names = NULL) {
    base_args <- list(file = f, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE, row.names = NULL, nrows = nrows)
    # read.table/read.delim treat an explicitly-passed colClasses=NULL
    # differently from the argument being omitted entirely (the former
    # errors: "cannot replicate NULL to a non-zero length") - so only
    # include it in the call when there's an actual value to give it.
    if (!is.null(keep_names)) {
      col_classes <- rep("NULL", length(header))
      col_classes[header %in% keep_names] <- NA
      col_classes[1] <- NA # always keep the leading ID/rowname-ish column
      base_args$colClasses <- col_classes
    }
    do.call(utils::read.delim, base_args)
  }
  keep_names <- NULL; h1 <- NULL; h2 <- NULL; common <- character(0)
  if (quick) {
    h1 <- gsub('^"|"$', "", peek_header(f1)); h2 <- gsub('^"|"$', "", peek_header(f2))
    common <- intersect(h1, h2)
    if (length(common) > max_cols) keep_names <- sample(common, max_cols)
  }
  d1 <- tryCatch(read_one(f1, h1, keep_names), error = function(e) NULL)
  d2 <- tryCatch(read_one(f2, h2, keep_names), error = function(e) NULL)
  if (is.null(d1) || is.null(d2)) {
    return(list(status = "FAIL", detail = "could not parse as delimited text"))
  }
  # write.csv(row.names=TRUE)'s unnamed leading ID column round-trips through
  # read.delim(row.names=NULL) as a column literally named "" (or
  # "row.names", depending on separator) - an empty string is a valid
  # `names()` entry but can't actually be indexed via `[` or `[[`
  # ("undefined columns selected" / silently returns NULL), so rename it to
  # something addressable before any column-name-based logic below runs.
  if (!nzchar(names(d1)[1]) || names(d1)[1] == "row.names") names(d1)[1] <- "..rowid"
  if (!nzchar(names(d2)[1]) || names(d2)[1] == "row.names") names(d2)[1] <- "..rowid"
  if (!is.null(ignore_cols)) {
    d1 <- d1[, setdiff(names(d1), ignore_cols), drop = FALSE]
    d2 <- d2[, setdiff(names(d2), ignore_cols), drop = FALSE]
  }

  # Align rows by an ID column before comparing, when one is available -
  # protects against spurious wholesale mismatches caused by differing row
  # order (e.g. sample order) rather than differing content. Tries the
  # explicit `align_by` arg first, else the first column, and only applies
  # the reorder when it looks safe to do so (unique values, same set on
  # both sides) - e.g. a wide gene-symbol-as-column matrix with no row ID at
  # all correctly skips this (first column wouldn't have unique values).
  align_col <- align_by %||% names(d1)[1]
  # An empty-string column name (R's default header for an unnamed
  # row-index column written by write.csv(row.names=TRUE)) looks like a
  # normal name in `names()` but `d1[[""]]` returns NULL rather than the
  # column - nzchar() guards against silently "aligning" on nothing, which
  # would otherwise reduce match() to integer(0) and the whole frame to 0 rows.
  align_attempted <- !is.null(align_col) && nzchar(align_col) && nrow(d1) > 1 &&
    align_col %in% names(d1) && align_col %in% names(d2)
  align_ok <- align_attempted && !anyDuplicated(d1[[align_col]]) && !anyDuplicated(d2[[align_col]]) &&
    setequal(d1[[align_col]], d2[[align_col]])
  if (align_ok) {
    d2 <- d2[match(d1[[align_col]], d2[[align_col]]), , drop = FALSE]
    rownames(d2) <- NULL
  }
  # In quick mode, "first N rows" only samples the SAME logical rows on
  # both sides if row order is identical between runs; if it isn't, the ID
  # sets of the two truncated samples may barely overlap, and every
  # downstream diff below is really a sampling artifact, not a real
  # mismatch - flag that explicitly rather than let it read as a huge
  # genuine discrepancy.
  align_caveat <- if (quick && align_attempted && !align_ok) paste0(
    " [WARNING: row-ID sets don't fully overlap within the sampled rows - ",
    "quick mode assumes consistent row order between runs; this result may ",
    "be a sampling artifact, not a real diff - rerun without --quick to confirm]"
  ) else ""

  tol <- tolerance
  quick_note <- if (quick) paste0(sprintf(" [quick mode: %s rows%s sampled, not exhaustive]",
                                    if (nrows > 0) sprintf("first %d", nrows) else "all",
                                    if (!is.null(keep_names)) sprintf(", %d/%d cols", length(keep_names), length(common)) else ""),
                                   align_caveat) else ""

  # Report a column-set mismatch explicitly rather than letting all.equal()
  # produce a confusing "Length mismatch: comparison on first N components"
  # message - and still compare whatever columns ARE shared, so a partial
  # column-set difference (whatever its cause) doesn't hide a real value
  # mismatch in the rest.
  common_cols <- intersect(names(d1), names(d2))
  only1 <- setdiff(names(d1), names(d2)); only2 <- setdiff(names(d2), names(d1))
  if (length(only1) > 0 || length(only2) > 0) {
    msg <- sprintf("column set differs: %d common, %d only-in-run1, %d only-in-run2",
                    length(common_cols), length(only1), length(only2))
    if (length(common_cols) == 0) return(list(status = "FAIL", detail = paste0(msg, quick_note)))
    check <- all.equal(d1[, common_cols, drop = FALSE], d2[, common_cols, drop = FALSE],
                        tolerance = tol, check.attributes = FALSE)
    detail <- if (isTRUE(check)) paste0(msg, "; common columns match", quick_note)
              else paste0(msg, "; common columns also differ: ", truncate_detail(check), quick_note)
    return(list(status = "FAIL", detail = detail))
  }

  check <- all.equal(d1, d2, tolerance = tol, check.attributes = FALSE)
  if (isTRUE(check)) list(status = "PASS", detail = paste0("matched", quick_note))
  else list(status = "FAIL", detail = paste0(truncate_detail(check), quick_note))
}

compare_csv <- function(f1, f2, ...) compare_delimited(f1, f2, ..., sep = ",")
compare_tsv <- function(f1, f2, ...) compare_delimited(f1, f2, ..., sep = "\t")

compare_binary_presence <- function(f1, f2, tolerance = DEFAULT_SIZE_DIFF_KB, ...) {
  # Absolute byte-diff, not relative: PDF/PNG size drift here is almost
  # always a few KB of incidental noise (embedded job_id/timestamp text),
  # which is roughly CONSTANT regardless of file size - a relative/percentage
  # bound can't be tuned correctly across a mix of small and large files (see
  # task_rules.yaml's global: comment for the measured distribution behind
  # this). `tolerance` is in KB here, unlike every other comparator's
  # numeric-epsilon `tolerance`.
  s1 <- file.size(f1); s2 <- file.size(f2)
  absdiff_kb <- abs(s1 - s2) / 1024
  ok <- absdiff_kb <= tolerance
  list(status = if (ok) "PASS" else "FAIL",
       detail = sprintf("size run1=%d run2=%d (diff=%.3f KB, allow=%.3f KB)", s1, s2, absdiff_kb, tolerance))
}

compare_binary_md5 <- function(f1, f2, ...) {
  ok <- unname(tools::md5sum(f1)) == unname(tools::md5sum(f2)) # md5sum returns a named vector keyed by path; compare values not names
  list(status = if (ok) "PASS" else "FAIL", detail = if (ok) "byte-identical" else "byte content differs")
}

compare_text_exact <- function(f1, f2, normalize_patterns = character(0), ...) {
  read_txt <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  t1 <- apply_patterns(read_txt(f1), normalize_patterns)
  t2 <- apply_patterns(read_txt(f2), normalize_patterns)
  ok <- identical(t1, t2)
  list(status = if (ok) "PASS" else "FAIL",
       detail = if (ok) "matched (after normalization)" else "text differs after normalization")
}

compare_structured <- function(f1, f2, tolerance = DEFAULT_TOLERANCE, format = "yaml",
                                only_keys = NULL, ignore_keys = NULL, ...) {
  load_one <- function(f) {
    if (format == "yaml") yaml::read_yaml(f)
    else jsonlite::fromJSON(f, simplifyVector = TRUE)
  }
  d1 <- tryCatch(load_one(f1), error = function(e) NULL)
  d2 <- tryCatch(load_one(f2), error = function(e) NULL)
  if (is.null(d1) || is.null(d2)) return(list(status = "FAIL", detail = sprintf("could not parse as %s", format)))
  quote_list <- function(x) paste(sprintf("'%s'", x), collapse = ", ")

  # only_keys: restrict comparison to specific top-level key(s) - the
  # complement of ignore_keys - for a file whose OTHER keys are known to
  # differ for reasons unrelated to this task (e.g. a whole-pipeline config
  # snapshot duplicated into every module folder, where only this module's
  # own sub-block is worth checking). Applied before ignore_keys so both can
  # be combined if ever needed.
  mode_note <- "full file"
  if (!is.null(only_keys)) {
    missing1 <- setdiff(only_keys, names(d1)); missing2 <- setdiff(only_keys, names(d2))
    if (length(missing1) || length(missing2)) {
      return(list(status = "FAIL", detail = sprintf(
        "only_keys not found: run1 missing [%s], run2 missing [%s]",
        paste(missing1, collapse = ","), paste(missing2, collapse = ","))))
    }
    d1 <- d1[only_keys]
    d2 <- d2[only_keys]
    mode_note <- sprintf("only_keys: %s", quote_list(only_keys))
  }
  # A key named in ignore_keys but absent from both files is a no-op (nothing
  # to drop) - usually a stale rule or a typo, so it's surfaced as a WARN
  # rather than silently doing nothing.
  stale_ignore <- character(0)
  if (!is.null(ignore_keys)) {
    stale_ignore <- setdiff(ignore_keys, union(names(d1), names(d2)))
    d1 <- d1[setdiff(names(d1), ignore_keys)]
    d2 <- d2[setdiff(names(d2), ignore_keys)]
    mode_note <- if (identical(mode_note, "full file")) sprintf("ignore_keys: %s", quote_list(ignore_keys))
                 else paste0(mode_note, "; ignore_keys: ", quote_list(ignore_keys))
  }
  n_keys <- length(union(names(d1), names(d2)))
  compared_note <- sprintf("%s (%d top-level key%s compared)", mode_note, n_keys, if (n_keys == 1) "" else "s")

  tol <- tolerance
  check <- all.equal(d1, d2, tolerance = tol, check.attributes = FALSE)
  ok <- isTRUE(check)
  detail <- if (ok) sprintf("matched - %s", compared_note) else sprintf("%s - %s", truncate_detail(check), compared_note)
  status <- if (ok) "PASS" else "FAIL"
  if (length(stale_ignore)) {
    detail <- sprintf("%s [WARN: ignore_keys not found in either file: %s]", detail, quote_list(stale_ignore))
    if (ok) status <- "WARN"
  }
  list(status = status, detail = detail)
}

# Generic comparator for serialized-R-object files (*.rda/*.rds). Byte-diff
# is unreliable here (R serialization embeds version-tagged headers even for
# identical content), so this loads each side and compares the object(s)
# themselves with all.equal().
compare_rdata_object <- function(f1, f2, tolerance = DEFAULT_TOLERANCE,
                                  ignore_keys = character(0), normalize_patterns = character(0), ...) {
  # .rds holds a single serialized object and must be read with readRDS()
  # (load() errors on it: "bad restore file magic number" - it's a
  # different container format from .rda/.RData's multi-object save()).
  load_env <- function(f) {
    e <- new.env()
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
      assign("object", readRDS(f), envir = e)
    } else {
      load(f, envir = e)
    }
    e
  }
  e1 <- load_env(f1); e2 <- load_env(f2)
  n1 <- ls(e1); n2 <- ls(e2)
  common <- intersect(n1, n2); only1 <- setdiff(n1, n2); only2 <- setdiff(n2, n1)
  tol <- tolerance
  strip_strings <- function(v) {
    if (is.character(v) && length(normalize_patterns) > 0) return(apply_patterns(v, normalize_patterns))
    v
  }
  mismatches <- character(0)
  mismatch_detail <- character(0)
  for (nm in common) {
    v1 <- get(nm, envir = e1); v2 <- get(nm, envir = e2)
    if (is.list(v1) && is.list(v2) && !is.null(names(v1))) {
      v1 <- v1[setdiff(names(v1), ignore_keys)]; v2 <- v2[setdiff(names(v2), ignore_keys)]
      v1 <- lapply(v1, strip_strings); v2 <- lapply(v2, strip_strings)
      subkeys <- union(names(v1), names(v2))
      bad <- Filter(function(k) !isTRUE(all.equal(v1[[k]], v2[[k]], tolerance = tol, check.attributes = FALSE)), subkeys)
      if (length(bad)) {
        mismatches <- c(mismatches, nm)
        mismatch_detail <- c(mismatch_detail, sprintf("%s$%s: %s -> %s", nm, bad,
          sapply(bad, function(k) paste(v1[[k]], collapse = ",")), sapply(bad, function(k) paste(v2[[k]], collapse = ","))))
      }
      next
    }
    v1 <- strip_strings(v1); v2 <- strip_strings(v2)
    ok <- isTRUE(all.equal(v1, v2, tolerance = tol, check.attributes = FALSE))
    if (!ok) { mismatches <- c(mismatches, nm); mismatch_detail <- c(mismatch_detail, nm) }
  }
  detail <- sprintf("objects compared=%s", paste(common, collapse = ","))
  if (length(only1)) detail <- paste0(detail, "; run1-only=", paste(only1, collapse = ","))
  if (length(only2)) detail <- paste0(detail, "; run2-only=", paste(only2, collapse = ","))
  if (length(mismatches)) detail <- paste0(detail, "; MISMATCH: ", truncate_detail(mismatch_detail))
  list(status = if (length(mismatches) == 0) "PASS" else "FAIL", detail = detail)
}

compare_parquet <- function(f1, f2, tolerance = DEFAULT_TOLERANCE, ...) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    return(list(status = "WARN", detail = "arrow package not installed; skipped content check (presence only)"))
  }
  d1 <- arrow::read_parquet(f1); d2 <- arrow::read_parquet(f2)
  tol <- tolerance
  check <- all.equal(as.data.frame(d1), as.data.frame(d2), tolerance = tol, check.attributes = FALSE)
  if (isTRUE(check)) list(status = "PASS", detail = "matched")
  else list(status = "FAIL", detail = truncate_detail(check))
}

# Maps a file extension to the generic comparator used for `strict`/`threshold`
# classifications (i.e. everything NOT dispatched to a `custom` function).
default_dispatch <- function(ext) {
  ext <- tolower(ext)
  switch(ext,
    gct = "gct", gctx = "gct",
    csv = "csv",
    tsv = "tsv", txt_tab = "tsv",
    txt = "text", log = "text",
    json = "json",
    yaml = "yaml", yml = "yaml",
    pdf = "binary_presence", png = "binary_presence",
    jpg = "binary_presence", jpeg = "binary_presence",
    html = "binary_presence", htm = "binary_presence",
    xlsx = "binary_presence",
    parquet = "parquet",
    "binary_md5"
  )
}

run_generic <- function(disp, f1, f2, tolerance, size_tol, normalize_patterns,
                         quick = FALSE, quick_sample_rows = DEFAULT_QUICK_SAMPLE_ROWS) {
  switch(disp,
    gct = compare_gct(f1, f2, tolerance = tolerance, quick = quick, quick_sample_rows = quick_sample_rows),
    csv = compare_csv(f1, f2, tolerance = tolerance, quick = quick, quick_sample_rows = quick_sample_rows),
    tsv = compare_tsv(f1, f2, tolerance = tolerance, quick = quick, quick_sample_rows = quick_sample_rows),
    json = compare_structured(f1, f2, tolerance = tolerance, format = "json"),
    yaml = compare_structured(f1, f2, tolerance = tolerance, format = "yaml"),
    text = compare_text_exact(f1, f2, normalize_patterns = normalize_patterns),
    # binary_presence's "tolerance" is an absolute KB allowance, not the
    # file's numeric `tolerance` - size_tol (already resolved for
    # strict-vs-threshold by the caller) is what's threaded into that slot.
    binary_presence = compare_binary_presence(f1, f2, tolerance = size_tol),
    parquet = compare_parquet(f1, f2, tolerance = tolerance),
    compare_binary_md5(f1, f2)
  )
}

# ---------------------------------------------------------------------------
# Parameter-file comparison (params-first, per module)
# ---------------------------------------------------------------------------

load_params <- function(path, format, key = NULL) {
  switch(format,
    yaml = yaml::read_yaml(path),
    json = jsonlite::fromJSON(path, simplifyVector = TRUE),
    rdata = {
      e <- new.env()
      load(path, envir = e)
      obj_names <- ls(e)
      if (length(obj_names) == 0) stop("no objects found in Rdata file: ", path)
      nm <- if (!is.null(key) && key %in% obj_names) key else obj_names[1]
      as.list(get(nm, envir = e))
    },
    stop("unsupported params format: ", format)
  )
}

diff_params <- function(p1, p2, ignore_keys = character(0), normalize_patterns = character(0)) {
  keys <- union(names(p1), names(p2))
  keys <- setdiff(keys, ignore_keys)
  diffs <- list(); new_in_1 <- character(0); new_in_2 <- character(0)
  for (k in keys) {
    in1 <- k %in% names(p1); in2 <- k %in% names(p2)
    if (in1 && !in2) { new_in_1 <- c(new_in_1, k); next }
    if (!in1 && in2) { new_in_2 <- c(new_in_2, k); next }
    v1 <- p1[[k]]; v2 <- p2[[k]]
    # Strip env-leakage (absolute paths with run/workspace UUIDs, cache
    # digests, timestamps) from string-valued params before comparing, so a
    # param whose only "difference" is where the pipeline happened to run
    # doesn't get flagged as a real config change.
    if (is.character(v1) && length(normalize_patterns) > 0) v1 <- apply_patterns(v1, normalize_patterns)
    if (is.character(v2) && length(normalize_patterns) > 0) v2 <- apply_patterns(v2, normalize_patterns)
    eq <- isTRUE(all.equal(v1, v2))
    if (!eq) diffs[[k]] <- list(run1 = p1[[k]], run2 = p2[[k]])
  }
  list(diffs = diffs, new_in_1 = new_in_1, new_in_2 = new_in_2, differs = length(diffs) > 0)
}

# Locate + compare a module's declared parameters file, per its `params:`
# rule (file = glob, format = rdata|yaml|json, key = Rdata object name,
# ignore_keys = keys to drop before comparing). Returns a result row plus
# a caveat note to attach to every other result in this module when
# parameters differ (so a real config difference reads as "expected
# divergence", not a silent bug).
resolve_params <- function(dir1, dir2, params_rule, module, normalize_patterns = character(0)) {
  if (is.null(params_rule) || isTRUE(params_rule$skip)) {
    return(list(result = NULL, note = "", relpath1 = NA_character_, relpath2 = NA_character_, reproducible = FALSE))
  }
  find_one <- function(dir) {
    hits <- list.files(dir, recursive = TRUE, pattern = utils::glob2rx(params_rule$file), full.names = FALSE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }
  rel1 <- find_one(dir1); rel2 <- find_one(dir2)
  if (is.na(rel1) || is.na(rel2)) {
    res <- new_result(module, params_rule$file, "ABSENT", "params",
                       note = sprintf("params file missing in %s", if (is.na(rel1)) "run1" else "run2"))
    return(list(result = res, note = "", relpath1 = rel1, relpath2 = rel2, reproducible = FALSE))
  }
  f1 <- file.path(dir1, rel1); f2 <- file.path(dir2, rel2)
  p1 <- tryCatch(load_params(f1, params_rule$format, params_rule$key), error = function(e) NULL)
  p2 <- tryCatch(load_params(f2, params_rule$format, params_rule$key), error = function(e) NULL)
  if (is.null(p1) || is.null(p2)) {
    res <- new_result(module, params_rule$file, "FAIL", "params", detail = "could not parse params file")
    return(list(result = res, note = "", relpath1 = rel1, relpath2 = rel2, reproducible = FALSE))
  }

  # Seed-aware reproducibility: for stochastic modules, if the RNG seed (or
  # whatever key(s) actually govern reproducibility) happens to be equal
  # between the two runs, output that's normally only "similar-enough"
  # comparable becomes expected to match exactly (modulo env leakage) - see
  # `reproducible:` in a module's rule set for the file patterns this gates.
  seed_keys <- params_rule$seed_keys
  reproducible <- FALSE
  if (!is.null(seed_keys) && all(seed_keys %in% names(p1)) && all(seed_keys %in% names(p2))) {
    reproducible <- all(vapply(seed_keys, function(k) isTRUE(all.equal(p1[[k]], p2[[k]])), logical(1)))
  }

  d <- diff_params(p1, p2, params_rule$ignore_keys, normalize_patterns)
  if (d$differs) {
    detail <- paste(sprintf("%s: %s -> %s", names(d$diffs),
                             sapply(d$diffs, function(x) paste(x$run1, collapse = ",")),
                             sapply(d$diffs, function(x) paste(x$run2, collapse = ","))),
                     collapse = "; ")
    if (length(d$new_in_1) || length(d$new_in_2)) {
      detail <- paste0(detail, sprintf(" [new params: run1-only=%s run2-only=%s]",
                                        paste(d$new_in_1, collapse = ","), paste(d$new_in_2, collapse = ",")))
    }
    res <- new_result(module, params_rule$file, "WARN", "params", detail = detail)
    note <- sprintf("NOTE: module parameters differ (%s) - results below may reflect an intentional config change, not a defect",
                     paste(names(d$diffs), collapse = ", "))
    return(list(result = res, note = note, relpath1 = rel1, relpath2 = rel2, reproducible = reproducible))
  }
  res <- new_result(module, params_rule$file, "PASS", "params",
                     detail = if (reproducible) "matched (seed identical - exact reproducibility expected)" else "matched")
  list(result = res, note = "", relpath1 = rel1, relpath2 = rel2, reproducible = reproducible)
}

# Recursively pair a whole module's contents level-by-level: files at each
# directory level are paired among themselves (pair_names on plain
# siblings, where a flat prefix-strip is always valid), and subdirectories
# are paired the same way and then recursed into - which correctly handles
# a dataset prefix embedded in a SUBDIRECTORY name too (e.g.
# "Tutorial_v1_6-so_nmf-cna_NMF_postprocess/" vs
# "Tutorial-BRCA-all-so_nmf-cna_NMF_postprocess/"), not just in filenames.
# An entire unpaired subdirectory contributes its full recursive file list
# to only_a/only_b (so each such file is reported ABSENT individually,
# rather than collapsing the whole subtree into one vague line).
# Returns list(pairs = list(list(rel=, f1=, f2=)), only_a = chr, only_b = chr).
pair_tree <- function(dir1, dir2, normalize_patterns, prefix = "") {
  qualify <- function(names, p) if (nzchar(p)) file.path(p, names) else names

  c1 <- list_children(dir1); c2 <- list_children(dir2)
  pairs <- list()

  fp <- pair_names(c1$files, c2$files, normalize_patterns)
  for (p_ in fp$pairs) {
    pairs[[length(pairs) + 1]] <- list(rel = qualify(p_$a, prefix),
                                        f1 = file.path(dir1, p_$a), f2 = file.path(dir2, p_$b))
  }
  only_a <- qualify(fp$only_a, prefix)
  only_b <- qualify(fp$only_b, prefix)

  dp <- pair_names(c1$dirs, c2$dirs, normalize_patterns)
  for (p_ in dp$pairs) {
    sub <- pair_tree(file.path(dir1, p_$a), file.path(dir2, p_$b), normalize_patterns, qualify(p_$a, prefix))
    pairs <- c(pairs, sub$pairs)
    only_a <- c(only_a, sub$only_a)
    only_b <- c(only_b, sub$only_b)
  }
  for (d in dp$only_a) {
    rel <- qualify(d, prefix)
    leftover <- list.files(file.path(dir1, d), recursive = TRUE, full.names = FALSE)
    only_a <- c(only_a, if (length(leftover)) file.path(rel, leftover) else rel)
  }
  for (d in dp$only_b) {
    rel <- qualify(d, prefix)
    leftover <- list.files(file.path(dir2, d), recursive = TRUE, full.names = FALSE)
    only_b <- c(only_b, if (length(leftover)) file.path(rel, leftover) else rel)
  }
  list(pairs = pairs, only_a = only_a, only_b = only_b)
}

# ---------------------------------------------------------------------------
# Module dispatch + recursive tree walk
# ---------------------------------------------------------------------------

get_check_fn <- function(name) {
  get(name, mode = "function")
}

run_check <- function(cls, f1, f2, relpath, module_key, cfg, note, quick = FALSE) {
  if (cls$kind == "ignore") return(new_result(module_key, relpath, "IGNORE", "ignore", detail = sprintf("matched ignore pattern '%s'", cls$pattern)))
  if (isTRUE(getOption("compare_results.verbose"))) message(sprintf("  [%s] %s (%s)", module_key, relpath, cls$kind))

  if (cls$kind == "custom") {
    fn <- tryCatch(get_check_fn(cls$spec$fn), error = function(e) NULL)
    if (is.null(fn)) {
      return(new_result(module_key, relpath, "FAIL", cls$spec$fn,
                         detail = sprintf("custom comparator '%s' not found", cls$spec$fn), note = note))
    }
    # Ambient context every custom comparator MIGHT want (a rule's own args:
    # win on key collision); comparators that don't declare a matching
    # parameter simply absorb it via `...` and ignore it.
    default_ctx <- list(normalize_patterns = cfg$global$normalize_patterns,
                         quick = quick, quick_sample_rows = cfg$global$quick_sample_rows)
    call_args <- c(list(f1, f2), modifyList(default_ctx, cls$spec$args %||% list()))
    if (!is.null(cls$spec$tolerance)) call_args$tolerance <- cls$spec$tolerance
    out <- tryCatch(do.call(fn, call_args),
                     error = function(e) list(status = "FAIL", detail = paste("custom check error:", conditionMessage(e))))
    return(new_result(module_key, relpath, out$status, cls$spec$fn, out$detail, note))
  }

  if (cls$kind == "UNCOVERED") {
    return(new_result(module_key, relpath, "INCOMPLETE_RULE", "-",
                       detail = "no strict/threshold/custom/ignore rule matches this file - add one to rules.yaml",
                       note = note))
  }

  ext <- tools::file_ext(relpath)
  disp <- default_dispatch(ext)
  # Two independent resolutions - only one is ever actually consumed
  # (run_generic's switch routes `tol` to numeric comparators and `size_tol`
  # to binary_presence), but both need their own "strict" meaning since
  # they're different units/concepts:
  #   - tol: cfg$global$tolerance is a FLOOR, not an aspirational "exact
  #     match" value - GCTs are rounded to 5 decimal places on write, and
  #     CSV/TSV derived from GCT data inherits the same rounding noise (see
  #     task_rules.yaml's global: comment) - so `strict:` resolves to this
  #     floor directly, same as before.
  #   - size_tol: binary_presence (*.pdf/*.png/etc) is genuinely
  #     byte-size-exact under `strict:` (0 KB allowance) - unlike numeric
  #     comparators, there's no realistic-precision reason to allow ANY
  #     drift for a strict binary file; `threshold:` is what admits the
  #     incidental job_id/timestamp-driven noise these files actually have.
  # Either way, an explicit per-glob threshold: [{glob: tol}] override
  # always wins outright over its corresponding default.
  tol <- cfg$global$tolerance
  size_tol <- if (cls$kind == "strict") 0 else cfg$global$size_diff_kb
  if (cls$kind == "threshold" && !is.null(cls$tolerance) && !is.na(cls$tolerance)) {
    tol <- cls$tolerance
    size_tol <- cls$tolerance
  }

  out <- run_generic(disp, f1, f2, tol, size_tol, cfg$global$normalize_patterns,
                      quick = quick, quick_sample_rows = cfg$global$quick_sample_rows)
  new_result(module_key, relpath, out$status, disp, out$detail, note)
}

# Compare one already-identified module directory pair using its rule set:
# params-first, then recursive file enumeration/pairing/classification.
# Dispatch one task instance. `task_name` selects the ruleset from
# cfg$tasks (task_rules.yaml); `label` is what shows up in the report and
# defaults to the task name, but callers pass the originating folder name
# instead when that's more useful for telling instances apart (e.g. so_nmf's
# 4 per-omic panoply_nmf instances all use task_name="panoply_nmf" but get
# distinct labels like "so_nmf_results/<omic>").
dispatch_module <- function(dir1, dir2, task_name, cfg, label = task_name, reproducible_override = NULL) {
  # --only-tasks excludes a task WITHOUT removing its task_rules.yaml entry
  # (cfg$tasks stays intact - see compare.R) specifically so this can be
  # told apart from a genuine missing-rule config gap below. Checked first,
  # before the progress header, so an excluded module never gets a
  # "now comparing..." line for something that's about to be skipped.
  if (!is.null(cfg$only_tasks) && !(task_name %in% cfg$only_tasks)) {
    return(new_result(label, ".", "SKIPPED", "-",
                       detail = sprintf("excluded via --only-tasks (kept: %s)", paste(cfg$only_tasks, collapse = ", "))))
  }
  verbose <- isTRUE(getOption("compare_results.verbose"))
  # Always-on header (not gated by -v) - the one thing a user watching a long
  # run needs to know is "it's still moving, and here's where": which task
  # instance/folder is being compared right now. -v layers a much noisier
  # per-file trace on top of this (see the message() call in run_check()).
  message(sprintf("==> %s (%s)", task_name, label))
  rules <- cfg$tasks[[task_name]]
  if (is.null(rules)) {
    return(new_result(label, ".", "INCOMPLETE_RULE", "-",
                       detail = sprintf("no task_rules.yaml entry for task '%s'", task_name)))
  }

  # Escape hatch for tasks too structurally irregular for generic
  # file-by-file pairing (e.g. clumps-ptm's dataset+grouping-variable-named
  # result folder) - hands the whole (dir1, dir2) pair to a bespoke function.
  if (!is.null(rules$custom_module)) {
    fn <- get_check_fn(rules$custom_module)
    out <- tryCatch(fn(dir1, dir2, rules, cfg, label),
                     error = function(e) new_result(label, ".", "FAIL", rules$custom_module,
                                                      detail = paste("custom_module error:", conditionMessage(e))))
    return(out)
  }

  results <- list()
  note <- ""
  params_relpaths <- character(0)
  reproducible <- FALSE
  if (!is.null(rules$params)) {
    p <- resolve_params(dir1, dir2, rules$params, label, cfg$global$normalize_patterns)
    if (!is.null(p$result)) results[[length(results) + 1]] <- p$result
    note <- p$note
    params_relpaths <- c(p$relpath1, p$relpath2)
    reproducible <- isTRUE(p$reproducible)
  }
  # A group dispatch (see dispatch_group) resolves the seed ONCE from a
  # sibling task instance and hands the verdict down here, overriding
  # whatever (if anything) this task's own params: block computed - this is
  # what lets a task with no seed-bearing file of its own (e.g.
  # panoply_nmf_postprocess) still get seed-gated reproducibility.
  if (!is.null(reproducible_override)) reproducible <- isTRUE(reproducible_override)

  normp <- c(cfg$global$normalize_patterns, rules$normalize)
  tree <- pair_tree(dir1, dir2, normp)
  tree$pairs <- Filter(function(p_) !(p_$rel %in% params_relpaths), tree$pairs)
  quick <- is_quick_module(task_name, cfg)

  # Intermittent (time-throttled, not per-file) progress within a large
  # instance - e.g. a task with a handful of multi-GB GCT/CSV files can sit
  # quiet for minutes otherwise. Skipped entirely under -v, which already
  # prints every file via run_check()'s own message() - no need for both.
  n_files <- length(tree$pairs)
  last_progress <- Sys.time()
  for (i in seq_along(tree$pairs)) {
    p_ <- tree$pairs[[i]]
    cls <- classify_file(p_$rel, rules, reproducible)
    results[[length(results) + 1]] <- run_check(cls, p_$f1, p_$f2, p_$rel, label, cfg, note, quick)
    if (!verbose && n_files > 1 && difftime(Sys.time(), last_progress, units = "secs") >= 2) {
      message(sprintf("    ... %d/%d files (currently: %s)", i, n_files, p_$rel))
      last_progress <- Sys.time()
    }
  }
  for (only1 in setdiff(tree$only_a, params_relpaths)) {
    results[[length(results) + 1]] <- new_result(label, only1, "ABSENT", "-", note = "present in run1 only")
  }
  for (only2 in setdiff(tree$only_b, params_relpaths)) {
    results[[length(results) + 1]] <- new_result(label, only2, "ABSENT", "-", note = "present in run2 only")
  }

  if (length(results) == 0) return(empty_results())
  do.call(rbind, results)
}

# List immediate children of a directory, split into dirs/files (relative names only).
list_children <- function(dir) {
  entries <- list.files(dir, full.names = FALSE, all.files = FALSE, no.. = TRUE)
  entries <- entries[entries != ".DS_Store"]
  is_dir <- file.info(file.path(dir, entries))$isdir
  list(dirs = entries[is_dir %in% TRUE], files = entries[!(is_dir %in% TRUE)])
}

# Resolve a folder_map entry (or its absence) into a dispatch action for one
# already-paired (sub1, sub2, label) directory. Shared by walk_compare and
# dispatch_each_child so `uses`/`each_child`/`group`/pass-through behave
# identically regardless of how deep in the tree - or how many each_child
# levels - they were reached.
dispatch_by_folder_map <- function(sub1, sub2, key, cfg, label, on_unmapped) {
  entry <- cfg$folder_map[[key]]
  if (!is.null(entry$uses)) {
    dispatch_module(sub1, sub2, entry$uses, cfg, label)
  } else if (!is.null(entry$each_child)) {
    dispatch_each_child(sub1, sub2, entry$each_child, cfg, label)
  } else if (!is.null(entry$group)) {
    dispatch_group(sub1, sub2, entry, cfg, label)
  } else {
    on_unmapped(sub1, sub2, label)
  }
}

# Recursively walk two directory trees, matching directory names against
# folder_map.yaml at whichever depth they're found. A folder mapped with
# `uses:` stops recursion and dispatches its task; one mapped with
# `each_child:` recurses one level and dispatches each child separately
# (see dispatch_each_child); an unmapped folder is a pass-through and
# recurses normally.
walk_compare <- function(dir1, dir2, cfg, path_label = "") {
  results <- list()
  c1 <- list_children(dir1); c2 <- list_children(dir2)

  dir_pairs <- pair_names(c1$dirs, c2$dirs, cfg$global$normalize_patterns)
  for (p_ in dir_pairs$pairs) {
    key <- p_$a
    sub1 <- file.path(dir1, p_$a); sub2 <- file.path(dir2, p_$b)
    label <- if (nzchar(path_label)) paste0(path_label, "/", key) else key
    results[[length(results) + 1]] <- dispatch_by_folder_map(
      sub1, sub2, key, cfg, label,
      on_unmapped = function(s1, s2, lbl) walk_compare(s1, s2, cfg, lbl)
    )
  }
  for (only1 in dir_pairs$only_a) {
    label <- if (nzchar(path_label)) paste0(path_label, "/", only1) else only1
    results[[length(results) + 1]] <- new_result(label, ".", "ABSENT", "module",
                                                   note = "entire module/directory present in run1 only")
  }
  for (only2 in dir_pairs$only_b) {
    label <- if (nzchar(path_label)) paste0(path_label, "/", only2) else only2
    results[[length(results) + 1]] <- new_result(label, ".", "ABSENT", "module",
                                                   note = "entire module/directory present in run2 only")
  }

  # Stray files at a non-task level: compare generically with a caveat
  # note, unless explicitly stray-ignored (e.g. the whole-pipeline
  # master-config dumps that sit directly at an ome-wrapper root, outside
  # any task's own folder - already ignored inside every task for the same
  # reason; always-expected-to-differ, not worth a FAIL every run).
  file_pairs <- pair_names(c1$files, c2$files, cfg$global$normalize_patterns)
  for (p_ in file_pairs$pairs) {
    label <- if (nzchar(path_label)) path_label else "(root)"
    if (glob_hit(p_$a, cfg$global$stray_ignore)) {
      results[[length(results) + 1]] <- new_result(label, p_$a, "IGNORE", "ignore",
        detail = sprintf("matched stray_ignore pattern '%s'", glob_which(p_$a, cfg$global$stray_ignore)))
      next
    }
    ext <- tools::file_ext(p_$a)
    disp <- default_dispatch(ext)
    out <- run_generic(disp, file.path(dir1, p_$a), file.path(dir2, p_$b),
                        cfg$global$tolerance, cfg$global$size_diff_kb,
                        cfg$global$normalize_patterns)
    results[[length(results) + 1]] <- new_result(label, p_$a, out$status, disp, out$detail,
                                                   note = "file found outside any recognized task instance")
  }
  for (only1 in file_pairs$only_a) {
    label <- if (nzchar(path_label)) path_label else "(root)"
    results[[length(results) + 1]] <- new_result(label, only1, "ABSENT", "-", note = "present in run1 only")
  }
  for (only2 in file_pairs$only_b) {
    label <- if (nzchar(path_label)) path_label else "(root)"
    results[[length(results) + 1]] <- new_result(label, only2, "ABSENT", "-", note = "present in run2 only")
  }

  if (length(results) == 0) return(empty_results())
  do.call(rbind, results)
}

# A folder_map `each_child:` entry: don't dispatch this folder itself -
# treat each of its immediate children as its own separate instance of
# `default_task` (e.g. so_nmf_results/<omic>_NMF_results, one per omic).
# A child that has its own, more specific folder_map entry uses that
# instead, so each_child composes with nested each_child/uses.
dispatch_each_child <- function(dir1, dir2, default_task, cfg, path_label) {
  results <- list()
  c1 <- list_children(dir1); c2 <- list_children(dir2)
  dir_pairs <- pair_names(c1$dirs, c2$dirs, cfg$global$normalize_patterns)
  for (p_ in dir_pairs$pairs) {
    sub1 <- file.path(dir1, p_$a); sub2 <- file.path(dir2, p_$b)
    label <- paste0(path_label, "/", p_$a)
    results[[length(results) + 1]] <- dispatch_by_folder_map(
      sub1, sub2, p_$a, cfg, label,
      on_unmapped = function(s1, s2, lbl) dispatch_module(s1, s2, default_task, cfg, lbl)
    )
  }
  for (only1 in dir_pairs$only_a) {
    label <- paste0(path_label, "/", only1)
    results[[length(results) + 1]] <- new_result(label, ".", "ABSENT", "instance", note = "present in run1 only")
  }
  for (only2 in dir_pairs$only_b) {
    label <- paste0(path_label, "/", only2)
    results[[length(results) + 1]] <- new_result(label, ".", "ABSENT", "instance", note = "present in run2 only")
  }

  # Loose files directly inside an each_child folder (not part of any child
  # instance) fall through to the same generic stray comparator walk_compare
  # uses - unusual in practice (today's only each_child cases have none),
  # kept for robustness rather than assuming it can't happen.
  file_pairs <- pair_names(c1$files, c2$files, cfg$global$normalize_patterns)
  for (p_ in file_pairs$pairs) {
    if (glob_hit(p_$a, cfg$global$stray_ignore)) {
      results[[length(results) + 1]] <- new_result(path_label, p_$a, "IGNORE", "ignore",
        detail = sprintf("matched stray_ignore pattern '%s'", glob_which(p_$a, cfg$global$stray_ignore)))
      next
    }
    ext <- tools::file_ext(p_$a)
    disp <- default_dispatch(ext)
    out <- run_generic(disp, file.path(dir1, p_$a), file.path(dir2, p_$b),
                        cfg$global$tolerance, cfg$global$size_diff_kb, cfg$global$normalize_patterns)
    results[[length(results) + 1]] <- new_result(path_label, p_$a, out$status, disp, out$detail,
                                                   note = "file found outside any recognized task instance")
  }
  for (only1 in file_pairs$only_a) {
    results[[length(results) + 1]] <- new_result(path_label, only1, "ABSENT", "-", note = "present in run1 only")
  }
  for (only2 in file_pairs$only_b) {
    results[[length(results) + 1]] <- new_result(path_label, only2, "ABSENT", "-", note = "present in run2 only")
  }

  if (length(results) == 0) return(empty_results())
  do.call(rbind, results)
}

# Resolves one `group.children` family entry (an ordinary {uses: task} or
# {each_child: task}, exactly as used outside a group) into a plain list of
# not-yet-dispatched instances: list(label=, dir1=, dir2=, task=). Kept
# separate from actually dispatching so dispatch_group can resolve every
# family first, correlate them by position, and inject a shared
# reproducible verdict before any instance actually runs.
# `only_a`/`only_b` attributes carry ABSENT children through for each_child
# families (a `uses` family is a single instance and has none).
resolve_family_instances <- function(dir1, dir2, family_entry, cfg, path_label) {
  if (!is.null(family_entry$uses)) {
    return(list(list(label = path_label, dir1 = dir1, dir2 = dir2, task = family_entry$uses)))
  }
  if (!is.null(family_entry$each_child)) {
    c1 <- list_children(dir1); c2 <- list_children(dir2)
    dp <- pair_names(c1$dirs, c2$dirs, cfg$global$normalize_patterns)
    # Sort by run1's name so that, across every family in the group, the
    # i-th instance corresponds to the same logical instance (e.g. the same
    # omic) - see dispatch_group's header comment for why this is safe.
    ord <- order(vapply(dp$pairs, function(p_) p_$a, character(1)))
    dp$pairs <- dp$pairs[ord]
    instances <- lapply(dp$pairs, function(p_) {
      list(label = paste0(path_label, "/", p_$a),
           dir1 = file.path(dir1, p_$a), dir2 = file.path(dir2, p_$b),
           task = family_entry$each_child)
    })
    attr(instances, "only_a") <- dp$only_a
    attr(instances, "only_b") <- dp$only_b
    return(instances)
  }
  list() # a family entry must be uses/each_child; anything else resolves to nothing
}

# A folder_map `group:` entry: several sibling "family" folders (e.g.
# so_nmf_results/so_nmf_figures/so_nmf_ssgsea) that don't nest inside one
# another but represent the SAME set of task instances split by content
# type (raw results vs postprocess vs ssGSEA). Resolves each family's
# instances independently, then correlates them BY POSITION (i-th instance
# of every family = the same logical instance, e.g. the same omic) - safe
# here because every family's instance names share an identical prefix up
# to where they diverge (".../<omic>_NMF_results" vs
# ".../<omic>_NMF_postprocess"), so sorting alphabetically groups same-omic
# entries at the same index in every family without needing to parse the
# omic name out of the string at all.
#
# For each correlated instance, resolves `group.seed_params` ONCE against
# the `group.seed_from` family's member for that instance, then dispatches
# every family's member for that instance with the resulting reproducible
# verdict - this is what lets panoply_nmf_postprocess (which has no
# seed-bearing file of its own) get seed-gated reproducibility without
# reaching into a sibling task's folder itself.
dispatch_group <- function(dir1, dir2, group_spec, cfg, path_label) {
  results <- list()
  children <- group_spec$children
  family_names <- names(children)
  family_instances <- list()

  for (fam in family_names) {
    f1 <- file.path(dir1, fam); f2 <- file.path(dir2, fam)
    if (!dir.exists(f1) || !dir.exists(f2)) {
      results[[length(results) + 1]] <- new_result(
        paste0(path_label, "/", fam), ".", "ABSENT", "group-member",
        note = sprintf("family folder missing in %s", if (!dir.exists(f1)) "run1" else "run2"))
      next
    }
    inst <- resolve_family_instances(f1, f2, children[[fam]], cfg, paste0(path_label, "/", fam))
    family_instances[[fam]] <- inst
    for (only1 in attr(inst, "only_a")) {
      results[[length(results) + 1]] <- new_result(paste0(path_label, "/", fam, "/", only1), ".", "ABSENT",
                                                     "instance", note = "present in run1 only")
    }
    for (only2 in attr(inst, "only_b")) {
      results[[length(results) + 1]] <- new_result(paste0(path_label, "/", fam, "/", only2), ".", "ABSENT",
                                                     "instance", note = "present in run2 only")
    }
  }

  seed_fam <- group_spec$group$seed_from
  seed_instances <- family_instances[[seed_fam]]
  n <- length(seed_instances)

  for (fam in family_names) {
    if (fam == seed_fam || is.null(family_instances[[fam]])) next
    m <- length(family_instances[[fam]])
    if (m != n) {
      results[[length(results) + 1]] <- new_result(
        paste0(path_label, "/", fam), ".", "WARN", "group",
        detail = sprintf(
          "instance count differs from seed family '%s' (%d vs %d) - positional seed-correlation may be misaligned for this family",
          seed_fam, n, m))
    }
  }

  for (i in seq_len(n)) {
    seed_inst <- seed_instances[[i]]
    p <- resolve_params(seed_inst$dir1, seed_inst$dir2, group_spec$group$seed_params, seed_inst$label,
                         cfg$global$normalize_patterns)
    reproducible_i <- isTRUE(p$reproducible)
    for (fam in family_names) {
      insts <- family_instances[[fam]]
      if (is.null(insts) || i > length(insts)) next
      inst <- insts[[i]]
      results[[length(results) + 1]] <- dispatch_module(inst$dir1, inst$dir2, inst$task, cfg, inst$label,
                                                          reproducible_override = reproducible_i)
    }
  }

  # A group only explicitly owns its declared family folders - anything
  # else sitting alongside them (another subdirectory folder_map maps
  # separately, e.g. mo_nmf_ssgsea living next to mo_nmf's grouped
  # results/figures; or loose files, e.g. mo_nmf's own top-level .html
  # reports) needs the SAME pass-through handling walk_compare gives an
  # ordinary folder, or it would silently vanish instead of being reported.
  c1 <- list_children(dir1); c2 <- list_children(dir2)
  other_dirs <- function(x) setdiff(x, family_names)
  dir_pairs <- pair_names(other_dirs(c1$dirs), other_dirs(c2$dirs), cfg$global$normalize_patterns)
  for (p_ in dir_pairs$pairs) {
    sub1 <- file.path(dir1, p_$a); sub2 <- file.path(dir2, p_$b)
    label <- paste0(path_label, "/", p_$a)
    results[[length(results) + 1]] <- dispatch_by_folder_map(
      sub1, sub2, p_$a, cfg, label,
      on_unmapped = function(s1, s2, lbl) walk_compare(s1, s2, cfg, lbl)
    )
  }
  for (only1 in dir_pairs$only_a) {
    label <- paste0(path_label, "/", only1)
    results[[length(results) + 1]] <- new_result(label, ".", "ABSENT", "module", note = "present in run1 only")
  }
  for (only2 in dir_pairs$only_b) {
    label <- paste0(path_label, "/", only2)
    results[[length(results) + 1]] <- new_result(label, ".", "ABSENT", "module", note = "present in run2 only")
  }

  file_pairs <- pair_names(c1$files, c2$files, cfg$global$normalize_patterns)
  for (p_ in file_pairs$pairs) {
    if (glob_hit(p_$a, cfg$global$stray_ignore)) {
      results[[length(results) + 1]] <- new_result(path_label, p_$a, "IGNORE", "ignore",
        detail = sprintf("matched stray_ignore pattern '%s'", glob_which(p_$a, cfg$global$stray_ignore)))
      next
    }
    ext <- tools::file_ext(p_$a)
    disp <- default_dispatch(ext)
    out <- run_generic(disp, file.path(dir1, p_$a), file.path(dir2, p_$b),
                        cfg$global$tolerance, cfg$global$size_diff_kb, cfg$global$normalize_patterns)
    results[[length(results) + 1]] <- new_result(path_label, p_$a, out$status, disp, out$detail,
                                                   note = "file found outside any recognized task instance")
  }
  for (only1 in file_pairs$only_a) {
    results[[length(results) + 1]] <- new_result(path_label, only1, "ABSENT", "-", note = "present in run1 only")
  }
  for (only2 in file_pairs$only_b) {
    results[[length(results) + 1]] <- new_result(path_label, only2, "ABSENT", "-", note = "present in run2 only")
  }

  if (length(results) == 0) return(empty_results())
  do.call(rbind, results)
}

# Flattens task_rules.yaml's `workflows: {wf_name: {task_name: {...}}}`
# nesting (grouping is for human readability only) into one global
# task-name -> ruleset lookup. A task name repeated under two different
# workflow groups is a config error, not something to silently merge.
flatten_task_rules <- function(workflows) {
  tasks <- list()
  for (wf_name in names(workflows)) {
    wf <- workflows[[wf_name]]
    for (task_name in names(wf)) {
      if (task_name %in% names(tasks)) {
        stop(sprintf(
          "task_rules.yaml: task '%s' is declared under more than one workflow (also found under '%s') - task names must be globally unique",
          task_name, wf_name))
      }
      tasks[[task_name]] <- wf[[task_name]]
    }
  }
  tasks
}

# ---------------------------------------------------------------------------
# Entry point used by compare.R: handles being pointed at a results/ or
# reports/ root, a single category folder, or a single task's raw output
# directory directly - by starting the same recursive walk at whatever root
# is given. Resolution order for the root itself: an explicit --task
# override (bypasses folder_map.yaml entirely - the point of "a task is the
# base unit of comparison" is that this works even for output that never
# went through panoply_unified_assemble_results at all) > a folder_map.yaml
# match on the root's own basename > the basename being a task name
# directly > fall back to walking.
# ---------------------------------------------------------------------------

compare_roots <- function(root1, root2, cfg, task_override = NULL) {
  if (!is.null(task_override)) {
    return(dispatch_module(root1, root2, task_override, cfg, label = task_override))
  }
  key <- basename(normalizePath(root1, mustWork = FALSE))
  entry <- cfg$folder_map[[key]]
  if (!is.null(entry$uses)) {
    dispatch_module(root1, root2, entry$uses, cfg, label = key)
  } else if (!is.null(entry$each_child)) {
    dispatch_each_child(root1, root2, entry$each_child, cfg, path_label = key)
  } else if (!is.null(entry$group)) {
    dispatch_group(root1, root2, entry, cfg, path_label = key)
  } else if (key %in% names(cfg$tasks)) {
    dispatch_module(root1, root2, key, cfg, label = key)
  } else {
    walk_compare(root1, root2, cfg)
  }
}
