#
# compare_results/module_checks.R
#
# Custom per-module comparators referenced by `custom:`/`custom_module:`
# entries in rules.yaml. Add a new one here (matching the contract
# `function(f1, f2, ...)` - declare each knob a rule's `args:` might set as
# its own named parameter with a real default, plus a trailing `...` to
# absorb ambient dispatcher context it doesn't use - for a `custom:`
# file-pattern check, or `function(dir1, dir2, rules, cfg)` for a whole-
# module `custom_module:` escape hatch) and reference it by name in
# rules.yaml - no changes to helpers.R or compare.R are needed.
#
# Depends on helpers.R being sourced first (new_result, diff_params,
# compare_binary_presence, etc.)
#

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# NMF / consensus-clustering: label-permutation-invariant partition matching
# ---------------------------------------------------------------------------

# contingency: rows = labels in set A, cols = labels in set B (same K).
# Small K (NMF/consensus clustering rarely exceeds ~8), so brute-force over
# all column permutations is cheap and exact - no need for a linear-
# assignment package dependency.
best_label_permutation <- function(contingency) {
  k <- nrow(contingency)
  if (ncol(contingency) != k) stop("contingency table must be square")
  perms <- function(v) {
    if (length(v) <= 1) return(list(v))
    out <- list()
    for (i in seq_along(v)) for (r in perms(v[-i])) out[[length(out) + 1]] <- c(v[i], r)
    out
  }
  best_score <- -Inf; best_perm <- seq_len(k)
  for (perm in perms(seq_len(k))) {
    score <- sum(diag(contingency[, perm, drop = FALSE]))
    if (score > best_score) { best_score <- score; best_perm <- perm }
  }
  list(perm = best_perm, agree = best_score, total = sum(contingency), frac_agree = best_score / sum(contingency))
}

# Compares two sample -> cluster membership tables where cluster LABELS are
# an arbitrary permutation between runs (confirmed empirically: NMF/
# consensus-clustering factor numbering is not stable across reruns even
# when the underlying sample partition is), by finding the best label
# alignment and reporting the fraction of samples that agree under it.
# id_col defaults to the first column; label_col defaults to autodetecting
# a "cluster|membership|NMF"-matching column, else the second column.
compare_cluster_membership <- function(f1, f2, id_col = NULL, label_col = NULL,
                                        agree_threshold = 0.9, warn_threshold = 0.75, ...) {
  # row.names=NULL is essential here: these TSVs are written with an
  # unnamed leading sample-ID column (one fewer header field than data
  # fields), and read.delim's default behavior auto-detects that pattern
  # and silently promotes the ID into rownames() instead of a real column -
  # which then makes every column-name-based lookup below pick the wrong
  # column entirely.
  d1 <- utils::read.delim(f1, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  d2 <- utils::read.delim(f2, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  # read.delim auto-names that leading unnamed column literally "row.names",
  # which merge(by=) then treats as a reference to actual row names (not a
  # same-named column) and errors on - rename it to something unambiguous.
  if (names(d1)[1] == "row.names") names(d1)[1] <- "..id"
  if (names(d2)[1] == "row.names") names(d2)[1] <- "..id"
  id_col <- id_col %||% names(d1)[1]
  if (is.null(label_col)) {
    cand <- grep("cluster|membership|NMF", setdiff(names(d1), id_col), ignore.case = TRUE, value = TRUE)
    label_col <- if (length(cand)) cand[1] else setdiff(names(d1), id_col)[1]
  }
  if (!(id_col %in% names(d1)) || !(label_col %in% names(d1)) ||
      !(id_col %in% names(d2)) || !(label_col %in% names(d2))) {
    return(list(status = "FAIL", detail = sprintf(
      "expected columns not found (id_col=%s, label_col=%s)", id_col, label_col)))
  }
  m <- merge(d1[, c(id_col, label_col)], d2[, c(id_col, label_col)], by = id_col, suffixes = c(".1", ".2"))
  if (nrow(m) == 0) return(list(status = "FAIL", detail = "no overlapping sample IDs between runs"))

  l1 <- factor(m[[paste0(label_col, ".1")]]); l2 <- factor(m[[paste0(label_col, ".2")]])
  if (nlevels(l1) != nlevels(l2)) {
    return(list(status = "FAIL", detail = sprintf(
      "cluster count differs: run1 K=%d, run2 K=%d - not a label-permutation issue",
      nlevels(l1), nlevels(l2))))
  }
  tab <- table(l1, l2)
  bp <- best_label_permutation(tab)
  status <- if (bp$frac_agree >= agree_threshold) "PASS" else if (bp$frac_agree >= warn_threshold) "WARN" else "FAIL"
  list(status = status, detail = sprintf(
    "n=%d matched samples, best label-aligned agreement = %.1f%% (%d/%d) across %d cluster labels",
    nrow(m), 100 * bp$frac_agree, bp$agree, bp$total, nlevels(l1)))
}

# ---------------------------------------------------------------------------
# CLUMPS-PTM: whole-module custom_module escape hatch
# ---------------------------------------------------------------------------
# The result folder name embeds BOTH the dataset name and the grouping
# variable tested (e.g. "Tutorial-BRCA-all_PIK3CA.mutation_clumps_ptm_full_
# results"), with no stable shared substring across runs when the grouping
# variable itself differs - generic file-pairing can't handle this, so the
# whole clumpsptm/ category is handed here instead.
compare_clumpsptm_module <- function(dir1, dir2, rules, cfg, label = "clumpsptm") {
  find_result_dir <- function(root) {
    kids <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    hits <- kids[grepl("_clumps_ptm_full_results$", basename(kids))]
    if (length(hits) == 0) NA_character_ else hits[1]
  }
  r1 <- find_result_dir(dir1); r2 <- find_result_dir(dir2)
  results <- list()
  if (is.na(r1) || is.na(r2)) {
    results[[length(results) + 1]] <- new_result(
      label, ".", "ABSENT", "custom_module",
      note = sprintf("clumps-ptm result folder missing in %s", if (is.na(r1)) "run1" else "run2"))
    return(do.call(rbind, results))
  }

  # 1. params.yaml - compare the analysis-defining keys first
  find_params <- function(root) {
    hits <- list.files(root, pattern = "^params\\.yaml$", recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }
  pf1 <- find_params(r1); pf2 <- find_params(r2)
  grouping_matches <- NA
  if (!is.na(pf1) && !is.na(pf2)) {
    p1 <- yaml::read_yaml(pf1); p2 <- yaml::read_yaml(pf2)
    # output_prefix is "<dataset-name>_<grouping-variable>" - the dataset
    # name varies release-to-release for reasons that have nothing to do
    # with clumps-ptm's config, so compare only the grouping variable
    # itself (its last token; the variable names seen so far -
    # "PIK3CA.mutation", "Proteome.Cluster" - use "." not "_", so the last
    # "_"-split segment reliably isolates it regardless of how the dataset
    # name is punctuated).
    grouping_var <- function(x) if (is.null(x)) NA_character_ else utils::tail(strsplit(x, "_")[[1]], 1)
    key1 <- list(grouping_var = grouping_var(p1$postprocess$output_prefix), pymol_gen = p1$postprocess$pymol_gen,
                 fdr_threshold = p1$postprocess$fdr_threshold)
    key2 <- list(grouping_var = grouping_var(p2$postprocess$output_prefix), pymol_gen = p2$postprocess$pymol_gen,
                 fdr_threshold = p2$postprocess$fdr_threshold)
    d <- diff_params(key1, key2)
    status <- if (d$differs) "WARN" else "PASS"
    detail <- if (d$differs) {
      paste(sprintf("%s: %s vs %s", names(d$diffs),
                     sapply(d$diffs, function(x) paste(x$run1, collapse = ",")),
                     sapply(d$diffs, function(x) paste(x$run2, collapse = ","))), collapse = "; ")
    } else "matched"
    detail <- sprintf("%s [output_prefix: %s vs %s]", detail, p1$postprocess$output_prefix, p2$postprocess$output_prefix)
    results[[length(results) + 1]] <- new_result(label, "params.yaml", status, "params", detail)
    grouping_matches <- !("grouping_var" %in% names(d$diffs))
  } else {
    results[[length(results) + 1]] <- new_result(label, "params.yaml", "ABSENT", "params")
  }

  # 2. group folder discovery (folders named "<index>_positive_results" / "<index>_negative_results")
  group_dirs <- function(root) {
    kids <- basename(list.dirs(root, recursive = FALSE))
    kids[grepl("^[0-9]+_(positive|negative)_results$", kids)]
  }
  g1 <- group_dirs(r1); g2 <- group_dirs(r2)
  idx1 <- sort(unique(as.integer(sub("_.*$", "", g1))))
  idx2 <- sort(unique(as.integer(sub("_.*$", "", g2))))
  results[[length(results) + 1]] <- new_result(
    label, "groups", "PASS", "info",
    detail = sprintf("run1 groups=%s (n=%d); run2 groups=%s (n=%d)",
                      paste(idx1, collapse = ","), length(idx1), paste(idx2, collapse = ","), length(idx2)))

  comparable <- isTRUE(grouping_matches) && length(idx1) == length(idx2) && length(idx1) > 0

  if (!comparable) {
    results[[length(results) + 1]] <- new_result(
      label, "groups", "WARN", "comparability",
      detail = paste(
        "runs used a different grouping variable and/or a different number of groups",
        "(see params.yaml diff above) - group-level statistical comparison skipped as not",
        "meaningful; only output schema was checked"))
    f1 <- list.files(r1, pattern = "clumpsptm_output\\.tsv$", recursive = TRUE, full.names = TRUE)
    f2 <- list.files(r2, pattern = "clumpsptm_output\\.tsv$", recursive = TRUE, full.names = TRUE)
    if (length(f1) && length(f2)) {
      h1 <- names(utils::read.delim(f1[1], sep = "\t", nrows = 0, check.names = FALSE))
      h2 <- names(utils::read.delim(f2[1], sep = "\t", nrows = 0, check.names = FALSE))
      ok <- length(h1) == length(h2)
      results[[length(results) + 1]] <- new_result(
        label, "clumpsptm_output.tsv", if (ok) "PASS" else "FAIL", "schema",
        detail = sprintf("run1 cols(%d)=%s | run2 cols(%d)=%s",
                          length(h1), paste(h1, collapse = ","), length(h2), paste(h2, collapse = ",")))
    }
  } else {
    for (i in seq_along(idx1)) {
      for (posneg in c("positive", "negative")) {
        gd1 <- file.path(r1, sprintf("%d_%s_results", idx1[i], posneg))
        gd2 <- file.path(r2, sprintf("%d_%s_results", idx2[i], posneg))
        out1 <- list.files(gd1, pattern = "clumpsptm_output\\.tsv$", full.names = TRUE)
        out2 <- list.files(gd2, pattern = "clumpsptm_output\\.tsv$", full.names = TRUE)
        if (length(out1) == 0 || length(out2) == 0) next
        d1 <- utils::read.delim(out1[1], sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
        d2 <- utils::read.delim(out2[1], sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
        key_cols <- intersect(c("geneSymbol", "variableSites"), intersect(names(d1), names(d2)))
        if (length(key_cols) == 0) next
        m <- merge(d1, d2, by = key_cols, suffixes = c(".1", ".2"))
        if (nrow(m) < 3 || !all(c("clumpsptm_pval.1", "clumpsptm_pval.2") %in% names(m))) next
        rho <- suppressWarnings(stats::cor(m$clumpsptm_pval.1, m$clumpsptm_pval.2, method = "spearman", use = "complete.obs"))
        sig1 <- if ("fdr_pass.1" %in% names(m)) m[[key_cols[1]]][m$fdr_pass.1 %in% TRUE] else character(0)
        sig2 <- if ("fdr_pass.2" %in% names(m)) m[[key_cols[1]]][m$fdr_pass.2 %in% TRUE] else character(0)
        jac <- if (length(union(sig1, sig2)) == 0) 1 else length(intersect(sig1, sig2)) / length(union(sig1, sig2))
        status <- if (is.na(rho)) "WARN" else if (rho >= 0.75 && jac >= 0.5) "PASS" else if (rho >= 0.5) "WARN" else "FAIL"
        results[[length(results) + 1]] <- new_result(
          label, sprintf("%d_%s_results/clumpsptm_output.tsv", idx1[i], posneg), status, "pval_rank_cor",
          detail = sprintf("n=%d spearman_rho(pval)=%.3f jaccard(fdr_pass)=%.3f", nrow(m), rho, jac))
      }
    }
  }

  # 3. summary figures - presence/size only
  fig1 <- file.path(r1, "figures"); fig2 <- file.path(r2, "figures")
  for (fn in c("clumpsptm_summary_figure.pdf", "clumpsptm_summary_figure.png")) {
    p1 <- file.path(fig1, fn); p2 <- file.path(fig2, fn)
    if (file.exists(p1) && file.exists(p2)) {
      out <- compare_binary_presence(p1, p2)
      results[[length(results) + 1]] <- new_result(label, file.path("figures", fn), out$status, "binary_presence", out$detail)
    } else {
      results[[length(results) + 1]] <- new_result(label, file.path("figures", fn), "ABSENT", "-")
    }
  }

  # 4. all_html_reports - normalize the "file vs directory-containing-one-file" quirk before comparing
  find_report_html <- function(category_root) {
    p <- file.path(category_root, "all_html_reports")
    if (dir.exists(p)) {
      hits <- list.files(p, pattern = "\\.html$", full.names = TRUE)
      if (length(hits)) return(hits[1])
    } else if (file.exists(p)) {
      return(p)
    }
    NA_character_
  }
  h1 <- find_report_html(dir1); h2 <- find_report_html(dir2)
  if (!is.na(h1) && !is.na(h2)) {
    out <- compare_binary_presence(h1, h2, tolerance = 0.8)
    results[[length(results) + 1]] <- new_result(label, "all_html_reports", out$status, "binary_presence", out$detail)
  } else {
    results[[length(results) + 1]] <- new_result(label, "all_html_reports", "ABSENT", "-")
  }

  # 5. per-protein parquet detail - deliberately not deep-compared, count only
  n_parq <- function(root) length(list.files(root, pattern = "clumpsptm\\.parquet$", recursive = TRUE))
  results[[length(results) + 1]] <- new_result(
    label, "*/clumpsptm/*.parquet", "IGNORE", "ignore",
    detail = sprintf("run1 n=%d run2 n=%d (per-protein permutation detail, not deep-compared)", n_parq(r1), n_parq(r2)))

  do.call(rbind, results)
}
