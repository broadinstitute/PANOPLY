#!/usr/bin/env Rscript
#
# compare_results/compare.R
#
# Compares a PANOPLY unified-run `results/` or `reports/` directory tree, a
# single task's raw output directory, or anything in between, between two
# runs. Comparison LOGIC lives in task_rules.yaml, keyed by task name;
# folder_map.yaml separately maps panoply_unified_assemble_results' on-disk
# folder names to those task names, so the same task's rules apply
# regardless of where the assemble step happens to put its output - see
# both files' header comments for the config schema, and helpers.R for the
# comparison engine.
#
# Usage:
#   Rscript compare.R path/to/run1 path/to/run2 [-o report.txt]
#   Rscript compare.R --run1 path/to/run1 --run2 path/to/run2 [-o report.txt]
#   Rscript compare.R <raw task dir> <raw task dir> --task panoply_nmf_postprocess
#
# run1/run2 may be given positionally (first two bare arguments) or via
# --run1/--run2 - not both: a named flag always wins over a positional for
# that slot, and any positional left over once both slots are filled is an
# error rather than a silent no-op.
#
if (!suppressMessages(requireNamespace("optparse", quietly = TRUE))) {
  stop("optparse package is required: install.packages('optparse')")
}
suppressMessages(library(optparse))

script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  "."
}
SCRIPT_DIR <- script_dir()
source(file.path(SCRIPT_DIR, "helpers.R"))
source(file.path(SCRIPT_DIR, "module_checks.R"))

option_list <- list(
  make_option(c("--run1"), dest = "run1", type = "character", default = NULL,
              help = "Path to the older/reference run's results/ or reports/ folder (or a single module folder). May also be given positionally (first bare argument) instead."),
  make_option(c("--run2"), dest = "run2", type = "character", default = NULL,
              help = "Path to the newer/test run's results/ or reports/ folder (same shape as run1). May also be given positionally (second bare argument) instead."),
  make_option(c("--task-rules"), dest = "task_rules_path", type = "character", default = NULL,
              help = sprintf("Path to task_rules.yaml [default: %s]", file.path(SCRIPT_DIR, "task_rules.yaml"))),
  make_option(c("--folder-map"), dest = "folder_map_path", type = "character", default = NULL,
              help = sprintf("Path to folder_map.yaml [default: %s]", file.path(SCRIPT_DIR, "folder_map.yaml"))),
  make_option(c("--task"), dest = "task", type = "character", default = NULL,
              help = "Force comparison as a single instance of this task, ignoring folder_map.yaml entirely - for comparing a task's raw output directly, independent of panoply_unified_assemble_results"),
  make_option(c("-o", "--output"), dest = "output", type = "character", default = "compare-report.txt",
              help = "Path to write the report [default: %default]"),
  make_option(c("--fail-on"), dest = "fail_on", type = "character", default = "FAIL",
              help = "Comma-separated statuses that cause a non-zero exit code [default: %default]"),
  make_option(c("--only-tasks"), dest = "only_tasks", type = "character", default = NULL,
              help = "Comma-separated task names to restrict comparison to (for fast iteration)"),
  make_option(c("--quick"), dest = "quick", action = "store_true", default = FALSE,
              help = "Sample a subset of rows from large tabular/GCT files instead of reading them in full - much faster on multi-GB files, at the cost of only spot-checking content"),
  make_option(c("--quick-sample-rows"), dest = "quick_sample_rows", type = "integer", default = NULL,
              help = sprintf("Rows sampled per file in quick mode [default: %d]", DEFAULT_QUICK_SAMPLE_ROWS)),
  make_option(c("--quick-except"), dest = "quick_except", type = "character", default = NULL,
              help = "Comma-separated task names to compare IN FULL even when --quick is set"),
  make_option(c("--quick-only"), dest = "quick_only", type = "character", default = NULL,
              help = "Comma-separated task names to run in quick mode; all others compared in full (overrides --quick/--quick-except)"),
  make_option(c("-v", "--verbose"), dest = "verbose", action = "store_true", default = FALSE,
              help = "Echo each task instance/file as it's checked"),
  make_option(c("--full-report"), dest = "full_report", action = "store_true", default = FALSE,
              help = "List every PASSed file (with its comparison method) in the report instead of collapsing PASS by module")
)
parsed <- parse_args(OptionParser(option_list = option_list, usage = "%prog [run1] [run2] [options]"),
                      positional_arguments = c(0, 2))
opt <- parsed$options
pos <- parsed$args
# Positional args only fill in whichever of run1/run2 wasn't given via flag -
# a named flag always wins, and anything left over once both slots are
# filled is an error rather than a silent no-op (e.g. --run1/--run2 given
# together with positional arguments, or one flag plus two positionals).
if (is.null(opt$run1) && length(pos) > 0) { opt$run1 <- pos[1]; pos <- pos[-1] }
if (is.null(opt$run2) && length(pos) > 0) { opt$run2 <- pos[1]; pos <- pos[-1] }
if (length(pos) > 0) {
  stop("Unexpected extra argument(s): ", paste(pos, collapse = ", "),
       " - run1/run2 are already set; pass them either positionally or via --run1/--run2, not both.")
}
options(compare_results.verbose = isTRUE(opt$verbose))

if (is.null(opt$run1) || is.null(opt$run2)) {
  stop("Both run1 and run2 are required (positionally or via --run1/--run2). Run with --help for usage.")
}
task_rules_path <- if (!is.null(opt$task_rules_path)) opt$task_rules_path else file.path(SCRIPT_DIR, "task_rules.yaml")
folder_map_path <- if (!is.null(opt$folder_map_path)) opt$folder_map_path else file.path(SCRIPT_DIR, "folder_map.yaml")
if (!file.exists(task_rules_path)) stop("task rules file not found: ", task_rules_path)
if (!file.exists(folder_map_path)) stop("folder map file not found: ", folder_map_path)
if (!dir.exists(opt$run1)) stop("run1 directory not found: ", opt$run1)
if (!dir.exists(opt$run2)) stop("run2 directory not found: ", opt$run2)

raw_task_rules <- yaml::read_yaml(task_rules_path)
raw_folder_map <- yaml::read_yaml(folder_map_path)
cfg <- list(
  global = modifyList(
    list(tolerance = DEFAULT_TOLERANCE, size_diff_kb = DEFAULT_SIZE_DIFF_KB, normalize_patterns = list()),
    raw_task_rules$global %||% list()
  ),
  tasks = flatten_task_rules(raw_task_rules$workflows %||% list()),
  folder_map = raw_folder_map$folder_map %||% list()
)
cfg$global$quick_default <- isTRUE(opt$quick)
cfg$global$quick_sample_rows <- if (!is.null(opt$quick_sample_rows)) opt$quick_sample_rows else DEFAULT_QUICK_SAMPLE_ROWS
if (!is.null(opt$quick_except)) cfg$global$quick_exclude <- trimws(strsplit(opt$quick_except, ",")[[1]])
if (!is.null(opt$quick_only)) cfg$global$quick_include <- trimws(strsplit(opt$quick_only, ",")[[1]])
if (!is.null(opt$only_tasks)) {
  # Recorded separately rather than filtering cfg$tasks itself, so a task
  # excluded here reports as SKIPPED (dispatch_module checks this first) -
  # not INCOMPLETE_RULE, which stays reserved for a genuine task_rules.yaml
  # gap on a task that WAS kept.
  cfg$only_tasks <- trimws(strsplit(opt$only_tasks, ",")[[1]])
}
if (!is.null(opt$task) && !(opt$task %in% names(cfg$tasks))) {
  stop(sprintf("--task '%s' has no entry in %s", opt$task, task_rules_path))
}

if (isTRUE(opt$verbose)) {
  message(sprintf("Loaded %d task rule(s) from %s, %d folder mapping(s) from %s",
                   length(cfg$tasks), task_rules_path, length(cfg$folder_map), folder_map_path))
}

results <- compare_roots(normalizePath(opt$run1), normalizePath(opt$run2), cfg, task_override = opt$task)
if (is.null(results) || nrow(results) == 0) {
  message("No files were compared - check that --run1/--run2 point at the right level of the tree.")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

fmt_rows <- function(df) {
  if (nrow(df) == 0) return("  (none)")
  paste(sprintf("  [%s] %s\n    check: %s\n    %s%s",
                df$module, df$path, df$check, df$detail,
                ifelse(nzchar(df$note), paste0("\n    ", df$note), "")),
        collapse = "\n")
}

fmt_collapsed <- function(df) {
  if (nrow(df) == 0) return("  (none)")
  agg <- as.data.frame(table(df$module))
  names(agg) <- c("module", "n")
  paste(sprintf("  %s: %d file(s)", agg$module, agg$n), collapse = "\n")
}

status_levels_present <- factor(results$status, levels = STATUS_ORDER)
tab <- table(results$module, status_levels_present)
summary_df <- as.data.frame.matrix(tab)
summary_df <- cbind(module = rownames(summary_df), summary_df)
totals <- c(module = "TOTAL", as.list(colSums(tab)))
# print.data.frame wraps onto stacked column-groups once a row's printed
# width exceeds getOption("width") (default 80) - a deeply-nested module
# path (e.g. so_nmf's per-omic instances) blows past that trivially, even
# though this is a flat text file with no real line-width constraint to
# respect. Widen it generously just for this one print call.
old_width <- options(width = 1000)
summary_lines <- c(
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"),
  "",
  sprintf("TOTAL: %s", paste(sprintf("%s=%d", names(colSums(tab)), colSums(tab)), collapse = ", "))
)
options(old_width)

fail_on <- trimws(strsplit(opt$fail_on, ",")[[1]])
n_fail_on <- sum(results$status %in% fail_on)

report_lines <- c(
  "PANOPLY Unified-Run Comparison Report",
  sprintf("run1 (older):  %s", opt$run1),
  sprintf("run2 (newer):  %s", opt$run2),
  sprintf("task rules:    %s", task_rules_path),
  sprintf("folder map:    %s", folder_map_path),
  "",
  "=== SUMMARY (by module) ===",
  summary_lines,
  "",
  sprintf("Exit status would be non-zero if any of [%s] present (see --fail-on). Count: %d",
          paste(fail_on, collapse = ","), n_fail_on),
  "",
  "=== FAIL ===", fmt_rows(results[results$status == "FAIL", , drop = FALSE]),
  "",
  "=== INCOMPLETE_RULE (no strict/threshold/custom/ignore rule matched this file) ===",
  fmt_rows(results[results$status == "INCOMPLETE_RULE", , drop = FALSE]),
  "",
  "=== WARN ===", fmt_rows(results[results$status == "WARN", , drop = FALSE]),
  "",
  "=== ABSENT (present in exactly one run) ===", fmt_rows(results[results$status == "ABSENT", , drop = FALSE]),
  "",
  if (isTRUE(opt$full_report)) {
    c("=== PASS (full file list) ===", fmt_rows(results[results$status == "PASS", , drop = FALSE]))
  } else {
    c("=== PASS (collapsed by module; pass --full-report for the full file list) ===",
      fmt_collapsed(results[results$status == "PASS", , drop = FALSE]))
  },
  "",
  if (isTRUE(opt$full_report)) {
    c("=== IGNORE (full file list) ===", fmt_rows(results[results$status == "IGNORE", , drop = FALSE]))
  } else {
    c("=== IGNORE (collapsed by module; pass --full-report for the full file list) ===",
      fmt_collapsed(results[results$status == "IGNORE", , drop = FALSE]))
  },
  "",
  "=== SKIPPED (excluded via --only-tasks) ===",
  fmt_rows(results[results$status == "SKIPPED", , drop = FALSE])
)

writeLines(report_lines, opt$output)
divider <- strrep("=", 60)
message("\n", divider)
message(sprintf("  Report written to: %s", opt$output))
message(divider, "\n")
message(paste(summary_lines, collapse = "\n"))

quit(status = if (n_fail_on > 0) 1 else 0)
