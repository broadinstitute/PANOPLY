#!/usr/bin/env Rscript
#
# compare_results/validate_rules.R
#
# Dev-level sanity check, two parts:
#   1. task_rules.yaml <-> real WDL: cross-references each task's `wdl_task`
#      (defaulting to the task's own key when not given) against the actual
#      panoply_* task definitions/calls found under hydrant/, so the config
#      doesn't silently drift out of sync as tasks are added/renamed/removed.
#   2. folder_map.yaml <-> task_rules.yaml: every `uses`/`each_child` target
#      must have a real task_rules.yaml entry - this is a HARD check (a
#      typo here means comparisons silently report INCOMPLETE_RULE at
#      runtime instead of failing loudly at config-load time).
#
# Part 1 is best-effort regex over the WDL files, not authoritative - some
# flagged gaps are expected (e.g. purely intermediate/passthrough tasks with
# no results/ subfolder of their own, or synthetic assemble-step
# identifiers with `wdl_task: null`); use judgement on that part of the
# output. Part 2 is exact and should always be clean.
#
# Usage: Rscript validate_rules.R [path/to/task_rules.yaml] [path/to/folder_map.yaml]
#
suppressMessages(library(yaml))

script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  "."
}
SCRIPT_DIR <- script_dir()
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
source(file.path(SCRIPT_DIR, "helpers.R")) # for flatten_task_rules() and %||%

args <- commandArgs(trailingOnly = TRUE)
task_rules_path <- if (length(args) >= 1) args[1] else file.path(SCRIPT_DIR, "task_rules.yaml")
folder_map_path <- if (length(args) >= 2) args[2] else file.path(SCRIPT_DIR, "folder_map.yaml")
if (!file.exists(task_rules_path)) stop("task_rules.yaml not found: ", task_rules_path)
if (!file.exists(folder_map_path)) stop("folder_map.yaml not found: ", folder_map_path)

raw_task_rules <- yaml::read_yaml(task_rules_path)
tasks <- tryCatch(
  flatten_task_rules(raw_task_rules$workflows %||% list()),
  error = function(e) { cat("!! ", conditionMessage(e), "\n\n", sep = ""); NULL }
)
folder_map <- yaml::read_yaml(folder_map_path)$folder_map %||% list()

cat("=== Part 1: task_rules.yaml vs repo WDL task inventory ===\n\n")

if (is.null(tasks)) {
  cat("Skipped (duplicate task name error above must be fixed first).\n\n")
} else {
  declared_tasks <- vapply(names(tasks), function(nm) {
    wt <- tasks[[nm]]$wdl_task
    # `wdl_task: false` marks a synthetic identifier with no real task to
    # check against (YAML's `null` is indistinguishable from "key absent"
    # once parsed, so it can't be used as this sentinel - false can).
    if (isFALSE(wt)) return(NA_character_)
    if (is.null(wt)) nm else wt
  }, character(1))
  declared_tasks <- unique(declared_tasks[!is.na(declared_tasks) & nzchar(declared_tasks)])

  wdl_files <- list.files(file.path(REPO_ROOT, "hydrant"), pattern = "\\.wdl$", recursive = TRUE, full.names = TRUE)
  if (length(wdl_files) == 0) stop("no .wdl files found under ", file.path(REPO_ROOT, "hydrant"))

  extract_names <- function(files, line_pattern, capture_pattern) {
    out <- character(0)
    for (f in files) {
      lines <- suppressWarnings(readLines(f, warn = FALSE))
      hits <- grep(line_pattern, lines, value = TRUE, perl = TRUE)
      if (length(hits) == 0) next
      mm <- regmatches(hits, regexpr(capture_pattern, hits, perl = TRUE))
      out <- c(out, mm[nzchar(mm)])
    }
    unique(out)
  }

  defined_tasks <- extract_names(
    wdl_files,
    "^\\s*task\\s+panoply_[A-Za-z0-9_]+\\s*\\{",
    "panoply_[A-Za-z0-9_]+"
  )
  called_tasks <- extract_names(
    wdl_files,
    "\\bcall\\s+[A-Za-z0-9_.]*panoply_[A-Za-z0-9_]+",
    "panoply_[A-Za-z0-9_]+(?=\\s*(\\{|as\\b|$))"
  )

  cat(sprintf("task_rules.yaml:                %s\n", task_rules_path))
  cat(sprintf("tasks declared:                 %d\n", length(tasks)))
  cat(sprintf("wdl_task values checked:        %d (synthetic identifiers with wdl_task: false excluded)\n", length(declared_tasks)))
  cat(sprintf("task{...} definitions found:    %d (under hydrant/)\n", length(defined_tasks)))
  cat(sprintf("distinct tasks called anywhere: %d (under hydrant/)\n\n", length(called_tasks)))

  stale <- sort(setdiff(declared_tasks, defined_tasks))
  if (length(stale) > 0) {
    cat("!! task_rules.yaml entries with NO matching task definition anywhere under hydrant/\n")
    cat("   (typo, or the task was renamed/removed):\n")
    cat(paste0("   - ", stale), sep = "\n")
    cat("\n")
  } else {
    cat("OK: every declared task matches a real task definition.\n\n")
  }

  missing <- sort(setdiff(called_tasks, declared_tasks))
  if (length(missing) > 0) {
    cat("Tasks called somewhere under hydrant/ with no task_rules.yaml entry\n")
    cat("(best-effort - some are legitimately intermediate/passthrough tasks with no\n")
    cat(" results/ subfolder of their own; use judgement before adding a rule for each):\n")
    cat(paste0("   - ", missing), sep = "\n")
  } else {
    cat("OK: every called panoply_* task has a task_rules.yaml entry.\n")
  }
}

cat("\n=== Part 2: folder_map.yaml vs task_rules.yaml (hard check) ===\n\n")

# An entry is either {uses: task}, {each_child: task}, or a `group:` entry
# whose `children:` map is itself a set of uses/each_child entries (one
# level of nesting only, matching folder_map.yaml's schema) - collect every
# task name reachable this way, tagged with a "path" string for reporting.
collect_targets <- function(entry, path) {
  if (!is.null(entry$uses)) return(setNames(entry$uses, path))
  if (!is.null(entry$each_child)) return(setNames(entry$each_child, path))
  if (!is.null(entry$group)) {
    out <- character(0)
    for (fam in names(entry$children)) {
      out <- c(out, collect_targets(entry$children[[fam]], paste0(path, "/", fam)))
    }
    return(out)
  }
  character(0)
}
targets_by_path <- unlist(lapply(names(folder_map), function(nm) collect_targets(folder_map[[nm]], nm)))
targets <- unname(targets_by_path)

cat(sprintf("folder_map.yaml:      %s\n", folder_map_path))
cat(sprintf("folder entries:       %d\n", length(folder_map)))
cat(sprintf("distinct task targets: %d\n\n", length(unique(targets))))

if (is.null(tasks)) {
  cat("Skipped (fix the duplicate-task-name error in Part 1 first).\n")
} else {
  unknown <- sort(unique(setdiff(targets, names(tasks))))
  if (length(unknown) > 0) {
    cat("!! folder_map.yaml entries pointing at a task with NO task_rules.yaml entry:\n")
    for (i in seq_along(targets_by_path)) {
      if (targets_by_path[i] %in% unknown) {
        cat(sprintf("   - %s -> %s\n", names(targets_by_path)[i], targets_by_path[i]))
      }
    }
    cat("\n")
  } else {
    cat("OK: every folder_map.yaml target has a task_rules.yaml entry.\n\n")
  }

  unused <- sort(setdiff(names(tasks), unique(targets)))
  if (length(unused) > 0) {
    cat("Tasks declared in task_rules.yaml but never referenced by folder_map.yaml\n")
    cat("(fine if intentionally reserved for --task standalone use; otherwise dead config):\n")
    cat(paste0("   - ", unused), sep = "\n")
  } else {
    cat("OK: every declared task is referenced by folder_map.yaml.\n")
  }
}
