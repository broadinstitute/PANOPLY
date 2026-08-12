# PANOPLY Workbench Notebook -- helper module entry point.
# Sourced by PANOPLY-workbench-notebook.ipynb via source("workbench-src/config.r").
# All relative paths below (here and in the sibling files) are relative to the notebook's
# own working directory, not to this file's location -- the notebook and workbench-src/ are
# assumed to be deployed side by side.

### ===
### Data-category / annotation constants (mirrors panda-src/build-config.r)
### ===

CAT_MAP <- c(
  "proteome", "phosphoproteome", "acetylome", "ubiquitylome", "methylation",
  "nglycoproteome", "rna", "cna", "metabolome", "annotation", "groups",
  "parameters", "ptmseaDB", "gseaDB"
)
PROTEOME_TYPES <- c("proteome", "phosphoproteome", "acetylome", "ubiquitylome", "methylation", "nglycoproteome")
REQUIRED_COLS  <- c("Sample.ID", "Type")
IGNORE_COLS    <- c("Experiment", "Channel", "Participant", "QC.status")

### ===
### GitHub / workflow targeting -- edit these to change defaults
### ===

GITHUB_REPO     <- "broadinstitute/PANOPLY"
GITHUB_REF      <- "issue-githubWDL"          # bump to a release branch/tag once one exists
TARGET_WORKFLOW <- "panoply_unified_workflow" # workflow to build inputs.json for

### ===
### Workbench path <-> S3 translation
### ===

wb_workbench_root <- function() path.expand("~/workbench")

wb_local_to_s3 <- function(local_path) {
  bucket  <- Sys.getenv("S3_BUCKET")
  project <- Sys.getenv("PROJECT_ID")
  if (!nzchar(bucket) || !nzchar(project)) {
    stop("S3_BUCKET and/or PROJECT_ID environment variables are not set; ",
         "cannot translate a local workbench path to S3.")
  }
  root     <- normalizePath(wb_workbench_root(), mustWork = FALSE)
  abs_path <- normalizePath(local_path, mustWork = FALSE)
  if (!startsWith(abs_path, root)) {
    stop(sprintf("Path '%s' is not under the workbench root '%s'; cannot translate to S3.",
                  local_path, root))
  }
  rel <- sub("^/+", "", substring(abs_path, nchar(root) + 1))
  paste0("s3://", bucket, "/research/projects/", project, "/", rel)
}

### ===
### Session state -- round-tripped as a single yaml file, replaces Terra's config.yaml restore
### ===

wb_state_path <- function() file.path(wb_workbench_root(), ".panoply-session.yaml")

wb_default_state <- function() {
  list(
    typemap = list(),
    groups_cols = character(0),
    groups_cols_continuous = character(0),
    groups_colors = list(),
    toggles = list(
      normalize_proteomics = FALSE, filter_proteomics = FALSE,
      run_ptmsea = FALSE, run_metab = FALSE, run_clumpsptm = FALSE,
      run_cmap = FALSE, run_mo_nmf = FALSE, run_so_nmf = TRUE
    ),
    cosmo_params = list(run_cosmo = FALSE, sample_label = ""),
    job_id = NULL,
    subsets = list(),
    github_ref = GITHUB_REF,
    target_workflow = TARGET_WORKFLOW
  )
}

wb_load_state <- function() {
  path <- wb_state_path()
  if (file.exists(path)) {
    message("Loaded existing session state from ", path)
    return(modifyList(wb_default_state(), yaml::read_yaml(path)))
  }
  wb_default_state()
}

wb_save_state <- function(state) {
  dir.create(wb_workbench_root(), showWarnings = FALSE, recursive = TRUE)
  yaml::write_yaml(state, wb_state_path())
  invisible(state)
}

### ===
### Small shared utilities
### ===

wb_run_cmd <- function(cmd, args = character(0)) {
  # system2() builds a shell command line without quoting args itself (e.g. a header value
  # like "Accept: application/vnd.github.raw" would otherwise be word-split on the space) --
  # shQuote() every element so each is treated as a single shell token.
  out <- suppressWarnings(system2(cmd, shQuote(args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("Command failed (%s %s):\n%s", cmd, paste(args, collapse = " "),
                 paste(out, collapse = "\n")))
  }
  out
}

# flush.console() after every status message -- Jupyter/IRkernel can otherwise buffer cat()
# output instead of showing it right away, which is especially misleading right before a
# step that takes a while (looks like the cell has silently hung).
wb_msg <- function(type, ...) { cat(sprintf("[%s] %s\n", type, paste0(...))); flush.console() }

wb_trim <- function(x) gsub("^\\s+|\\s+$", "", x)

# Recognized at any wb_smart_readline() prompt (case-insensitive) to back out of the current
# step cleanly, instead of being stuck in a validation loop with no escape but a kernel
# interrupt. Ported from build-config.r's exit_commands/valid_choice()/smart_readline().
WB_EXIT_COMMANDS <- c("q", "quit", "exit", "cancel")

wb_smart_readline <- function(prompt, valid = NULL, allow_empty = FALSE) {
  # readline() that re-prompts until the response is valid, and lets the user type an exit
  # command to cancel out at any point -- returns NULL in that case, so callers can just
  # check is.null(result) rather than each needing their own escape hatch.
  #
  # `valid`, if given, is called on the trimmed input and should return TRUE (accept),
  # FALSE (reject with a generic message), or a character string (reject with THAT specific
  # message) -- e.g. valid = function(x) if (x %in% choices) TRUE else "Not a valid choice."
  repeat {
    choice <- wb_trim(readline(prompt))
    flush.console()
    if (tolower(choice) %in% WB_EXIT_COMMANDS) {
      wb_msg("CANCELLED", "No changes made.")
      return(NULL)
    }
    if (!allow_empty && !nzchar(choice)) {
      cat("Input cannot be empty (or type 'quit' to cancel). Please try again.\n")
      next
    }
    result <- if (is.null(valid)) TRUE else valid(choice)
    if (isTRUE(result)) return(choice)
    cat(if (is.character(result)) result else "Invalid input (or type 'quit' to cancel).", "\n")
  }
}

wb_confirm <- function(prompt) {
  choice <- wb_smart_readline(
    paste0(prompt, " (y/n): "),
    valid = function(ch) if (tolower(ch) %in% c("y", "yes", "n", "no")) TRUE else "Please answer y or n (or 'quit' to cancel)."
  )
  if (is.null(choice)) return(FALSE)  # quitting a y/n question is treated as declining
  tolower(choice) %in% c("y", "yes")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

### ===
### Vendored proteomics-Rutil scripts (map_id, set_annot_colors) -- see workbench-src/r-utils/
### ===

wb_source_rutil_vendor <- function(dir = "workbench-src/r-utils") {
  tryCatch({
    source(file.path(dir, "color-mod-utils.r"))
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(dir)
    source("map-to-genes.r")
    TRUE
  }, error = function(e) {
    wb_msg("WARNING", "Could not load vendored proteomics-Rutil scripts (", conditionMessage(e), "). ",
           "Gene/protein-ID conversion and default color assignment will be unavailable until this is resolved.")
    FALSE
  })
}

### ===
### One-time environment setup -- installs any missing packages, then loads the vendored
### proteomics-Rutil scripts. Called explicitly by the notebook's Setup cell (wb_setup()),
### not automatically on source(), so a fresh environment gets a chance to install packages
### before anything that depends on them runs.
### ===

wb_setup <- function() {
  # cat() output can sit in Jupyter/IRkernel's console buffer instead of displaying right
  # away -- particularly unhelpful right before a step that can take a while, since it looks
  # like the cell is silently hanging rather than working. flush.console() forces whatever's
  # been printed so far to actually show up before moving on (same reason the old
  # panda-src/build-config.r sprinkled it after every user-facing print).
  cat("Checking installed packages...\n"); flush.console()

  cran_pkgs <- c("yaml", "jsonlite", "RColorBrewer", "dplyr", "khroma", "qs", "BiocManager")
  bioc_pkgs <- c("cmapR", "org.Hs.eg.db", "EnsDb.Hsapiens.v79")
  all_pkgs  <- c(cran_pkgs, bioc_pkgs)

  is_missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  missing_pkgs <- is_missing(all_pkgs)

  if (length(missing_pkgs) > 0) {
    cat("Installing missing packages (first run only -- this can take a while):\n -",
        paste(missing_pkgs, collapse = ", "), "\n")
    flush.console()

    conda_bin <- Sys.which("mamba"); if (!nzchar(conda_bin)) conda_bin <- Sys.which("conda")
    if (nzchar(conda_bin)) {
      # Target the conda env that's actually backing THIS running R session explicitly --
      # installing without --prefix depends on ambient CONDA_PREFIX/activation state, which
      # doesn't reliably match the Jupyter kernel's environment and can silently install
      # packages where this R session will never see them (no error, just doesn't help).
      r_home_prefix <- dirname(dirname(R.home()))
      conda_prefix <- if (dir.exists(file.path(r_home_prefix, "conda-meta"))) {
        r_home_prefix
      } else {
        Sys.getenv("CONDA_PREFIX")
      }
      conda_names <- ifelse(missing_pkgs %in% bioc_pkgs,
                            paste0("bioconductor-", tolower(missing_pkgs)),
                            paste0("r-", tolower(missing_pkgs)))
      # --freeze-installed keeps every already-installed package (R itself included) as a
      # hard constraint, so the solver can only ever *add* the requested package -- never
      # upgrade R or anything else to satisfy it. Without this, an unpinned "install
      # r-qs" can be solved by bumping R itself (and everything built against it) to
      # whatever newer R version has the newest builds, silently trashing the environment
      # the notebook is actually running in.
      args <- c("install", "-y", "--freeze-installed", "-c", "conda-forge", "-c", "bioconda")
      if (nzchar(conda_prefix)) args <- c(args, "--prefix", conda_prefix)
      # One call per package, not one batched call for all of them: if any single name
      # isn't resolvable on these channels (e.g. some CRAN-only packages were never
      # published to conda-forge/bioconda under an "r-<pkg>" name), mamba aborts the
      # *whole* transaction and installs nothing -- silently blocking every other
      # package in the batch that would otherwise have resolved fine on its own.
      for (name in conda_names) system2(conda_bin, c(args, name))
      missing_pkgs <- is_missing(all_pkgs)  # re-check -- don't trust exit status alone
    }

    still_cran <- intersect(missing_pkgs, cran_pkgs)
    still_bioc <- intersect(missing_pkgs, bioc_pkgs)
    if (length(still_cran) > 0) install.packages(still_cran)
    if (length(still_bioc) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(still_bioc, update = FALSE, ask = FALSE)
    }
  }

  still_missing <- is_missing(all_pkgs)
  if (length(still_missing) > 0) {
    warning("Setup incomplete -- still missing: ", paste(still_missing, collapse = ", "),
            ". Some features may not work until these are installed manually and wb_setup() ",
            "is re-run.", call. = FALSE)
  } else {
    cat("All required packages are available.\n")
  }
  flush.console()

  invisible(wb_source_rutil_vendor())
  invisible(length(still_missing) == 0)
}

### ===
### Load the rest of the module
### ===

source("workbench-src/inputs.r")
source("workbench-src/groups.r")
source("workbench-src/subsets.r")
source("workbench-src/parameters.r")
source("workbench-src/wdl.r")
