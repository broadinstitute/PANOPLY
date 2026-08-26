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
  "parameters", "ptmseaDB", "gseaDB", "clumpsFASTA"
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

# Inverse of wb_local_to_s3(): given an s3:// URI, returns the corresponding local path under
# wb_workbench_root() IF that URI is for this session's own bucket/project -- i.e. it's a
# well-formed local upload just quoted back in its S3 form (Manifold shows both). Returns NULL
# (not an error) for anything else -- a different project's URI, a malformed one, or
# S3_BUCKET/PROJECT_ID not being set -- since callers use this as a best-effort "maybe it's
# actually local" check, not a real S3 fetch (this never touches S3 itself).
wb_s3_to_local <- function(s3_uri) {
  bucket  <- Sys.getenv("S3_BUCKET")
  project <- Sys.getenv("PROJECT_ID")
  if (!nzchar(bucket) || !nzchar(project)) return(NULL)
  prefix <- paste0("s3://", bucket, "/research/projects/", project, "/")
  if (!startsWith(s3_uri, prefix)) return(NULL)
  file.path(wb_workbench_root(), substring(s3_uri, nchar(prefix) + 1))
}

### ===
### Sessions -- a session is a self-contained folder (mapped-file copies, subsets,
### master-parameters.yaml, and once named/saved, the built inputs.json) under sessions/
### alongside the deployed notebook/workbench-src (the notebook's own working directory -- see
### the file header comment). "current-session" is always the live, actively-edited session;
### naming and saving one (wb_save_session(), see sessions.r) snapshots it into its own named
### folder.
### This is also exactly the folder deploy-workbench.sh's --delete mirrors from the git repo:
### current-session/ is fair game for that (it only ever exists on the deployed side), but
### named sessions are explicitly excluded from that sync (see deploy-workbench.sh) so saving
### one actually protects it from a redeploy.
### ===

wb_sessions_root <- function() file.path(getwd(), "sessions")
wb_session_dir <- function(name = "current-session") file.path(wb_sessions_root(), name)

### ===
### Session state -- round-tripped as a single yaml file, replaces Terra's config.yaml restore
### ===

wb_state_path <- function() file.path(wb_session_dir(), ".panoply-session.yaml")

wb_default_state <- function() {
  list(
    typemap = list(),
    typemap_originals = list(),
    active_named_session = NULL,
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
    param_overrides = list(),
    subsets = list(),
    github_ref = GITHUB_REF,
    target_workflow = TARGET_WORKFLOW
  )
}

# Reads current-session/.panoply-session.yaml if present, without gating on file.exists()
# first: on some Manifold home directory mounts, file.exists() can still report TRUE for a
# just-deleted file (stale directory-listing metadata) even though the actual open then
# fails. Catching the read itself means "missing" and "exists but unreadable" are handled
# identically. Returns NULL (not an error) either way -- callers decide what that means.
wb_try_read_state <- function(path = wb_state_path()) {
  suppressWarnings(tryCatch(yaml::read_yaml(path), error = function(e) NULL))
}

wb_load_state <- function() {
  current <- wb_try_read_state()
  saved_names <- wb_list_saved_sessions()

  # Nothing to choose between -- a new session is the only sensible outcome, so skip the
  # menu entirely rather than making the user pick the one option that exists.
  if (is.null(current) && length(saved_names) == 0) {
    wb_msg("INFO", "No existing session found -- starting a new one.")
    wb_done()
    return(wb_default_state())
  }

  options <- character(0)
  if (!is.null(current)) options["resume"] <- "Use the current session (resume where you left off)"
  options["new"] <- "Start a new session (current-session/ will be reset)"
  if (length(saved_names) > 0) {
    options["load"] <- sprintf("Load a saved session (%s)", paste(saved_names, collapse = ", "))
  }

  # Falling back to "resume" (if a current session exists) or plain in-memory defaults is
  # never destructive -- used for both an outright quit and a declined destructive confirm.
  fall_back <- function() {
    if (!is.null(current)) {
      wb_msg("INFO", "Keeping the existing current session.")
      result <- modifyList(wb_default_state(), current)
    } else {
      result <- wb_default_state()
    }
    wb_done()
    result
  }

  menu <- paste0(
    "How would you like to start?\n",
    paste(sprintf("  %d) %s", seq_along(options), unname(options)), collapse = "\n"),
    "\n> "
  )
  idx <- wb_smart_readline(menu, valid = function(ch) {
    n <- suppressWarnings(as.integer(ch))
    if (is.na(n) || n < 1 || n > length(options)) {
      sprintf("Please enter a number from 1 to %d.", length(options))
    } else TRUE
  })
  if (is.null(idx)) return(fall_back())
  choice <- names(options)[as.integer(idx)]

  if (choice == "resume") {
    wb_msg("INFO", sprintf("Resuming the current session (%s).", wb_state_path()))
    result <- modifyList(wb_default_state(), current)
    wb_done()
    return(result)
  }

  if (choice == "new") {
    if (!is.null(current) &&
        !wb_confirm("This will erase current-session/ (mapped files, subsets, outputs). Continue?")) {
      return(fall_back())
    }
    if (unlink(wb_session_dir(), recursive = TRUE) != 0) {
      stop(sprintf("Failed to fully clear '%s' -- check for locked/in-use files and try again.", wb_session_dir()))
    }
    # dir.create() returns FALSE both on a genuine failure AND when the directory already
    # exists -- check dir.exists() too so the (harmless) "already there" case isn't mistaken
    # for a real error. This situation is unlikely after unlink() but technically possible.
    if (!dir.create(wb_session_dir(), recursive = TRUE) && !dir.exists(wb_session_dir())) {
      stop(sprintf("Failed to create a fresh '%s'.", wb_session_dir()))
    }
    wb_msg("INFO", "Starting a new session.")
    wb_done()
    return(wb_default_state())
  }

  # choice == "load"
  name <- wb_select_from_list(
    "Saved sessions:", saved_names,
    "Which saved session -- name or number: "
  )
  if (is.null(name)) return(fall_back())
  if (!is.null(current) &&
      !wb_confirm(sprintf("This will overwrite current-session/ with the contents of '%s'. Continue?", name))) {
    return(fall_back())
  }
  wb_msg("INFO", "Copying session files -- this can take a while for large GCTs, please wait...")
  wb_copy_session_tree(wb_session_dir(name), wb_session_dir(), exclude = SESSION_FINALIZED_FILES)
  loaded <- wb_try_read_state()
  wb_msg("INFO", sprintf("Loaded saved session '%s'.", name))
  result <- modifyList(wb_default_state(), loaded)
  wb_done()
  result
}

wb_save_state <- function(state, done = TRUE) {
  dir.create(dirname(wb_state_path()), showWarnings = FALSE, recursive = TRUE)
  wb_write_verified(function() yaml::write_yaml(state, wb_state_path()), wb_state_path())
  if (done) wb_done()
  invisible(state)
}

### ===
### Small shared utilities
### ===

# Shortens an absolute path under wherever sessions actually live (e.g.
# ".../sessions/current-session/inputs/foo.gct") down to the portion relative to that
# ("sessions/current-session/inputs/foo.gct") for display -- easier to read in a printed
# listing without losing the (still-unambiguous, since it's always relative to the same root)
# information. Deliberately derived from wb_sessions_root() rather than hardcoding getwd()
# directly -- they're the same in the real deployment, but this stays correct if that ever
# changes. Falls back to the path unchanged if it isn't under that root for some reason.
wb_display_path <- function(path) {
  root <- paste0(dirname(wb_sessions_root()), "/")
  if (startsWith(path, root)) return(substring(path, nchar(root) + 1))
  path
}

# General-purpose resolver for ANY user-supplied file path -- typed at a wb_smart_readline()
# prompt, or passed directly as a function argument (e.g. fasta_path=, zip_path=) -- so every
# such input point behaves the same way. s3:// URIs are the only case that need special
# handling -- translated via wb_s3_to_local() (the inverse of wb_local_to_s3()) if it's this
# project's own bucket/path (as Manifold shows for an uploaded file, e.g.
# "s3://<bucket>/research/projects/<id>/inputs/foo.fasta" -> "~/workbench/inputs/foo.fasta");
# NA is returned if that can't resolve it (different bucket/project, malformed URI,
# S3_BUCKET/PROJECT_ID unset), so the caller can give a clear "use a local path instead"
# message rather than silently treating "s3://..." itself as a (nonexistent) local path.
# Anything else is already a normal local path -- path.expand() handles "~", and relative/
# absolute paths need no further wrangling (file.exists() resolves them against the cwd as-is).
wb_resolve_user_path <- function(p) {
  if (startsWith(p, "s3://")) {
    local <- wb_s3_to_local(p)
    return(if (is.null(local)) NA_character_ else local)
  }
  path.expand(p)
}

# Generic wb_smart_readline(valid=...) check for a user-supplied file path: resolves it (see
# wb_resolve_user_path()), confirms it exists, and -- if `extensions` is given -- that it ends
# in one of them (case-insensitive, without the leading '.'). Returns TRUE to accept, or a
# message string for wb_smart_readline() to show and re-prompt.
wb_validate_user_file <- function(raw, extensions = NULL) {
  resolved <- wb_resolve_user_path(raw)
  if (is.na(resolved)) {
    return(paste(
      "Cloud (s3://) paths aren't supported directly here unless they're this project's own",
      "workbench bucket/path -- please upload the file under ~/workbench/ and enter its local",
      "path instead, try again."
    ))
  }
  if (!file.exists(resolved)) return(sprintf("No file found at '%s', try again.", resolved))
  if (!is.null(extensions)) {
    pattern <- paste0("\\.(", paste(extensions, collapse = "|"), ")$")
    if (!grepl(pattern, resolved, ignore.case = TRUE)) {
      return(sprintf("'%s' doesn't look like a .%s file, try again.", resolved, paste(extensions, collapse = "/.")))
    }
  }
  TRUE
}

### ===
### Standardized list-selection -- one mechanism for every "print a list, pick from it" prompt
### in workbench-src/. Accepts a comma-separated mix of literal values, 1-based indexes, and
### (where multiple selections make sense) index ranges like "3:5" -- freely intermixed, e.g.
### "Stage,2,4:6". An exact literal match always wins over index parsing, so a value that
### happens to look numeric (e.g. an annotation value literally "2") is never misread as an
### index -- this keeps the rule a single, predictable precedence order rather than
### context-dependent. Matching is case-sensitive, matching every existing literal-value
### validator elsewhere in this codebase (all plain %in% checks).
### ===

# Resolves one already-trimmed token against `options`: an exact literal match wins outright;
# otherwise a plain integer is a 1-based index, and (only when multi=TRUE) "a:b" is a range of
# indexes. Returns the 1-based indexes it names (length 1 for a literal/single-index match,
# length >1 for a range), or NULL if the token doesn't resolve to anything at all.
wb_resolve_selection_token <- function(token, options, multi) {
  if (token %in% options) return(which(options == token)[1])
  if (multi && grepl("^[0-9]+:[0-9]+$", token)) {
    bounds <- as.integer(strsplit(token, ":")[[1]])
    return(bounds[1]:bounds[2])
  }
  if (grepl("^[0-9]+$", token)) return(as.integer(token))
  NULL
}

# Resolves a full (comma-separated) selection string against `options`. Returns
# list(indices=) on success (1-based, in the order given, deduplicated) or list(error=) with a
# message describing the first problem found -- never both.
wb_resolve_selection <- function(raw, options, multi) {
  tokens <- wb_trim(strsplit(raw, ",")[[1]])
  if (!multi && length(tokens) > 1) {
    return(list(error = "Only one selection is allowed here, try again."))
  }
  indices <- integer(0)
  for (token in tokens) {
    resolved <- wb_resolve_selection_token(token, options, multi)
    if (is.null(resolved)) {
      return(list(error = sprintf("'%s' isn't a listed value, index, or range, try again.", token)))
    }
    if (any(resolved < 1 | resolved > length(options))) {
      return(list(error = sprintf("Index out of range (1-%d), try again.", length(options))))
    }
    indices <- c(indices, resolved)
  }
  list(indices = unique(indices))
}

wb_print_numbered_list <- function(header, options) {
  cat(header, "\n", sep = "")
  for (i in seq_along(options)) cat(sprintf("  %2d: %s\n", i, options[i]))
  flush.console()
}

# Prompts for a selection from `options` that the caller has ALREADY displayed (e.g. via a
# richer listing than a plain numbered one, like the color editor's swatch preview) -- see the
# section header above for accepted input forms. `extra_valid`, if given, is called with the
# resolved character vector of selected option(s) (length 1 unless multi=TRUE) and should
# return TRUE or an error message, layered on top of the generic value/index/range validation.
# Returns the selected option(s) as a character vector (length 1 unless multi=TRUE), or NULL if
# cancelled.
wb_prompt_selection <- function(options, prompt, multi = FALSE, cancel_msg = "No changes made.", extra_valid = NULL) {
  raw <- wb_smart_readline(
    prompt,
    valid = function(ch) {
      result <- wb_resolve_selection(ch, options, multi)
      if (!is.null(result$error)) return(result$error)
      if (!is.null(extra_valid)) {
        check <- extra_valid(options[result$indices])
        if (!isTRUE(check)) return(check)
      }
      TRUE
    },
    cancel_msg = cancel_msg
  )
  if (is.null(raw)) return(NULL)
  options[wb_resolve_selection(raw, options, multi)$indices]
}

# Prints a numbered list of `options` under `header`, then prompts via wb_prompt_selection() --
# the common case; use wb_prompt_selection() directly when the list has already been shown.
wb_select_from_list <- function(header, options, prompt, multi = FALSE, cancel_msg = "No changes made.", extra_valid = NULL) {
  wb_print_numbered_list(header, options)
  wb_prompt_selection(options, prompt, multi = multi, cancel_msg = cancel_msg, extra_valid = extra_valid)
}

### ===
### master-parameters.yaml overrides -- a single, standardized mechanism for "the uploaded data
### is fine, the CONFIGURED name/type/value was just wrong" corrections (PTM-SEA's
### seqwin_column, Clumps-PTM's accession_number_colname/FASTA_sep_type, MetaboAnalyst's
### meta_id_col/meta_id_type, etc.), replacing the earlier per-module pattern of rewriting the
### uploaded GCT in place to match the default. Overrides are recorded into `state` (there's no
### master-parameters.yaml file on disk yet at the point most of these run -- Finalize
### Parameters, which actually builds it, comes later in the notebook) and applied by
### wb_build_master_parameters_yaml() at build time, the same way state$toggles/state$cosmo_params/
### state$groups_cols already are.
### ===

# Records that master-parameters.yaml's nested key `path` (e.g.
# c("panoply_preprocess_gct", "seqwin_column")) should be set to `value` when the YAML is next
# built. `value = NULL` is a legitimate override (e.g. MetaboAnalyst's "use NULL for rid"
# convention), not a no-op -- see wb_set_nested_yaml_value() for how that's preserved through to
# the written file rather than silently deleting the key.
wb_set_param_override <- function(state, path, value) {
  state$param_overrides[[paste(path, collapse = ".")]] <- list(path = path, value = value)
  state
}

# Sets a nested key (`path`, a character vector of keys from root to leaf) inside a
# yaml-shaped nested list, creating intermediate levels as needed. Uses single-bracket
# `x[key] <- list(value)` for the final assignment rather than `x[[key]] <- value` --
# the latter DELETES the key entirely when value is NULL, instead of setting it to NULL, which
# would silently discard an explicit "use NULL for rid"-style override.
wb_set_nested_yaml_value <- function(x, path, value) {
  if (length(path) == 1) {
    x[path[1]] <- list(value)
    return(x)
  }
  child <- if (path[1] %in% names(x)) x[[path[1]]] else list()
  x[[path[1]]] <- wb_set_nested_yaml_value(child, path[-1], value)
  x
}

# Applies every recorded state$param_overrides entry onto a yaml-shaped nested list (see
# wb_set_param_override()). Order follows insertion order, so if two overrides ever target the
# same path within one session, the later one wins -- consistent with how re-running any other
# selection step in this notebook (groups, toggles, etc.) already overwrites earlier choices.
wb_apply_param_overrides <- function(yaml_list, param_overrides) {
  for (entry in param_overrides) {
    yaml_list <- wb_set_nested_yaml_value(yaml_list, entry$path, entry$value)
  }
  yaml_list
}

### ===
### s3fs-mount write verification -- ~/workbench is a FUSE (s3fs) mount, which can report a
### local write as successful (the VFS-level write()/close() returns 0) while the underlying
### S3 PUT happens asynchronously and fails silently, with no error surfaced back to R. This bit
### an inputs.json write: jsonlite::write_json() returned normally, but the file never actually
### existed in S3. Every write/copy that targets a path under the workbench mount should go
### through wb_write_verified() (or wb_retry() directly, for shapes write_verified() doesn't
### fit) rather than trusting a write call's own return value alone.
### ===

# Generic retry helper -- calls `fn()` up to max_attempts times. `fn` should throw (via stop())
# to signal a retryable failure; any other return is treated as success and returned
# immediately.
wb_retry <- function(fn, max_attempts = 3, retry_delay = 1, context = "operation") {
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(list(ok = TRUE, value = fn()), error = function(e) {
      last_error <<- e
      list(ok = FALSE, value = NULL)
    })
    if (result$ok) return(result$value)
    if (attempt < max_attempts) Sys.sleep(retry_delay)
  }
  stop(sprintf(
    "Failed %s after %d attempt(s)%s. ", context, max_attempts,
    if (!is.null(last_error)) sprintf(" (last error: %s)", conditionMessage(last_error)) else ""
  ), "This can happen transiently on the s3fs-backed workbench mount.")
}

# A bare file.exists() after a write can't distinguish "just written" from "stale file already
# there from before" when OVERWRITING an existing path -- exactly the inputs.json case, which
# almost always already exists. This additionally requires the mtime to have advanced past
# `since` (a small tolerance absorbs mtime-rounding/clock-skew on the mount, not genuine
# staleness).
wb_write_landed <- function(path, since) {
  file.exists(path) && file.mtime(path) >= (since - 1)
}

# Relative file paths under `dir`, recursively, sorted -- used to structurally verify a
# directory copy actually landed in full (see wb_copy_session_tree() in sessions.r), since
# directory mtimes on the s3fs mount aren't a reliable freshness signal the way a plain file's
# mtime is.
wb_recursive_relpaths <- function(dir) {
  sort(list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = FALSE))
}

# Wraps a single-file write with retry + landing verification (wb_write_landed()) -- `write_fn`
# performs the write to `path` as a side effect; its return value is ignored, only whether it
# throws and whether `path` actually shows a fresh mtime afterward.
wb_write_verified <- function(write_fn, path, max_attempts = 3, retry_delay = 1) {
  wb_retry(function() {
    since <- Sys.time()
    write_fn()
    if (!wb_write_landed(path, since)) {
      stop("write call reported success, but the file's modification time never advanced")
    }
    path
  }, max_attempts = max_attempts, retry_delay = retry_delay, context = sprintf("writing '%s'", path))
}

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

# Printed by every top-level wb_*() function the notebook calls directly, right before it
# returns, so it's unambiguous in the cell output that the function actually finished (as
# opposed to still running, or having silently stopped partway through).
wb_done <- function() { cat("=== DONE ===\n"); flush.console() }

wb_trim <- function(x) gsub("^\\s+|\\s+$", "", x)

# Recognized at any wb_smart_readline() prompt (case-insensitive) to back out of the current
# step cleanly, instead of being stuck in a validation loop with no escape but a kernel
# interrupt. Ported from build-config.r's exit_commands/valid_choice()/smart_readline().
WB_EXIT_COMMANDS <- c("q", "quit", "exit", "cancel")

wb_smart_readline <- function(prompt, valid = NULL, allow_empty = FALSE, cancel_msg="No changes made.") {
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
      wb_msg("CANCELLED", cancel_msg)
      return(NULL)
    }
    if (!allow_empty && !nzchar(choice)) {
      cat("Input cannot be empty (or type 'quit' to cancel). Please try again.\n")
      flush.console()
      next
    }
    result <- if (is.null(valid)) TRUE else valid(choice)
    if (isTRUE(result)) return(choice)
    cat(if (is.character(result)) result else "Invalid input (or type 'quit' to cancel).", "\n")
    flush.console()
  }
}

wb_confirm <- function(prompt, ...) {
  choice <- wb_smart_readline(
    paste0(prompt, " (y/n): "),
    valid = function(ch) if (tolower(ch) %in% c("y", "yes", "n", "no")) TRUE else "Please answer y or n (or 'quit' to cancel).",
    ...
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
  result <- length(still_missing) == 0
  wb_done()
  invisible(result)
}

### ===
### Load the rest of the module
### ===

source("workbench-src/sessions.r")
source("workbench-src/inputs.r")
source("workbench-src/groups.r")
source("workbench-src/subsets.r")
source("workbench-src/parameters.r")
source("workbench-src/wdl.r")
