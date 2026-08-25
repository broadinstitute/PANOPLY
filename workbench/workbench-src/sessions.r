# Session save/load/copy mechanics. "current-session" is the live, actively-edited
# session (see config.r for wb_sessions_root()/wb_session_dir()/wb_state_path()); named
# sessions under sessions/<name>/ are point-in-time snapshots of it.

# inputs.json (see wb_update_inputs_json_for_subset()) is only ever written directly into a
# *named* session -- never into current-session/. Excluded here (on the load direction
# specifically -- see wb_load_state()) so loading a saved session never pulls a copy of it back
# into current-session/, where it would just be a stale, untracked file that happens to look
# live. master-parameters.yaml is NOT in this list -- unlike inputs.json, it's built directly
# in current-session/ (see wb_build_master_parameters_yaml()) and carried into a named session
# the same way subsets/mapped inputs are, so it should round-trip on load just like those do.
SESSION_FINALIZED_FILES <- c("inputs.json", "inputs.json.bak")

wb_copy_session_tree <- function(from_dir, to_dir, exclude = character(0), max_attempts = 3, retry_delay = 1) {
  if (!dir.exists(from_dir)) stop(sprintf("Session directory not found: %s", from_dir))
  parent <- dirname(to_dir)
  dir.create(parent, showWarnings = FALSE, recursive = TRUE)
  entries <- list.files(from_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  entries <- entries[!basename(entries) %in% exclude]

  # Copy into a temp directory *alongside* to_dir (same filesystem, so the final
  # file.rename() is atomic) before touching the real destination -- a mid-copy failure
  # (large GCTs, disk full, an interrupted kernel) then never leaves to_dir half-overwritten.
  # The whole copy+rename is retried (not just re-checked) on a landing failure, since -- like
  # any write to the s3fs-backed mount -- file.copy()/file.rename() reporting success doesn't
  # guarantee the underlying S3 write actually happened (see wb_retry(), config.r). Verifying
  # BEFORE deleting an existing to_dir also means a failed/incomplete copy never destroys a
  # good, already-saved session in the process of trying to replace it.
  wb_retry(function() {
    staging_dir <- tempfile("session-copy-", tmpdir = parent)
    on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)
    dir.create(staging_dir)
    if (length(entries) > 0) {
      ok <- file.copy(entries, staging_dir, recursive = TRUE)
      if (!all(ok)) stop(sprintf("Failed to copy into the session: %s", paste(basename(entries[!ok]), collapse = ", ")))
    }
    staged <- wb_recursive_relpaths(staging_dir)
    if (length(entries) > 0 && length(staged) == 0) {
      stop("staged copy is empty despite file.copy() reporting success")
    }
    if (dir.exists(to_dir)) unlink(to_dir, recursive = TRUE)
    if (!file.rename(staging_dir, to_dir)) {
      stop(sprintf("Failed to move the copied session into place at '%s'.", to_dir))
    }
    if (!identical(wb_recursive_relpaths(to_dir), staged)) {
      stop("moved session tree doesn't match what was staged")
    }
  }, max_attempts = max_attempts, retry_delay = retry_delay, context = sprintf("copying the session into '%s'", to_dir))

  invisible(to_dir)
}

wb_copy_into_session <- function(source_path, max_attempts = 3, retry_delay = 1) {
  dest_dir <- file.path(wb_session_dir(), "inputs")
  if (!dir.exists(dest_dir)) dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE) # only create dir if dir exists
  dest_path <- file.path(dest_dir, basename(source_path))
  wb_write_verified(function() {
    if (!suppressWarnings(file.copy(source_path, dest_path, overwrite = TRUE))) {
      stop("file.copy() reported failure")
    }
  }, dest_path, max_attempts = max_attempts, retry_delay = retry_delay)
  dest_path
}

wb_list_saved_sessions <- function() {
  root <- wb_sessions_root()
  if (!dir.exists(root)) return(character(0))
  setdiff(list.dirs(root, recursive = FALSE, full.names = FALSE), "current-session")
}

wb_save_session <- function(state, name = NULL) {
  valid_name <- function(ch) {
    if (identical(ch, "current-session")) return("'current-session' is reserved -- choose a different name.")
    if (!grepl("^[A-Za-z0-9_.-]+$", ch)) return("Session names may only contain letters, numbers, '-', '_', and '.', try again.")
    TRUE
  }
  if (is.null(name)) {
    name <- wb_smart_readline("Name for this session: ", valid = valid_name)
    if (is.null(name)) {
      wb_msg("CANCELLED", "Session not saved.")
      wb_done()
      return(state)
    }
  } else if (!isTRUE(valid_name(name))) {
    stop(valid_name(name))
  }
  if (dir.exists(wb_session_dir(name)) &&
      !wb_confirm(sprintf("A saved session named '%s' already exists. Overwrite it?", name))) {
    wb_msg("CANCELLED", "Session not saved.")
    wb_done()
    return(state)
  }
  state$active_named_session <- name
  state <- wb_save_state(state, done = FALSE)
  wb_msg("INFO", "Copying session files -- this can take a while for large GCTs, please wait...")
  wb_copy_session_tree(wb_session_dir(), wb_session_dir(name))
  wb_msg("INFO", sprintf("Session saved as '%s' (%s).", name, wb_session_dir(name)))
  wb_done()
  state
}

# A fast alternative to wb_load_state()'s "load a saved session" option, for when all you want
# is to regenerate inputs.json for an already-saved session -- that's the one artifact that's
# always written directly into the named session regardless of current-session/'s contents (see
# wb_update_inputs_json_for_subset()), so nothing needs to be copied anywhere to do it. Reads
# the saved session's own state file directly, without touching current-session/ at all.
#
# IMPORTANT: only use the returned state to call wb_update_inputs_json_for_subset() (or inspect
# it). Anything that writes based on wb_session_dir() with no argument -- wb_create_subset(),
# wb_load_and_map_inputs(), wb_build_master_parameters_yaml(), wb_select_groups(), etc. -- always
# targets current-session/, which this deliberately leaves untouched, so using this state for
# further editing would silently write to the wrong place. For ongoing editing of a saved
# session, use wb_load_state()'s "load a saved session" option instead, which fully (and more
# slowly) copies it into current-session/ first.
wb_open_saved_session <- function(name = NULL) {
  saved_names <- wb_list_saved_sessions()
  if (length(saved_names) == 0) stop("No saved sessions found -- run wb_save_session() first.")
  if (is.null(name)) {
    name <- wb_smart_readline(
      sprintf("Which saved session? (%s): ", paste(saved_names, collapse = ", ")),
      valid = function(ch) if (ch %in% saved_names) TRUE else "Not a known saved session, try again."
    )
    if (is.null(name)) {
      wb_msg("CANCELLED", "No session opened.")
      wb_done()
      return(invisible(NULL))
    }
  } else if (!(name %in% saved_names)) {
    stop(sprintf("'%s' is not a known saved session (%s).", name, paste(saved_names, collapse = ", ")))
  }
  loaded <- wb_try_read_state(path = file.path(wb_session_dir(name), ".panoply-session.yaml"))
  if (is.null(loaded)) stop(sprintf("Could not read the state file for saved session '%s'.", name))
  state <- modifyList(wb_default_state(), loaded)
  wb_msg("INFO", sprintf(
    paste("Opened '%s' directly, without copying it into current-session/ -- use this only to",
          "regenerate inputs.json for it. For ongoing editing, load it via wb_load_state() instead."),
    name
  ))
  wb_done()
  state
}
