# Session save/load/copy mechanics. "current-session" is the live, actively-edited
# session (see config.r for wb_sessions_root()/wb_session_dir()/wb_state_path()); named
# sessions under sessions/<name>/ are point-in-time snapshots of it.

wb_copy_session_tree <- function(from_dir, to_dir) {
  if (!dir.exists(from_dir)) stop(sprintf("Session directory not found: %s", from_dir))
  parent <- dirname(to_dir)
  dir.create(parent, showWarnings = FALSE, recursive = TRUE)
  # Copy into a temp directory *alongside* to_dir (same filesystem, so the final
  # file.rename() is atomic) before touching the real destination -- a mid-copy failure
  # (large GCTs, disk full, an interrupted kernel) then never leaves to_dir half-overwritten.
  staging_dir <- tempfile("session-copy-", tmpdir = parent)
  dir.create(staging_dir)
  entries <- list.files(from_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(entries) > 0) {
    ok <- file.copy(entries, staging_dir, recursive = TRUE)
    if (!all(ok)) {
      unlink(staging_dir, recursive = TRUE)
      stop(sprintf("Failed to copy into the session: %s", paste(basename(entries[!ok]), collapse = ", ")))
    }
  }
  if (dir.exists(to_dir)) unlink(to_dir, recursive = TRUE)
  if (!file.rename(staging_dir, to_dir)) {
    stop(sprintf("Failed to move the copied session into place at '%s'.", to_dir))
  }
  invisible(to_dir)
}

wb_copy_into_session <- function(source_path, max_attempts = 3, retry_delay = 1) {
  dest_dir <- file.path(wb_session_dir(), "inputs")
  if (!dir.exists(dest_dir)) dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE) # only create dir if dir exists
  dest_path <- file.path(dest_dir, basename(source_path))
  # attempt copy `max_attempt` times
  for (attempt in seq_len(max_attempts)) {
    ok <- suppressWarnings(file.copy(source_path, dest_path, overwrite = TRUE))
    if (ok) return(dest_path) # return if successful
    if (attempt < max_attempts) Sys.sleep(retry_delay) # try again after small delay
  }
  # print failure message if 
  stop(sprintf(
    "Failed to copy '%s' into the session (destination: '%s') after %d attempt(s). ",
    source_path, dest_path, max_attempts),
    "This can happen transiently on the s3fs-backed workbench mount.")
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
