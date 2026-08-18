# GitHub fetch, WDL input parsing, and inputs.json build/update.

wb_gh_fetch <- function(path, ref = GITHUB_REF, raw = FALSE, repo = GITHUB_REPO) {
  api_path <- sprintf("repos/%s/contents/%s?ref=%s", repo, path, utils::URLencode(ref, reserved = TRUE))

  if (nzchar(Sys.which("gh"))) {
    result <- tryCatch({
      args <- if (raw) c("api", "-H", "Accept: application/vnd.github.raw", api_path) else c("api", api_path)
      wb_run_cmd("gh", args)
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
    wb_msg("WARNING", "gh api call failed; falling back to an unauthenticated HTTPS request.")
  }

  url <- if (raw) {
    sprintf("https://raw.githubusercontent.com/%s/%s/%s", repo, ref, path)
  } else {
    sprintf("https://api.github.com/repos/%s/contents/%s?ref=%s", repo, path, utils::URLencode(ref, reserved = TRUE))
  }
  readLines(url, warn = FALSE)
}

wb_gh_fetch_binary <- function(path, out_path, ref = GITHUB_REF, repo = GITHUB_REPO) {
  # For binary files (e.g. .qs), wb_gh_fetch()'s text-line stdout capture would corrupt the
  # content -- redirect straight to a file instead, which preserves bytes exactly.
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  if (nzchar(Sys.which("gh"))) {
    api_path <- sprintf("repos/%s/contents/%s?ref=%s", repo, path, utils::URLencode(ref, reserved = TRUE))
    args <- shQuote(c("api", "-H", "Accept: application/vnd.github.raw", api_path))
    status <- tryCatch(system2("gh", args, stdout = out_path, stderr = FALSE), error = function(e) 1L)
    if (identical(status, 0L) && file.exists(out_path) && file.info(out_path)$size > 0) return(out_path)
    wb_msg("WARNING", "gh api call failed; falling back to an unauthenticated HTTPS request.")
  }

  url <- sprintf("https://raw.githubusercontent.com/%s/%s/%s", repo, ref, path)
  utils::download.file(url, destfile = out_path, mode = "wb", quiet = TRUE)
  out_path
}

wb_list_github_workflows <- function(ref = GITHUB_REF, repo = GITHUB_REPO) {
  listing <- wb_gh_fetch("hydrant/workflows", ref = ref, raw = FALSE, repo = repo)
  entries <- jsonlite::fromJSON(paste(listing, collapse = "\n"), simplifyDataFrame = FALSE)
  names <- vapply(entries, function(e) if (identical(e$type, "dir")) e$name else NA_character_, character(1))
  result <- sort(names[!is.na(names)])
  wb_done()
  result
}

wb_fetch_workflow_wdl <- function(workflow_name, ref = GITHUB_REF, repo = GITHUB_REPO, use_cache = TRUE) {
  cache_dir <- file.path(wb_workbench_root(), ".wdl_cache", ref)
  cache_path <- file.path(cache_dir, paste0(workflow_name, ".wdl"))
  if (use_cache && file.exists(cache_path)) {
    wb_msg("INFO", sprintf(
      "Using cached %s.wdl from a previous fetch (%s). Delete this file, or pass use_cache=FALSE, to force a fresh fetch.",
      workflow_name, cache_path
    ))
    return(paste(readLines(cache_path, warn = FALSE), collapse = "\n"))
  }

  path <- sprintf("hydrant/workflows/%s/%s.wdl", workflow_name, workflow_name)
  text <- paste(wb_gh_fetch(path, ref = ref, raw = TRUE, repo = repo), collapse = "\n")

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(text, cache_path)
  text
}

### ===
### WDL input parsing: brace-depth scan (draft-2) + input{} block (version 1.x)
### ===

wb_parse_wdl_inputs <- function(wdl_text, workflow_name) {
  lines <- strsplit(gsub("\r", "", wdl_text), "\n")[[1]]
  strip_comment <- function(line) sub("#.*$", "", line)

  wf_start <- grep(sprintf("^\\s*workflow\\s+%s\\s*\\{", workflow_name), lines)
  if (length(wf_start) == 0) {
    stop(sprintf("Could not find 'workflow %s {' in the WDL source.", workflow_name))
  }
  wf_start <- wf_start[1]

  depth <- 0
  # mode: "scan" (draft-2 top-level scan) -> becomes "input_block" once an explicit `input {`
  # is seen, then "done" the moment that block closes. Declarations are only ever collected
  # in "scan" (depth 1) or "input_block" (depth 2) mode -- "done" means stop looking entirely,
  # since a versioned WDL's `output { }` block sits at that same depth-2 nesting level and
  # must not be mistaken for more inputs once the real input{} block has already closed.
  mode <- "scan"
  results <- list()

  for (i in seq(wf_start, length(lines))) {
    line <- strip_comment(lines[i])
    trimmed <- trimws(line)
    opens  <- lengths(regmatches(line, gregexpr("\\{", line)))
    closes <- lengths(regmatches(line, gregexpr("\\}", line)))
    pre_depth <- depth

    if (i == wf_start) {
      depth <- depth + opens - closes
      next
    }
    if (mode == "done") break

    if (mode == "scan" && pre_depth == 1 && grepl("^input\\s*\\{", trimmed)) {
      mode <- "input_block"
      depth <- depth + opens - closes
      next
    }

    target_depth <- if (mode == "input_block") 2 else 1
    if (pre_depth == target_depth && nzchar(trimmed) &&
        !grepl("^(call|scatter|if|output|meta|parameter_meta|\\})\\b", trimmed)) {
      m <- regmatches(trimmed, regexec(
        "^([A-Za-z_][A-Za-z0-9_]*(?:\\[[^]]*\\])?\\??)\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(=.*)?$",
        trimmed, perl = TRUE
      ))[[1]]
      if (length(m) == 4) {
        wdl_type <- m[2]
        name <- m[3]
        has_default <- nzchar(m[4])
        is_input <- if (mode == "input_block") TRUE else !has_default
        if (is_input) {
          results[[length(results) + 1]] <- list(
            name = name, wdl_type = wdl_type,
            optional = grepl("\\?$", wdl_type),
            is_file = identical(sub("\\?$", "", wdl_type), "File"),
            has_default = has_default
          )
        }
      }
    }

    depth <- depth + opens - closes

    if (mode == "input_block" && depth == 1) mode <- "done"
    if (mode == "scan" && depth <= 0) break
  }

  if (length(results) == 0) {
    return(data.frame(name = character(0), wdl_type = character(0), optional = logical(0),
                      is_file = logical(0), has_default = logical(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(results, function(r) {
    data.frame(name = r$name, wdl_type = r$wdl_type, optional = r$optional,
               is_file = r$is_file, has_default = r$has_default, stringsAsFactors = FALSE)
  }))
}

### ===
### Semantic-role aliasing, for mapping WDL input names -> known file/parameter roles
### ===

INPUT_ALIASES <- list(
  yaml = c("yaml", "yaml_file"),
  job_id = c("job_id", "job_identifier", "output_prefix", "label"),
  prote_ome = "prote_ome", phospho_ome = "phospho_ome", acetyl_ome = "acetyl_ome",
  ubiquityl_ome = "ubiquityl_ome", nglyco_ome = "nglyco_ome", methyl_ome = "methyl_ome",
  metabol_ome = "metabol_ome", rna_data = "rna_data", cna_data = "cna_data",
  groups_file = "groups_file", groups_file_nmf = "groups_file_nmf",
  groups_file_metaboanlayst = "groups_file_metaboanlayst",
  groups_file_clumpsptm = "groups_file_clumpsptm",
  geneset_db = "geneset_db", ptm_db = "ptm_db"
)

ROLE_TO_SUBSET_CATEGORY <- c(
  prote_ome = "proteome", phospho_ome = "phosphoproteome", acetyl_ome = "acetylome",
  ubiquityl_ome = "ubiquitylome", nglyco_ome = "nglycoproteome", methyl_ome = "methylation",
  metabol_ome = "metabolome", rna_data = "rna", cna_data = "cna",
  groups_file = "groups", groups_file_clumpsptm = "groups_clumpsptm"
)
ROLE_TO_STATIC_CATEGORY <- c(geneset_db = "gseaDB", ptm_db = "ptmseaDB")

TOGGLE_ALIASES <- c(
  run_cmap = "run_cmap", run_mo_nmf = "run_mo_nmf", run_so_nmf = "run_so_nmf",
  run_ptmsea = "run_ptmsea", run_clumps = "run_clumpsptm", run_metab = "run_metab",
  normalizeProteomics = "normalize_proteomics", filterProteomics = "filter_proteomics"
)

wb_map_semantic_role <- function(spec_name) {
  for (role in names(INPUT_ALIASES)) if (spec_name %in% INPUT_ALIASES[[role]]) return(role)
  NA_character_
}

### ===
### inputs.json build + surgical update
### ===

wb_build_inputs_json <- function(state, subset_name, workflow_name = state$target_workflow %||% TARGET_WORKFLOW,
                                  github_ref = state$github_ref %||% GITHUB_REF) {
  specs <- wb_parse_wdl_inputs(wb_fetch_workflow_wdl(workflow_name, github_ref), workflow_name)
  subset_files <- wb_subset_files(state, subset_name)
  master_params_path <- file.path(wb_session_dir(state$active_named_session), "master-parameters.yaml")

  inputs <- list()
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    role <- wb_map_semantic_role(spec$name)
    key <- paste0(workflow_name, ".", spec$name)

    if (isTRUE(spec$is_file)) {
      local_path <- NULL
      if (!is.na(role)) {
        if (role == "yaml") {
          local_path <- master_params_path
        } else if (role %in% names(ROLE_TO_SUBSET_CATEGORY)) {
          local_path <- subset_files[[ROLE_TO_SUBSET_CATEGORY[[role]]]]
        } else if (role %in% names(ROLE_TO_STATIC_CATEGORY)) {
          local_path <- state$typemap[[ROLE_TO_STATIC_CATEGORY[[role]]]]
        }
      }
      if (!is.null(local_path) && !is.na(local_path) && file.exists(local_path)) {
        inputs[[key]] <- wb_local_to_s3(local_path)
      } else if (!spec$optional) {
        wb_msg("WARNING", sprintf("Required File input '%s' could not be resolved -- fill it in manually.", spec$name))
      }
      next
    }

    if (!is.na(role) && role == "job_id") {
      inputs[[key]] <- as.character(state$job_id %||% subset_name)
      next
    }

    toggle_name <- unname(TOGGLE_ALIASES[spec$name])
    if (!is.na(toggle_name)) {
      value <- state$toggles[[toggle_name]]
      if (!is.null(value)) {
        inputs[[key]] <- if (grepl("^Boolean", spec$wdl_type)) as.logical(value) else tolower(as.character(value))
      }
    }
  }

  inputs
}

wb_update_inputs_json_for_subset <- function(state, subset_name,
                                              workflow_name = state$target_workflow %||% TARGET_WORKFLOW,
                                              github_ref = state$github_ref %||% GITHUB_REF,
                                              existing_inputs_path = NULL,
                                              out_path = NULL) {
  # Requires a named session (not current-session) so this file's S3 paths -- once baked
  # into a submitted job -- can't be silently invalidated by later, unrelated work in
  # current-session. See wb_save_session() in sessions.r.
  if (is.null(state$active_named_session)) {
    stop("No named session found -- run `state <- wb_save_session(state, \"your-name\")` first ",
         "(see the Sessions section) before generating inputs.json.")
  }
  existing_inputs_path <- existing_inputs_path %||% file.path(wb_session_dir(state$active_named_session), "inputs.json")
  out_path <- out_path %||% existing_inputs_path
  fresh <- wb_build_inputs_json(state, subset_name, workflow_name, github_ref)

  if (is.null(existing_inputs_path) || !file.exists(existing_inputs_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    jsonlite::write_json(fresh, out_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
    wb_msg("INFO", sprintf("Wrote a new inputs.json to %s", out_path))
    wb_done()
    return(out_path)
  }

  backup_path <- paste0(existing_inputs_path, ".bak")
  if (!file.copy(existing_inputs_path, backup_path, overwrite = TRUE)) {
    stop(sprintf("Failed to back up '%s' to '%s' -- aborting without touching it.", existing_inputs_path, backup_path))
  }
  existing <- jsonlite::fromJSON(existing_inputs_path, simplifyVector = FALSE)

  specs <- wb_parse_wdl_inputs(wb_fetch_workflow_wdl(workflow_name, github_ref), workflow_name)
  file_specs <- specs[specs$is_file, , drop = FALSE]
  for (i in seq_len(nrow(file_specs))) {
    key <- paste0(workflow_name, ".", file_specs$name[i])
    existing[[key]] <- fresh[[key]]
  }

  jsonlite::write_json(existing, out_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  wb_msg("INFO", sprintf("Updated file-path inputs in %s for subset '%s' (backup at %s.bak)",
                        out_path, subset_name, existing_inputs_path))
  wb_done()
  out_path
}
