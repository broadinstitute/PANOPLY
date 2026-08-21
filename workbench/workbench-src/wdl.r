# GitHub fetch, WDL input parsing, and inputs.json build/update.

wb_gh_fetch <- function(path, ref = GITHUB_REF, raw = FALSE, repo = GITHUB_REPO) {
  # An anonymous request works for any public repo (which is all this notebook ever targets)
  # and, for raw = TRUE, hits raw.githubusercontent.com's CDN rather than the much more tightly
  # rate-limited api.github.com -- so it's tried first, with 'gh' (which most notebook users
  # won't have authenticated) only as a fallback for the rare case the anonymous request fails
  # (rate limiting, a private fork/branch, etc).
  url <- if (raw) {
    sprintf("https://raw.githubusercontent.com/%s/%s/%s", repo, ref, path)
  } else {
    sprintf("https://api.github.com/repos/%s/contents/%s?ref=%s", repo, path, utils::URLencode(ref, reserved = TRUE))
  }
  result <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
  if (!is.null(result)) return(result)

  if (!nzchar(Sys.which("gh"))) {
    stop(sprintf("Failed to fetch '%s' from GitHub (anonymous request failed, and no 'gh' CLI is available to retry with).", path))
  }
  wb_msg("WARNING", "Anonymous GitHub request failed; retrying via the 'gh' CLI (rate-limited, or a private repo/branch?).")
  api_path <- sprintf("repos/%s/contents/%s?ref=%s", repo, path, utils::URLencode(ref, reserved = TRUE))
  args <- if (raw) c("api", "-H", "Accept: application/vnd.github.raw", api_path) else c("api", api_path)
  wb_run_cmd("gh", args)
}

wb_gh_fetch_binary <- function(path, out_path, ref = GITHUB_REF, repo = GITHUB_REPO) {
  # For binary files (e.g. .qs), wb_gh_fetch()'s text-line stdout capture would corrupt the
  # content -- redirect straight to a file instead, which preserves bytes exactly. Anonymous
  # first, 'gh' only as a fallback -- see wb_gh_fetch() above for why.
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

  url <- sprintf("https://raw.githubusercontent.com/%s/%s/%s", repo, ref, path)
  ok <- tryCatch({
    utils::download.file(url, destfile = out_path, mode = "wb", quiet = TRUE)
    file.exists(out_path) && file.info(out_path)$size > 0
  }, error = function(e) FALSE)
  if (ok) return(out_path)

  if (!nzchar(Sys.which("gh"))) {
    stop(sprintf("Failed to fetch '%s' from GitHub (anonymous request failed, and no 'gh' CLI is available to retry with).", path))
  }
  wb_msg("WARNING", "Anonymous GitHub request failed; retrying via the 'gh' CLI (rate-limited, or a private repo/branch?).")
  api_path <- sprintf("repos/%s/contents/%s?ref=%s", repo, path, utils::URLencode(ref, reserved = TRUE))
  args <- shQuote(c("api", "-H", "Accept: application/vnd.github.raw", api_path))
  status <- system2("gh", args, stdout = out_path, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(out_path) || file.info(out_path)$size == 0) {
    stop(sprintf("Failed to fetch '%s' from GitHub via both an anonymous request and the 'gh' CLI.", path))
  }
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

# On Manifold, workflows are targeted by GitHub branch (GITHUB_REF), not a pinned release --
# the WDL at a given (ref, workflow_name) can change at any time as that branch is pushed to,
# so a fetch is fresh by default. use_cache = TRUE is an explicit opt-in for anyone who wants
# to avoid repeated fetches (e.g. iterating quickly, or working offline against a copy already
# on disk) and is willing to accept it may not reflect the branch's current state.
wb_workflow_wdl_path <- function(workflow_name) sprintf("hydrant/workflows/%s/%s.wdl", workflow_name, workflow_name)

wb_fetch_workflow_wdl <- function(workflow_name, ref = GITHUB_REF, repo = GITHUB_REPO, use_cache = FALSE) {
  cache_dir <- file.path(wb_workbench_root(), ".wdl_cache", ref)
  cache_path <- file.path(cache_dir, paste0(workflow_name, ".wdl"))
  if (use_cache && file.exists(cache_path)) {
    wb_msg("INFO", sprintf(
      "Using cached %s.wdl from a previous fetch (%s). Delete this file, or pass use_cache=FALSE, to force a fresh fetch.",
      workflow_name, cache_path
    ))
    return(paste(readLines(cache_path, warn = FALSE), collapse = "\n"))
  }

  path <- wb_workflow_wdl_path(workflow_name)
  text <- paste(wb_gh_fetch(path, ref = ref, raw = TRUE, repo = repo), collapse = "\n")

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(text, cache_path)
  text
}

# A simpler, generic cached fetch for the recursive nested-input scan below (wb_gh_fetch_binary()
# and wb_fetch_workflow_wdl() above have their own, path-shape-specific caching) -- keyed by the
# full repo-relative path (sanitized for use as a filename) since sub-fetches range over
# arbitrarily-nested workflow/task WDLs, not just hydrant/workflows/<name>/<name>.wdl. Fresh by
# default, same reasoning as wb_fetch_workflow_wdl(): branches can change underfoot.
wb_fetch_wdl_cached <- function(path, ref = GITHUB_REF, repo = GITHUB_REPO, use_cache = FALSE) {
  cache_dir <- file.path(wb_workbench_root(), ".wdl_cache", ref, "nested")
  cache_path <- file.path(cache_dir, gsub("[^A-Za-z0-9_.-]", "_", path))
  if (use_cache && file.exists(cache_path)) return(paste(readLines(cache_path, warn = FALSE), collapse = "\n"))
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
### Nested (non-top-level) required-input discovery: a workflow's own input{} block only
### covers what IT declares -- if it calls a sub-workflow without binding one of that
### sub-workflow's own required inputs, Cromwell still needs a value for it, addressable via a
### dot-scoped key (e.g. "panoply_unified_workflow.nmf.run_ssgsea": panoply_nmf_workflow's own
### run_ssgsea toggle, never re-declared by panoply_unified_workflow itself). A flat scan of
### the top workflow's input{} block alone misses these entirely.
### ===

# import "REL_PATH" as ALIAS  ->  c(ALIAS = "REL_PATH", ...)
wb_parse_wdl_imports <- function(wdl_text) {
  lines <- strsplit(gsub("\r", "", wdl_text), "\n")[[1]]
  m <- regmatches(lines, regexec('^\\s*import\\s+"([^"]+)"\\s+as\\s+([A-Za-z_][A-Za-z0-9_]*)', lines))
  hits <- m[lengths(m) == 3]
  if (length(hits) == 0) return(character(0))
  setNames(vapply(hits, `[[`, character(1), 2), vapply(hits, `[[`, character(1), 3))
}

# WDL import paths are relative to the importing file's own directory, same as a filesystem
# path -- resolves e.g. "../panoply_nmf_workflow/panoply_nmf_workflow.wdl" against
# "hydrant/workflows/panoply_unified_workflow/panoply_unified_workflow.wdl" to
# "hydrant/workflows/panoply_nmf_workflow/panoply_nmf_workflow.wdl", collapsing ".." components.
wb_resolve_wdl_import_path <- function(importer_path, import_path) {
  parts <- strsplit(file.path(dirname(importer_path), import_path), "/")[[1]]
  stack <- character(0)
  for (p in parts) {
    if (p == "." || p == "") next
    if (p == "..") stack <- utils::head(stack, -1) else stack <- c(stack, p)
  }
  paste(stack, collapse = "/")
}

# Every `call ALIAS.CALLEE (as CALL_ALIAS)? { input: ... }` block within a workflow's own body,
# with: the import alias, callee name, effective call alias (defaults to the callee name when
# there's no explicit "as"), the set of parameter names bound at this call site, and whether
# the call sits inside a scatter block (a scattered call's inputs can't be individually
# addressed via inputs.json -- there's no single override that applies to every iteration, so
# these are flagged for the caller to skip).
wb_parse_wdl_calls <- function(wdl_text, workflow_name) {
  lines <- strsplit(gsub("\r", "", wdl_text), "\n")[[1]]
  wf_start <- grep(sprintf("^\\s*workflow\\s+%s\\s*\\{", workflow_name), lines)
  if (length(wf_start) == 0) return(list())
  wf_start <- wf_start[1]

  depth <- 0
  scatter_depth <- 0  # depth at which the innermost open scatter(...) began; 0 = not in one
  calls <- list()
  n <- length(lines)

  for (i in seq(wf_start, n)) {
    line <- lines[i]
    trimmed <- trimws(line)
    opens  <- lengths(regmatches(line, gregexpr("\\{", line)))
    closes <- lengths(regmatches(line, gregexpr("\\}", line)))

    if (i > wf_start && grepl("^scatter\\s*\\(", trimmed) && scatter_depth == 0) scatter_depth <- depth + 1

    if (i > wf_start) {
      call_m <- regmatches(trimmed, regexec(
        "^call\\s+([A-Za-z_][A-Za-z0-9_]*)\\.([A-Za-z_][A-Za-z0-9_]*)(?:\\s+as\\s+([A-Za-z_][A-Za-z0-9_]*))?",
        trimmed
      ))[[1]]
      if (length(call_m) == 4 && nzchar(call_m[1])) {
        call_alias <- if (nzchar(call_m[4])) call_m[4] else call_m[3]
        # Scan forward from this line to the call block's own matching closing brace,
        # collecting bound parameter names from its "input:" section along the way.
        call_depth <- opens - closes
        bound <- character(0)
        j <- i
        while (call_depth > 0 && j < n) {
          j <- j + 1
          jline <- lines[j]
          pm <- regmatches(trimws(jline), regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*=", trimws(jline)))[[1]]
          if (length(pm) == 2) bound <- c(bound, pm[2])
          call_depth <- call_depth + lengths(regmatches(jline, gregexpr("\\{", jline))) -
            lengths(regmatches(jline, gregexpr("\\}", jline)))
        }
        calls[[length(calls) + 1]] <- list(
          import_alias = call_m[2], callee = call_m[3], call_alias = call_alias,
          bound = bound, in_scatter = scatter_depth > 0
        )
      }
    }

    depth <- depth + opens - closes
    if (scatter_depth > 0 && depth < scatter_depth) scatter_depth <- 0
    if (i > wf_start && depth <= 0) break
  }
  calls
}

# Recursively walks workflow-to-workflow calls (leaf tasks are out of scope -- checking every
# task's own inputs would multiply the fetch count for comparatively little value, since tasks
# are typically fully parameterized by their immediate caller) looking for inputs that are
# required in some CALLED workflow but never bound at the call site. Returns a list of specs:
# list(dot_key=, param_name=, wdl_type=, optional=, has_default=, is_file=, context=) where
# dot_key is the full inputs.json key (already including the top workflow's own name) and
# context is a human-readable "top -> alias (callee)" path for prompts/logs.
wb_discover_nested_specs <- function(workflow_name, wdl_path, wdl_text, ref, repo,
                                      key_prefix = workflow_name, context = workflow_name,
                                      visited = character(0), max_fetches = 60, fetched = 0,
                                      use_cache = FALSE, notify_scatter = TRUE) {
  # One blanket notice per top-level call (not per scattered call encountered, which could
  # otherwise repeat many times across a large tree) -- scatter-called sub-workflows generally
  # shouldn't have their own required inputs beyond what the top workflow already wires anyway,
  # and a scattered call's inputs can't be individually addressed via inputs.json regardless.
  # notify_scatter = FALSE for wb_build_inputs_json()'s call -- it re-runs this same discovery
  # purely to re-derive dot-key structure for values state$toggles/state$typemap already has
  # (there's no other way to know a nested key's shape without parsing the call graph again),
  # not to check or prompt for anything, so the notice would just be repeated noise there.
  if (notify_scatter && length(visited) == 0) {
    wb_msg("INFO", "Sub-workflows called inside a scatter block are not scanned for their own required inputs.")
  }
  if (workflow_name %in% visited) return(list())  # import-cycle guard
  visited <- c(visited, workflow_name)

  imports <- wb_parse_wdl_imports(wdl_text)
  calls <- wb_parse_wdl_calls(wdl_text, workflow_name)
  specs <- list()

  for (call in calls) {
    if (call$in_scatter) next
    # imports is a named character vector, not a list -- [[ ]] on a missing name throws
    # ("subscript out of bounds") rather than returning NULL, so index with [ ] instead.
    import_path <- unname(imports[call$import_alias])
    if (is.na(import_path)) next
    if (fetched >= max_fetches) {
      wb_msg("WARNING", "Reached the WDL-fetch limit while scanning for nested required inputs; some may be missed.")
      break
    }
    callee_path <- wb_resolve_wdl_import_path(wdl_path, import_path)
    callee_text <- tryCatch(wb_fetch_wdl_cached(callee_path, ref, repo, use_cache), error = function(e) NULL)
    fetched <- fetched + 1
    if (is.null(callee_text)) next

    is_wf <- any(grepl(sprintf("^\\s*workflow\\s+%s\\s*\\{", call$callee), strsplit(callee_text, "\n")[[1]]))
    if (!is_wf) next  # a task -- out of scope, see header comment

    callee_specs <- wb_parse_wdl_inputs(callee_text, call$callee)
    dot_prefix <- paste0(key_prefix, ".", call$call_alias)
    this_context <- sprintf("%s -> %s (%s)", context, call$call_alias, call$callee)

    for (i in seq_len(nrow(callee_specs))) {
      s <- callee_specs[i, ]
      if (s$name %in% call$bound) next  # already satisfied at the call site
      specs[[length(specs) + 1]] <- list(
        dot_key = paste0(dot_prefix, ".", s$name), param_name = s$name, wdl_type = s$wdl_type,
        optional = s$optional, has_default = s$has_default, is_file = s$is_file, context = this_context
      )
    }

    specs <- c(specs, wb_discover_nested_specs(
      call$callee, callee_path, callee_text, ref, repo,
      key_prefix = dot_prefix, context = this_context, visited = visited,
      max_fetches = max_fetches, fetched = fetched, use_cache = use_cache, notify_scatter = notify_scatter
    ))
  }

  specs
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
  geneset_db = "geneset_db", ptm_db = "ptm_db",
  fasta_ref = "FASTA_ref_file"
)

ROLE_TO_SUBSET_CATEGORY <- c(
  prote_ome = "proteome", phospho_ome = "phosphoproteome", acetyl_ome = "acetylome",
  ubiquityl_ome = "ubiquitylome", nglyco_ome = "nglycoproteome", methyl_ome = "methylation",
  metabol_ome = "metabolome", rna_data = "rna", cna_data = "cna",
  groups_file = "groups", groups_file_clumpsptm = "groups_clumpsptm"
)
# fasta_ref (panoply_clumps_ptm_workflow's own FASTA_ref_file) isn't a top-level input of
# panoply_unified_workflow at all -- it's only ever reached via wb_discover_nested_specs()'s
# recursive scan, using the same role-mapping machinery as top-level static files.
ROLE_TO_STATIC_CATEGORY <- c(geneset_db = "gseaDB", ptm_db = "ptmseaDB", fasta_ref = "clumpsFASTA")

TOGGLE_ALIASES <- c(
  run_cmap = "run_cmap", run_mo_nmf = "run_mo_nmf", run_so_nmf = "run_so_nmf",
  run_ptmsea = "run_ptmsea", run_clumps = "run_clumpsptm", run_metab = "run_metab",
  normalizeProteomics = "normalize_proteomics", filterProteomics = "filter_proteomics"
)

# Toggle fields already driven by an earlier, more specific (and data-aware) step elsewhere in
# the pipeline -- wb_select_workflow_toggles() below skips these rather than re-asking a plainer
# version of a question that's already been answered conditionally (e.g. PTM-SEA is only ever
# asked about if phosphoproteome data was actually mapped; asking again here, unconditionally,
# for every workflow would be both redundant and nonsensical when that data isn't present).
TOGGLES_HANDLED_ELSEWHERE <- c("normalize_proteomics", "filter_proteomics", "run_ptmsea", "run_metab", "run_clumpsptm")

wb_map_semantic_role <- function(spec_name) {
  for (role in names(INPUT_ALIASES)) if (spec_name %in% INPUT_ALIASES[[role]]) return(role)
  NA_character_
}

# Prompts for whichever required toggles the CURRENT target workflow actually declares -- both
# its own top-level ones (via TOGGLE_ALIASES) and ones buried in a called sub-workflow that
# were never re-declared at the top (e.g. panoply_nmf_workflow's own run_ssgsea, unbound in
# panoply_unified_workflow's call to it -- see wb_discover_nested_specs()) -- skipping anything
# already handled by one of the more specific steps above. Keeps this generic across workflows
# instead of hardcoding a fixed set of toggle names that only make sense for one of them.
wb_select_workflow_toggles <- function(state, workflow_name = state$target_workflow %||% TARGET_WORKFLOW,
                                        github_ref = state$github_ref %||% GITHUB_REF) {
  wdl_path <- wb_workflow_wdl_path(workflow_name)
  wdl_text <- wb_fetch_workflow_wdl(workflow_name, github_ref)
  specs <- wb_parse_wdl_inputs(wdl_text, workflow_name)

  ask_toggle <- function(label, store_key) {
    current <- state$toggles[[store_key]]
    hint <- if (!is.null(current)) sprintf(" (currently %s)", toupper(as.character(current))) else ""
    state$toggles[[store_key]] <<- wb_confirm(sprintf("Run %s%s?", label, hint))
  }

  asked <- FALSE

  # Top-level: required (no default), name matches a known toggle -- works for Boolean- or
  # String-typed WDL toggles alike (e.g. run_cmap is declared String, not Boolean, in
  # panoply_unified_workflow.wdl; filtering by WDL type alone silently missed it before).
  required <- specs[!specs$optional & !specs$has_default & specs$name %in% names(TOGGLE_ALIASES), , drop = FALSE]
  for (i in seq_len(nrow(required))) {
    wdl_name <- required$name[i]
    state_field <- unname(TOGGLE_ALIASES[wdl_name])
    if (state_field %in% TOGGLES_HANDLED_ELSEWHERE) next
    asked <- TRUE
    ask_toggle(wdl_name, state_field)
  }

  # Nested: required Boolean inputs from a called sub-workflow, unbound at the call site.
  # Anything required but not Boolean-typed (or a File -- handled by wb_build_inputs_json()'s
  # file-role mapping instead) is reported rather than guessed at, since there's no generically
  # safe way to collect an arbitrary String/Int value interactively.
  wb_msg("INFO", "Scanning called sub-workflows for their own required inputs -- this may take a moment...")
  nested <- tryCatch(
    wb_discover_nested_specs(workflow_name, wdl_path, wdl_text, github_ref, GITHUB_REPO),
    error = function(e) {
      wb_msg("WARNING", sprintf("Could not fully scan nested sub-workflow inputs (%s); some may be missed.", conditionMessage(e)))
      list()
    }
  )
  for (s in Filter(function(s) !isTRUE(s$optional) && !isTRUE(s$has_default), nested)) {
    if (isTRUE(s$is_file)) next
    toggle_name <- unname(TOGGLE_ALIASES[s$param_name])
    store_key <- if (!is.na(toggle_name)) toggle_name else s$dot_key
    if (!is.na(toggle_name) && toggle_name %in% TOGGLES_HANDLED_ELSEWHERE) next
    if (!grepl("^Boolean", s$wdl_type)) {
      wb_msg("WARNING", sprintf(
        "Required input '%s' (%s) in %s has no default and isn't auto-promptable -- set it manually in inputs.json (key: %s).",
        s$param_name, s$wdl_type, s$context, s$dot_key
      ))
      next
    }
    asked <- TRUE
    ask_toggle(sprintf("%s [%s]", s$param_name, s$context), store_key)
  }

  if (!asked) wb_msg("INFO", sprintf("No additional required toggles found for '%s'.", workflow_name))

  wb_save_state(state)
}

### ===
### inputs.json build + surgical update
### ===

# state$typemap/state$subsets[[*]]$dir always point at paths under current-session/ -- that's
# where wb_load_and_map_inputs()/wb_create_subset() write, so ongoing edits keep working there
# regardless of whether a session's been named yet. But current-session/ is explicitly the
# *unprotected*, still-mutable one -- baking its paths into a submitted job's inputs.json would
# defeat the entire point of naming and saving a session (a later edit, or clearing the session,
# could invalidate a path the job already depends on). wb_save_session() copies the whole tree,
# so the named session has the same file at the same relative path -- a prefix swap gives the
# right one.
wb_in_named_session <- function(state, local_path) {
  if (is.null(local_path) || is.na(local_path)) return(local_path)
  current_dir <- wb_session_dir()
  named_dir <- wb_session_dir(state$active_named_session)
  if (startsWith(local_path, current_dir)) {
    return(paste0(named_dir, substring(local_path, nchar(current_dir) + 1)))
  }
  local_path
}

# Returns list(inputs=, always_keys=, toggle_keys=) -- see wb_update_inputs_json_for_subset()'s
# surgical-update path for how these two key sets get treated differently on a regeneration.
wb_build_inputs_json <- function(state, subset_name, workflow_name = state$target_workflow %||% TARGET_WORKFLOW,
                                  github_ref = state$github_ref %||% GITHUB_REF) {
  wdl_path <- wb_workflow_wdl_path(workflow_name)
  wdl_text <- wb_fetch_workflow_wdl(workflow_name, github_ref)
  specs <- wb_parse_wdl_inputs(wdl_text, workflow_name)
  subset_files <- wb_subset_files(state, subset_name)  # already resolved against the named session
  master_params_path <- wb_in_named_session(state, file.path(wb_session_dir(), "master-parameters.yaml"))

  resolve_file_role <- function(role) {
    if (is.na(role)) return(NULL)
    if (role == "yaml") return(master_params_path)
    if (role %in% names(ROLE_TO_SUBSET_CATEGORY)) return(subset_files[[ROLE_TO_SUBSET_CATEGORY[[role]]]])
    if (role %in% names(ROLE_TO_STATIC_CATEGORY)) return(wb_in_named_session(state, state$typemap[[ROLE_TO_STATIC_CATEGORY[[role]]]]))
    NULL
  }

  inputs <- list()
  # always_keys: recomputed fresh every time from the subset/session (file paths, job_id) --
  # always safe, and correct, to overwrite on a regeneration, since the whole point of
  # regenerating is to get the right ones for the (possibly new) subset.
  # toggle_keys: user-preference values from state$toggles -- these *could* have been
  # hand-edited directly in an existing inputs.json since the last build, so
  # wb_update_inputs_json_for_subset() only overwrites them if the user opts in.
  always_keys <- character(0)
  toggle_keys <- character(0)

  set_file_input <- function(key, param_name, optional, has_default, role) {
    local_path <- resolve_file_role(role)
    if (!is.null(local_path) && !is.na(local_path) && file.exists(local_path)) {
      inputs[[key]] <<- wb_local_to_s3(local_path)
      always_keys <<- c(always_keys, key)
    } else if (!optional && !has_default) {
      # Genuinely required (no WDL default to fall back on) and we have nothing to offer --
      # a spec with its own default (e.g. panoply_clumps_ptm_workflow's PDB_manifest/
      # UNIPROT_SWISSPROT/SIFTS_DB, each declared "File X = gs://...") is left alone entirely:
      # no warning, and critically not added to always_keys, so a surgical update never nulls
      # out whatever's already there (Cromwell's own default applies if nothing's set at all).
      wb_msg("WARNING", sprintf("Required File input '%s' could not be resolved -- fill it in manually.", param_name))
    }
  }

  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    role <- wb_map_semantic_role(spec$name)
    key <- paste0(workflow_name, ".", spec$name)

    if (isTRUE(spec$is_file)) {
      set_file_input(key, spec$name, spec$optional, spec$has_default, role)
      next
    }

    if (!is.na(role) && role == "job_id") {
      # Combine the named session with the subset, e.g. "odg-v4-ODG" -- not just the subset
      # name, which alone can't distinguish output labeled with the same subset name across
      # different named sessions/runs. state$job_id, if explicitly set, still overrides this.
      inputs[[key]] <- as.character(state$job_id %||% paste0(state$active_named_session, "-", subset_name))
      always_keys <- c(always_keys, key)
      next
    }

    toggle_name <- unname(TOGGLE_ALIASES[spec$name])
    if (!is.na(toggle_name)) {
      value <- state$toggles[[toggle_name]]
      if (!is.null(value)) {
        inputs[[key]] <- if (grepl("^Boolean", spec$wdl_type)) as.logical(value) else tolower(as.character(value))
        toggle_keys <- c(toggle_keys, key)
      }
    }
  }

  # Nested (non-top-level) specs -- unbound at their call site inside some sub-workflow, so
  # never appear in the top workflow's own input{} block (e.g. panoply_clumps_ptm_workflow's
  # own FASTA_ref_file, or panoply_nmf_workflow's own run_ssgsea). See wb_discover_nested_specs().
  # notify_scatter = FALSE -- this is re-deriving dot-key structure to wire values already
  # decided by wb_select_workflow_toggles(), not checking or prompting for anything, so the
  # scatter-skip notice (genuinely useful there) would just be repeated noise here.
  nested <- tryCatch(
    wb_discover_nested_specs(workflow_name, wdl_path, wdl_text, github_ref, GITHUB_REPO, notify_scatter = FALSE),
    error = function(e) {
      wb_msg("WARNING", sprintf("Could not fully scan nested sub-workflow inputs (%s); some may be missed.", conditionMessage(e)))
      list()
    }
  )
  for (s in nested) {
    role <- wb_map_semantic_role(s$param_name)
    if (isTRUE(s$is_file)) {
      set_file_input(s$dot_key, s$param_name, s$optional, s$has_default, role)
      next
    }
    toggle_name <- unname(TOGGLE_ALIASES[s$param_name])
    value <- if (!is.na(toggle_name)) state$toggles[[toggle_name]] else state$toggles[[s$dot_key]]
    if (!is.null(value)) {
      inputs[[s$dot_key]] <- if (grepl("^Boolean", s$wdl_type)) as.logical(value) else tolower(as.character(value))
      toggle_keys <- c(toggle_keys, s$dot_key)
    }
  }

  list(inputs = inputs, always_keys = always_keys, toggle_keys = toggle_keys)
}

wb_update_inputs_json_for_subset <- function(state, subset_name = NULL,
                                              workflow_name = state$target_workflow %||% TARGET_WORKFLOW,
                                              github_ref = state$github_ref %||% GITHUB_REF,
                                              existing_inputs_path = NULL,
                                              out_path = NULL) {
  # Requires a named session (not current-session) so this file's S3 paths -- once baked
  # into a submitted job -- can't be silently invalidated by later, unrelated work in
  # current-session. See wb_save_session() in sessions.r.
  if (is.null(state$active_named_session)) {
    stop("No named session found -- run `state <- wb_save_session(state)` first ",
         "(see the Sessions section) before generating inputs.json.")
  }

  if (is.null(subset_name)) {
    subset_names <- wb_list_subsets(state)
    if (length(subset_names) == 0) stop("No subsets found -- run wb_create_subset() first.")
    cat("Subsets:\n")
    for (i in seq_along(subset_names)) cat(sprintf("  %2d: %s\n", i, subset_names[i]))
    flush.console()
    idx <- wb_smart_readline(
      "Select a subset to generate inputs.json for (or 'quit' to cancel): ",
      valid = function(ch) {
        n <- suppressWarnings(as.integer(ch))
        if (is.na(n) || n < 1 || n > length(subset_names)) sprintf("Please enter a number from 1 to %d.", length(subset_names)) else TRUE
      }
    )
    if (is.null(idx)) {
      wb_msg("CANCELLED", "inputs.json not generated.")
      wb_done()
      return(invisible(NULL))
    }
    subset_name <- subset_names[as.integer(idx)]
  }

  # existing_inputs_path still defaults to the standard per-session path -- but if something's
  # already there, offer the option to point at a different copy instead (e.g. one you moved,
  # renamed, or have been hand-editing) rather than always updating that default in place.
  default_existing <- file.path(wb_session_dir(state$active_named_session), "inputs.json")
  if (is.null(existing_inputs_path)) {
    existing_inputs_path <- default_existing
    if (file.exists(default_existing) &&
        !wb_confirm(sprintf("Update the existing inputs.json at '%s'?", default_existing))) {
      alt <- wb_smart_readline(
        "Path to the inputs.json you'd like to update instead (local or this project's s3://, or 'quit' to use the default): ",
        valid = function(ch) wb_validate_user_file(ch, extensions = "json")
      )
      if (!is.null(alt)) existing_inputs_path <- wb_resolve_user_path(alt)
    }
  }
  out_path <- out_path %||% existing_inputs_path
  built <- wb_build_inputs_json(state, subset_name, workflow_name, github_ref)
  fresh <- built$inputs

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

  # File paths and job_id always get refreshed -- they're recomputed from the (possibly new)
  # subset/session every time, so there's nothing to preserve. Toggles are different: there's
  # no reliable way to tell whether you've re-run wb_select_workflow_toggles() with a genuine
  # change since this file was last built, versus this being the only signal that you've since
  # hand-edited one of those same keys directly in the JSON -- so ask, rather than guess.
  keys_to_refresh <- built$always_keys
  if (length(built$toggle_keys) > 0 &&
      wb_confirm(paste(
        "Refresh required toggle/parameter values with your current settings?",
        "(This will only impact parameters set in the previous cell; ",
        "other parameters in inputs.json will be left as-is.)"
      ))) {
    keys_to_refresh <- c(keys_to_refresh, built$toggle_keys)
  }
  for (key in keys_to_refresh) existing[[key]] <- fresh[[key]]

  jsonlite::write_json(existing, out_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  wb_msg("INFO", sprintf("Updated %d input(s) in %s for subset '%s' (backup at %s.bak)",
                        length(keys_to_refresh), out_path, subset_name, existing_inputs_path))
  wb_done()
  out_path
}
