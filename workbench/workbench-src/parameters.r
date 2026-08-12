# master-parameters.yaml build: the stable, staged default (pushed at deploy time), with a
# live-fetch-from-GitHub backup, merged with notebook-collected runtime overrides.

wb_load_default_master_parameters <- function(github_ref = GITHUB_REF) {
  # Staged copy first -- populated by deploy-workbench.sh directly from the canonical
  # src/panoply_common/ right before upload, so this is the stable, version-pinned copy that
  # shipped with this deployment. Dev-mode fallback for running straight out of a full repo
  # checkout. Only if neither is present do we live-fetch from GitHub as a backup -- with a
  # warning, since a live fetch may not match the version actually pinned at deployment.
  staged_candidates <- c(
    file.path("workbench-src", "defaults", "master-parameters.yaml"),
    file.path("..", "src", "panoply_common", "master-parameters.yaml")
  )
  staged <- staged_candidates[file.exists(staged_candidates)][1]
  if (!is.na(staged)) return(yaml::read_yaml(staged))

  wb_msg("WARNING", sprintf(
    paste("No staged master-parameters.yaml found (expected at %s); live-fetching from",
          "GitHub (%s@%s) as a backup -- this may not match the stable, version-pinned copy",
          "normally staged when the workbench is deployed."),
    paste(staged_candidates, collapse = " or "), GITHUB_REPO, github_ref
  ))
  tryCatch({
    text <- wb_gh_fetch("src/panoply_common/master-parameters.yaml", ref = github_ref, raw = TRUE)
    yaml::yaml.load(paste(text, collapse = "\n"))
  }, error = function(e) {
    stop("No staged master-parameters.yaml found, and the live GitHub fetch failed too (",
         conditionMessage(e), ").")
  })
}

wb_build_master_parameters_yaml <- function(state,
                                             out_path = file.path(wb_workbench_root(), "master-parameters.yaml"),
                                             github_ref = state$github_ref) {
  defaults <- if (!is.null(state$typemap$parameters)) {
    yaml::read_yaml(state$typemap$parameters)
  } else {
    wb_load_default_master_parameters(github_ref)
  }

  overrides <- list()
  overrides[["normalize.proteomics"]]   <- state$toggles$normalize_proteomics
  overrides[["filter.proteomics"]]      <- state$toggles$filter_proteomics
  overrides[["run.ptmsea"]]             <- state$toggles$run_ptmsea
  overrides[["run.metab"]]              <- state$toggles$run_metab
  overrides[["run.clumpsptm"]]          <- state$toggles$run_clumpsptm
  overrides[["cosmo.params"]]           <- state$cosmo_params
  overrides[["groups.cols"]]            <- state$groups_cols
  overrides[["groups.cols.continuous"]] <- state$groups_cols_continuous
  overrides[["groups.colors"]]          <- state$groups_colors

  merged <- modifyList(defaults, overrides)

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  yaml::write_yaml(merged, out_path, handlers = list(logical = function(x) {
    structure(ifelse(x, "TRUE", "FALSE"), class = "verbatim")
  }))
  wb_done()
  out_path
}
