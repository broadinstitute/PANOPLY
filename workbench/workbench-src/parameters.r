# master-parameters.yaml build: canonical defaults (fetched live, with a bundled fallback)
# merged with notebook-collected runtime overrides.

wb_load_default_master_parameters <- function(github_ref = GITHUB_REF) {
  local_fallback <- file.path("workbench-src", "defaults", "master-parameters.yaml")
  tryCatch({
    text <- wb_gh_fetch("src/panoply_common/master-parameters.yaml", ref = github_ref, raw = TRUE)
    yaml::yaml.load(paste(text, collapse = "\n"))
  }, error = function(e) {
    wb_msg("WARNING", sprintf(
      "Could not fetch master-parameters.yaml from GitHub (%s); using bundled fallback copy.",
      conditionMessage(e)
    ))
    yaml::read_yaml(local_fallback)
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
  out_path
}
