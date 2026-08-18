# Local-folder sample subsetting -- replaces Terra sample sets.

wb_create_subset <- function(state, name, filter_col = NULL, filter_vals = NULL,
                              out_root = file.path(wb_session_dir(), "subsets")) {
  annot <- read.csv(state$typemap$annotation, stringsAsFactors = FALSE, quote = '"')
  sample_ids <- if (is.null(filter_col)) {
    annot$Sample.ID
  } else {
    if (!(filter_col %in% colnames(annot))) stop(sprintf("'%s' is not a column in the annotation table.", filter_col))
    annot$Sample.ID[annot[[filter_col]] %in% filter_vals]
  }
  if (length(sample_ids) == 0) stop("No samples matched the given filter.")

  subset_dir <- file.path(out_root, name)
  dir.create(subset_dir, showWarnings = FALSE, recursive = TRUE)

  gct_categories <- intersect(names(state$typemap), c(PROTEOME_TYPES, "rna", "cna", "metabolome"))
  for (cat_name in gct_categories) {
    gct <- cmapR::parse_gctx(state$typemap[[cat_name]])
    keep <- intersect(sample_ids, gct@cid)
    if (length(keep) == 0) {
      wb_msg("WARNING", sprintf("No samples from this subset found in '%s'; skipping.", cat_name))
      next
    }
    sub <- cmapR::subset_gct(gct, cid = keep)
    wb_write_gct_atomic(sub, file.path(subset_dir, paste0(cat_name, ".gct")))
  }

  csv_categories <- intersect(names(state$typemap), c("annotation", "groups", "groups_clumpsptm"))
  for (cat_name in csv_categories) {
    csv <- read.csv(state$typemap[[cat_name]], stringsAsFactors = FALSE, quote = '"')
    write.csv(csv[csv$Sample.ID %in% sample_ids, , drop = FALSE],
              file.path(subset_dir, paste0(cat_name, ".csv")), row.names = FALSE, quote = TRUE)
  }

  state$subsets[[name]] <- list(filter_col = filter_col, filter_vals = filter_vals,
                                dir = subset_dir, n_samples = length(sample_ids))
  wb_msg("INFO", sprintf("Subset '%s' created with %d sample(s) at %s", name, length(sample_ids), subset_dir))
  wb_save_state(state)
}

wb_list_subsets <- function(state) names(state$subsets)

wb_subset_files <- function(state, name) {
  if (is.null(state$subsets[[name]])) stop(sprintf("No subset named '%s' -- run wb_create_subset() first.", name))
  files <- list.files(state$subsets[[name]]$dir, full.names = TRUE)
  setNames(as.list(files), sub("\\.(gct|csv)$", "", basename(files)))
}
