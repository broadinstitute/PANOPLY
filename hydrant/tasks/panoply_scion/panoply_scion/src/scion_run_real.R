#!/usr/bin/env Rscript
# Runs the real (unpermuted) SCION network and saves everything a downstream
# permutation shard needs: the processed target/regulator matrices, the fixed
# cluster assignment, and the exact params used (so permutations stay
# consistent with the real network -- see compute_fdr_threshold()).

suppressMessages({
  library(optparse)
  library(SCION)
})

option_list <- list(
  make_option("--target_data_file", type = "character"),
  make_option("--reg_data_file", type = "character"),
  make_option("--target_genes_file", type = "character", default = NULL),
  make_option("--reg_genes_file", type = "character", default = NULL),
  make_option("--gene_list_header", type = "logical", default = TRUE),
  make_option("--clustering_data_file", type = "character", default = NULL),
  make_option("--format", type = "character", default = "csv"),
  make_option("--clustering_method", type = "character", default = "none"),
  make_option("--clustering_threshold", type = "double", default = 0.5),
  make_option("--clusters_file", type = "character", default = NULL),
  make_option("--connect_hubs", type = "logical", default = TRUE),
  make_option("--weightthreshold", type = "double", default = 0),
  make_option("--normalize", type = "logical", default = TRUE),
  make_option("--num_cores", type = "integer", default = 1),
  make_option("--engine", type = "character", default = "randomForest"),
  make_option("--ptm_sep", type = "character", default = "."),
  make_option("--seed", type = "integer", default = 2020),
  make_option("--nb_trees", type = "integer", default = 10000),
  make_option("--out_dir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

result <- run_scion(
  target_data_file = opt$target_data_file,
  reg_data_file = opt$reg_data_file,
  target_genes_file = opt$target_genes_file,
  reg_genes_file = opt$reg_genes_file,
  gene_list_header = opt$gene_list_header,
  clustering_data_file = opt$clustering_data_file,
  format = opt$format,
  clustering_method = opt$clustering_method,
  clustering_threshold = opt$clustering_threshold,
  clusters_file = opt$clusters_file,
  connect_hubs = opt$connect_hubs,
  weightthreshold = opt$weightthreshold,
  normalize = opt$normalize,
  num.cores = opt$num_cores,
  engine = opt$engine,
  ptm_sep = opt$ptm_sep,
  seed = opt$seed,
  permute = FALSE,
  nb.trees = opt$nb_trees
)

write_scion_network(result$network, file.path(opt$out_dir, "network.tsv"))
saveRDS(result$network, file.path(opt$out_dir, "network.rds"))
saveRDS(result$target, file.path(opt$out_dir, "target.rds"))
saveRDS(result$reg, file.path(opt$out_dir, "reg.rds"))
saveRDS(result$cluster_assignment, file.path(opt$out_dir, "cluster_assignment.rds"))
saveRDS(c(result$params, list(nb.trees = opt$nb_trees)), file.path(opt$out_dir, "params.rds"))

cat("Done. Real network:", nrow(result$network), "edges.\n")
