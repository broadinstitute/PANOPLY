#!/usr/bin/env Rscript
# Runs exactly ONE permutation, by index -- meant to be called once per shard
# in a Terra/WDL scatter (or one HPC job-array task). Reuses the fixed
# target/reg/cluster_assignment/params saved by scion_run_real.R so this
# permutation is directly comparable to the real network.

suppressMessages({
  library(optparse)
  library(SCION)
})

option_list <- list(
  make_option("--target_rds", type = "character"),
  make_option("--reg_rds", type = "character"),
  make_option("--cluster_assignment_rds", type = "character"),
  make_option("--params_rds", type = "character"),
  make_option("--index", type = "integer"),
  make_option("--base_seed", type = "integer", default = 0),
  make_option("--permute_dim", type = "character", default = "col"),
  make_option("--num_cores", type = "integer", default = 1),
  make_option("--out_dir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

target <- readRDS(opt$target_rds)
reg <- readRDS(opt$reg_rds)
cluster_assignment <- readRDS(opt$cluster_assignment_rds) # NULL when no clustering was used
params <- readRDS(opt$params_rds)

permuted <- permute_network(
  target, reg,
  cluster_assignment = cluster_assignment,
  indices = opt$index,
  permute_dim = opt$permute_dim,
  base_seed = opt$base_seed,
  num.cores = opt$num_cores,
  weightthreshold = params$weightthreshold,
  normalize = params$normalize,
  connect_hubs = params$connect_hubs,
  ptm_sep = params$ptm_sep,
  nb.trees = params$nb.trees
)

save_permutation(permuted[[1]], opt$index, dir = opt$out_dir)

cat("Done. Permutation", opt$index, "-", nrow(permuted[[1]]), "edges.\n")
