#!/usr/bin/env Rscript
# Gathers all per-permutation .rds files (from scion_run_permutation.R shards)
# and the real network, computes the FDR threshold, and writes the final
# thresholded network + diagnostic plots.

suppressMessages({
  library(optparse)
  library(SCION)
})

option_list <- list(
  make_option("--network_rds", type = "character"),
  make_option("--network_tsv", type = "character"),
  make_option("--permutation_dir", type = "character"),
  make_option("--target_fdr", type = "double", default = 0.05),
  make_option("--out_dir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

real_network <- readRDS(opt$network_rds)
permuted_networks <- load_permutations(opt$permutation_dir)

fdr_result <- compute_fdr_threshold(real_network, permuted_networks, target_fdr = opt$target_fdr)
saveRDS(fdr_result, file.path(opt$out_dir, "fdr_result.rds"))

write_scion_network(fdr_result$thresholded_network, file.path(opt$out_dir, "thresholded_network.tsv"))

ggplot2::ggsave(file.path(opt$out_dir, "fdr_curve.png"),
                plot_fdr_curve(fdr_result, type = "curve"), width = 6, height = 5)
ggplot2::ggsave(file.path(opt$out_dir, "weight_comparison.png"),
                plot_fdr_curve(fdr_result, type = "weight_comparison"), width = 6, height = 5)

if (nrow(fdr_result$thresholded_network) > 0) {
  grDevices::png(file.path(opt$out_dir, "network_plot.png"), width = 2000, height = 2000, res = 300)
  plot_network(fdr_result$thresholded_network)
  grDevices::dev.off()
}

# bundle everything into one tar, matching how panoply_cmap_analysis/
# blacksheep/immune_analysis/so_nmf each always produce a single tar output
bundle_files <- c("thresholded_network.tsv", "fdr_curve.png", "weight_comparison.png", "fdr_result.rds")
if (file.exists(file.path(opt$out_dir, "network_plot.png"))) {
  bundle_files <- c(bundle_files, "network_plot.png")
}
file.copy(opt$network_tsv, file.path(opt$out_dir, "network.tsv"))
bundle_files <- c(bundle_files, "network.tsv")

starting_dir <- getwd()
on.exit(setwd(starting_dir), add = TRUE)
setwd(opt$out_dir)
utils::tar("scion_results.tar.gz", files = bundle_files, compression = "gzip")
setwd(starting_dir)

cat("Done. Threshold:", fdr_result$threshold,
    "| kept", nrow(fdr_result$thresholded_network), "of", nrow(real_network), "edges.\n")
