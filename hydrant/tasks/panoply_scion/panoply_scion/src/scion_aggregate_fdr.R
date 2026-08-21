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

cat("Done. Threshold:", fdr_result$threshold,
    "| kept", nrow(fdr_result$thresholded_network), "of", nrow(real_network), "edges.\n")
