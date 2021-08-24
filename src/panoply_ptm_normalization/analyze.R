#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Evaluates PTM normalization methods based on across and in-sample correlations of proteome, PTM, normalized PTM
# Khoi Pham (Munchic)
#
# Summarize metrics into tables and distributions to analyze results.
#
#########################################################################################################

library("pacman")
p_load(cmapR)
p_load(tidyr)
p_load(ggplot2)
library(ggpubr)


summarize_metrics <- function(data_dir, models = NULL) {
  out_dir <- file.path(data_dir, "out")
  if (is.null(models)) {
    models <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE)
  }

  # placeholder for summary dataframes
  corr_acr_samples_summary <- NULL
  corr_in_samples_summary <- NULL
  ptm_stats_summary <- NULL
  ptm_stats_per_group_summary <- NULL
  
  for (model in models) {
    model_path <- file.path(out_dir, model)
    # correlation metrics
    corr_acr_samples <- read.csv(file.path(model_path, "corr_across_samples.csv"))
    corr_in_samples <- read.csv(file.path(model_path, "corr_in_samples.csv"))
    corr_acr_samples <- corr_acr_samples[, -which(names(corr_acr_samples) == "id.x")]
    corr_in_samples <- corr_in_samples[, -which(names(corr_in_samples) == "id.y")]
    
    corr_acr_samples_mean <- sapply(corr_acr_samples, mean, na.rm = T)
    corr_acr_samples_mean <- c(model_stat = paste0(model, ".mean"), corr_acr_samples_mean)
    corr_acr_samples_std <- sapply(corr_acr_samples, sd, na.rm = T)
    corr_acr_samples_std <- c(model_stat = paste0(model, ".std"), corr_acr_samples_std)
    corr_acr_samples_summary <- rbind(corr_acr_samples_summary, corr_acr_samples_mean, corr_acr_samples_std)
    
    corr_in_samples_mean <- sapply(corr_in_samples, mean, na.rm = T)
    corr_in_samples_mean <- c(model_stat = paste0(model, ".mean"), corr_in_samples_mean)
    corr_in_samples_std <- sapply(corr_in_samples, sd, na.rm = T)
    corr_in_samples_std <- c(model_stat = paste0(model, ".std"), corr_in_samples_std)
    corr_in_samples_summary <- rbind(corr_in_samples_summary, corr_in_samples_mean, corr_in_samples_std)
    
    # log fold change metrics
    ptm_stats <- read.csv(file.path(model_path, "ptm_log_fold_stats.csv"))
    ptm_stats[["model_stat"]] <- add_prefix_to_series(model, ptm_stats[["model_stat"]], sep = ".")
    ptm_stats_summary <- rbind(ptm_stats_summary, ptm_stats)
    
    ptm_stats_per_group <- read.csv(file.path(model_path, "ptm_log_fold_stats_per_group.csv"))
    ptm_stats_per_group[["model_stat"]] <- add_prefix_to_series(model, ptm_stats_per_group[["model_stat"]], sep = ".")
    ptm_stats_per_group_summary <- rbind(ptm_stats_per_group_summary, ptm_stats_per_group)
  }
  
  rownames(corr_acr_samples_summary) <- NULL
  rownames(corr_in_samples_summary) <- NULL
  rownames(ptm_stats_summary) <- NULL
  rownames(ptm_stats_per_group_summary) <- NULL
  
  write.csv(corr_acr_samples_summary, file.path(out_dir, "corr_acr_samples_summary.csv"))
  write.csv(corr_in_samples_summary, file.path(out_dir, "corr_in_samples_summary.csv"))
  write.csv(ptm_stats_summary, file.path(out_dir, "ptm_stats_summary.csv"))
  write.csv(ptm_stats_per_group_summary, file.path(out_dir, "ptm_stats_per_group_summary.csv"))
}

metric_distr_across_models <- function(data_dir, metric, models = NULL) {
  out_dir <- file.path(data_dir, "out")
  if (is.null(models)) {
    models <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE)
  }
  
  metric_store_acr_samples <- list()
  means_store <- list()
  for (model in models) {
    model_path <- file.path(out_dir, model)
    corr_acr_samples <- read.csv(file.path(model_path, "corr_across_samples.csv"))
    metric_store_acr_samples[[model]] <- corr_acr_samples[[metric]]
  }
  
  metric_acr_samples <- stack(metric_store_acr_samples)
  colnames(metric_acr_samples) <- c("pearson_r", "model")
  
  medians <- metric_acr_samples %>%
    group_by(model) %>%
    summarize(median = median(pearson_r, na.rm=TRUE))

  distr_acr <- ggplot(metric_acr_samples, aes(x = model, y = pearson_r)) +
               geom_boxplot(outlier.shape = NA) + stat_boxplot(geom = "errorbar") +
               geom_violin(aes(fill = model, color = model), alpha = 0.5) +
               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(file.path(out_dir, paste0("distr_across_models-", metric, ".png")), plot = distr_acr)
}