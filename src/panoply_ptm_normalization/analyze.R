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


summarize_models <- function(data_dir, models = NULL) {
  out_dir <- file.path(data_dir, "out")
  if (is.null(models)) {
    models <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE)
  }

  # placeholder for summary dataframes
  corr_acr_samples_summary <- NULL
  corr_in_samples_summary <- NULL
  
  for (model in models) {
    model_path <- file.path(out_dir, model)
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
  }
  
  rownames(corr_acr_samples_summary) <- NULL
  rownames(corr_in_samples_summary) <- NULL
  write.csv(corr_acr_samples_summary, file.path(out_dir, "corr_acr_samples_summary.csv"))
  write.csv(corr_in_samples_summary, file.path(out_dir, "corr_in_samples_summary.csv"))
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
  
  distr_acr <- ggplot(metric_acr_samples, aes(x = pearson_r, color = model)) +
    geom_histogram(fill = "white", alpha = 0, position = "identity", bins = 50) +
    geom_vline(data = medians, aes(xintercept = median, color = model), size = 0.5, alpha = 0.8)
  dist_acr_log <- distr_acr + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2))
  comb_fig <- ggarrange(distr_acr, dist_acr_log, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right")
  
  return(comb_fig)
}