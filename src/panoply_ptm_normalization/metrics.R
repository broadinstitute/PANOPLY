#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Evaluates performance of PTM normalization methods with various parameter settings
# Khoi Pham (Munchic)
#
# In-sample and across sample metrics to compare PTM/protein correlation in and across samples
#
#########################################################################################################


library("pacman")
p_load(cmapR)
p_load(tidyr)
p_load(dplyr)
p_load(purrr)

source("panoply_ptm_normalization/helper.R")


run_all_metrics <- function(
  prot_gct_path,
  ptm_gct_path,
  ptm_norm_gct_path,
  groups_colname = NULL,  # name of column indicating grouping (e.g., "pert_time")
  subset_cond = NULL,
  min_n_values = 4,
  method = "pearson")
{
  # obtain across samples correlations for prot and PTM
  prot_ptm_unnorm_corr <- prot_ptm_corr_across_samples(prot_gct_path, ptm_gct_path, groups_colname, subset_cond, min_n_values, method)
  colnames(prot_ptm_unnorm_corr) <- add_prefix_to_colnames("prot_ptm_unnorm", prot_ptm_unnorm_corr, except = "id.x")
  
  # obtain across samples correlations for prot and normalized PTM
  prot_ptm_norm_corr <- prot_ptm_corr_across_samples(prot_gct_path, ptm_norm_gct_path, groups_colname, subset_cond, min_n_values, method)
  colnames(prot_ptm_norm_corr) <- add_prefix_to_colnames("prot_ptm_norm", prot_ptm_norm_corr, except = "id.x")
  
  # obtain across samples correlations for PTM and normalized PTM
  ptm_unnorm_ptm_norm_corr <- ptm_corr_across_samples(ptm_gct_path, ptm_norm_gct_path, groups_colname, subset_cond, min_n_values, method)
  colnames(ptm_unnorm_ptm_norm_corr) <- add_prefix_to_colnames("ptm_unnorm_ptm_norm", ptm_unnorm_ptm_norm_corr, except = "id.x")
  
  comb_df <- merge(prot_ptm_unnorm_corr, prot_ptm_norm_corr, by = "id.x")
  comb_df <- merge(comb_df, ptm_unnorm_ptm_norm_corr, by = "id.x")
  
  return(comb_df)
}

prot_ptm_corr_across_samples <- function(
  prot_gct_path,
  ptm_gct_path,
  groups_colname = NULL,  # name of column indicating grouping (e.g., "pert_time")
  subset_cond = NULL,
  min_n_values = 4,
  method = "pearson")
{
  prot <- parse_gctx(prot_gct_path)
  temp_ptm <- parse_gctx(ptm_gct_path)
  ptm <- match_ptm_to_proteome(temp_ptm, prot)
  
  comb_ptm_prot <- merge_ptm_prot_df(ptm, prot)
  if (!is.null(subset_groups)) {
    comb_ptm_prot <- comb_ptm_prot %>% filter(!!rlang::parse_expr(subset_cond))
  }
  corr_store <- list()
  
  # get correlations for all PTMs
  corr_store[["all"]] = list()
  for (ptm_id in unique(comb_ptm_prot$id.x)) {
    acr_samples <- comb_ptm_prot %>% filter(id.x == ptm_id)
    corr_store[["all"]][[ptm_id]] <- run_corr(acr_samples$value.prot, acr_samples$value, min_n_values, method = method)
  }
  
  # get correlations for PTMs by grouping
  if (!is.null(groups_colname)) {
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
    
    for (group in sample_groups) {
      group_name <- paste0(groups_colname, ":", group)
      corr_store[[group_name]] = list()
      
      for (ptm_id in unique(comb_ptm_prot$id.x)) {
        acr_samples <- comb_ptm_prot %>% filter(id.x == ptm_id) %>% filter(!!rlang::sym(groups_colname) == group)
        corr_store[[group_name]][[ptm_id]] <- run_corr(acr_samples$value.prot, acr_samples$value, min_n_values, method = method)
      }
    }
  }
  
  corr_df <- corr_df_from_list(corr_store)
  return(corr_df)
}

ptm_corr_across_samples <- function(ptm_gct_path, ptm_normalized_gct_path, groups_colname = NULL, subset_cond = NULL, min_n_values = 4, method = "pearson") {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_norm <- parse_gctx(ptm_normalized_gct_path)
  ptm_vals <- melt_gct(ptm)[, c("id.y", "id.x", "value")]
  ptm_norm_vals <- melt_gct(ptm_norm)[, c("id.y", "id.x", "value", ..groups_colname)]
  
  comb_ptm <- merge(ptm_vals, ptm_norm_vals, by = c("id.x", "id.y"))
  colnames(comb_ptm) <- c("id.x", "id.y", "value", "value.norm", groups_colname)  # TODO: this hard-coded order is quite dangerous...
  
  if (!is.null(comb_ptm)) {
    comb_ptm <- comb_ptm %>% filter(!!rlang::parse_expr(subset_cond))
  }
  corr_store <- list()
  
  # get correlations for all PTMs
  corr_store[["all"]] = list()
  for (ptm_id in unique(comb_ptm$id.x)) {
    acr_samples <- comb_ptm %>% filter(id.x == ptm_id)
    corr_store[["all"]][[ptm_id]] <- run_corr(acr_samples$value, acr_samples$value.norm, min_n_values, method = method)
  }
  
  # get correlations for PTMs by grouping
  if (!is.null(groups_colname)) {
    sample_groups <- unique(comb_ptm[[groups_colname]])
    
    for (group in sample_groups) {
      group_name <- paste0(groups_colname, ":", group)
      corr_store[[group_name]] = list()
      
      for (ptm_id in unique(comb_ptm$id.x)) {
        acr_samples <- comb_ptm %>% filter(id.x == ptm_id) %>% filter(!!rlang::sym(groups_colname) == group)
        corr_store[[group_name]][[ptm_id]] <- run_corr(acr_samples$value, acr_samples$value.norm, min_n_values, method = method)
      }
    }
  }
  
  corr_df <- corr_df_from_list(corr_store)
  return(corr_df)
}

prot_ptm_corr_per_sample <- function(prot_gct_path, ptm_gct_path) {
  prot <- parse_gctx(prot_gct_path)
  temp_ptm <- parse_gctx(ptm_gct_path)
  ptm <- match_ptm_to_proteome(temp_ptm, prot)
  
  comb_ptm_prot <- merge_ptm_prot_df(ptm, prot)
  corr_store <- list()
  for (sample_name in unique(comb_ptm_prot$id.y)) {
    sample <- comb_ptm_prot %>% filter(id.y == sample_name)
    corr_store[sample_name] <- cor(sample$value.prot, sample$value, method = "pearson")
  }
  
  return(corr_store)
}

ptm_corr_per_sample <- function(ptm_gct_path, ptm_normalized_gct_path) {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_norm <- parse_gctx(ptm_normalized_gct_path)
  ptm_vals <- melt_gct(ptm)[, c("id.y", "id.x", "value")]
  ptm_norm_vals <- melt_gct(ptm_norm)[, c("id.y", "id.x", "value")]
  
  comb_ptm <- merge(ptm_vals, ptm_norm_vals, by = c("id.x", "id.y"))
  colnames(comb_ptm) <- c("id.x", "id.y", "value", "value.norm")
  
  corr_store <- list()
  for (sample_name in unique(comb_ptm$id.y)) {
    sample <- comb_ptm %>% filter(id.y == sample_name)
    corr_store[sample_name] <- cor(sample$value, sample$value.norm, method = "pearson")
  }
  
  return(corr_store)
}

run_corr <- function(list_1, list_2, min_n_values = 4, method = "pearson") {
  if (sum(!(is.na(list_1) | is.na(list_2))) >= min_n_values) {
    corr <- cor(list_1, list_2, method = method)
  } else {
    corr <- NA
  }
  
  return(corr)
}
