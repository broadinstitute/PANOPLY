#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Evaluates performance of PTM normalization methods with various parameter settings
# Khoi Pham (Munchic)
#
# Various metrics to compare PTM/protein correlation in and across samples
#
#########################################################################################################


library("pacman")
p_load(cmapR)
p_load(tidyr)
p_load(dplyr)
p_load(purrr)

prot_ptm_corr_per_sample <- function(prot_gct_path, ptm_gct_path) {
  prot <- parse_gctx(prot_gct_path)
  temp_ptm <- parse_gctx(ptm_gct_path)
  ptm <- match_ptm_to_proteome(temp_ptm, prot)
  
  combined_data <- merge_ptm_prot_df(ptm, prot)
  corr_store <- list()
  for (sample_name in unique(combined_data$id.y)) {
    sample <- combined_data %>% filter(id.y == sample_name)
    corr_store[sample_name] <- cor(sample$value.prot, sample$value, method = "pearson")
  }
  
  return(corr_store)
}

ptm_corr_per_sample <- function(ptm_gct_path, ptm_normalized_gct_path) {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_norm <- parse_gctx(ptm_normalized_gct_path)
  ptm_vals <- melt_gct(ptm)[, c("id.y", "id.x", "value")]
  ptm_norm_vals <- melt_gct(ptm_norm)[, c("id.y", "id.x", "value")]
  
  combined_ptm <- merge(ptm_vals, ptm_norm_vals, by = c("id.x", "id.y"))
  colnames(combined_ptm) <- c("id.x", "id.y", "value", "value.norm")
  
  corr_store <- list()
  for (sample_name in unique(combined_ptm$id.y)) {
    sample <- combined_ptm %>% filter(id.y == sample_name)
    corr_store[sample_name] <- cor(sample$value, sample$value.norm, method = "pearson")
  }
  
  return(corr_store)
}

prot_ptm_corr_across_samples <- function(prot_gct_path, ptm_gct_path, use_all_samples = FALSE, sample_groups = c("06h", "24h", "96h"), min_values_present = 4) {
  prot <- parse_gctx(prot_gct_path)
  temp_ptm <- parse_gctx(ptm_gct_path)
  ptm <- match_ptm_to_proteome(temp_ptm, prot)
  
  combined_data <- merge_ptm_prot_df(ptm, prot)
  corr_store <- list()
  
  if (use_all_samples) {
    sample_groups = c("")  # just so that it matches to everything
  }
  
  for (ptm_id in unique(combined_data$id.x)) {
    for (group in sample_groups) {
      acr_samples <- combined_data %>% filter(id.x == ptm_id) %>% filter(grepl(group, id.y))
      if (nrow(acr_samples) >=  min_values_present) {
        corr <- cor(acr_samples$value.prot, acr_samples$value, method = "pearson")
      } else {
        corr <- NA
      }
      
      store_id <- ifelse(group == "", ptm_id, paste0(ptm_id, " <", group, ">"))
      corr_store[store_id] <- corr
    }
  }
  
  return(corr_store)
}

ptm_corr_per_across_samples <- function(ptm_gct_path, ptm_normalized_gct_path, use_all_samples = FALSE, sample_groups = c("06h", "24h", "96h"), min_values_present = 4) {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_norm <- parse_gctx(ptm_normalized_gct_path)
  ptm_vals <- melt_gct(ptm)[, c("id.y", "id.x", "value")]
  ptm_norm_vals <- melt_gct(ptm_norm)[, c("id.y", "id.x", "value")]
  
  combined_ptm <- merge(ptm_vals, ptm_norm_vals, by = c("id.x", "id.y"))
  colnames(combined_ptm) <- c("id.x", "id.y", "value", "value.norm")
  corr_store <- list()
  
  if (use_all_samples) {
    sample_groups = c("")  # just so that it matches to everything
  }
  
  for (ptm_id in unique(combined_ptm$id.x)) {
    for (group in sample_groups) {
      acr_samples <- combined_ptm %>% filter(id.x == ptm_id) %>% filter(grepl(group, id.y))
      if (nrow(acr_samples) >=  min_values_present) {
        corr <- cor(acr_samples$value, acr_samples$value.norm, method = "pearson")
      } else {
        corr <- NA
      }
      
      store_id <- ifelse(group == "", ptm_id, paste0(ptm_id, " <", group, ">"))
      corr_store[store_id] <- corr
    }
  }
  
  return(corr_store)
}

