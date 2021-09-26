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


# TODO: how many sites at selected gene have |log| > 1 (all, for 6 and 24), average of each timepoint

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
  method = "pearson",
  center_data = TRUE)
{
  prot_gct <- parse_gctx(prot_gct_path)
  ptm_gct <- parse_gctx(ptm_gct_path)
  ptm_norm_gct <- parse_gctx(ptm_norm_gct_path)
  
  # TODO: in the future, store this centered proteome when we normalize
  if (center_data) {
    prot_gct@mat <- sweep(prot_gct@mat, 2, apply(prot_gct@mat, 2, function(col) { median(col, na.rm = TRUE) }), "-")
  }
  
  # obtain across samples correlations for prot and PTM
  prot_ptm_unnorm_corr <- prot_ptm_corr_across_samples(prot_gct, ptm_gct, groups_colname, subset_cond, min_n_values, method)
  colnames(prot_ptm_unnorm_corr) <- add_prefix_to_colnames("prot_ptm_unnorm", prot_ptm_unnorm_corr, except = "id.x")

  # obtain across samples correlations for prot and normalized PTM
  prot_ptm_norm_corr <- prot_ptm_corr_across_samples(prot_gct, ptm_norm_gct, groups_colname, subset_cond, min_n_values, method)
  colnames(prot_ptm_norm_corr) <- add_prefix_to_colnames("prot_ptm_norm", prot_ptm_norm_corr, except = "id.x")

  # obtain across samples correlations for PTM and normalized PTM
  ptm_unnorm_ptm_norm_corr <- ptm_corr_across_samples(ptm_gct, ptm_norm_gct, groups_colname, subset_cond, min_n_values, method)
  colnames(ptm_unnorm_ptm_norm_corr) <- add_prefix_to_colnames("ptm_unnorm_ptm_norm", ptm_unnorm_ptm_norm_corr, except = "id.x")

  comb_acr_corr <- merge(prot_ptm_unnorm_corr, prot_ptm_norm_corr, by = "id.x", all = TRUE)
  comb_acr_corr <- merge(comb_acr_corr, ptm_unnorm_ptm_norm_corr, by = "id.x", all = TRUE)

  save_folder <- dirname(ptm_norm_gct_path)
  acr_corr_path <- file.path(save_folder, "corr_across_samples.csv")
  write.csv(comb_acr_corr, acr_corr_path, row.names = FALSE)
  print(paste0("Per PTM site, across sample correlations saved at: ", acr_corr_path))
  
  prot_ptm_unnorm_corr_in_sample <- prot_ptm_corr_per_sample(prot_gct, ptm_gct)
  colnames(prot_ptm_unnorm_corr_in_sample) <- add_prefix_to_colnames("prot_ptm_unnorm", prot_ptm_unnorm_corr_in_sample, except = "id.y")
  
  prot_ptm_norm_corr_in_sample <- prot_ptm_corr_per_sample(prot_gct, ptm_norm_gct)
  colnames(prot_ptm_norm_corr_in_sample) <- add_prefix_to_colnames("prot_ptm_norm", prot_ptm_norm_corr_in_sample, except = "id.y")
  
  ptm_unnorm_ptm_norm_corr_in_sample <- ptm_corr_per_sample(ptm_gct, ptm_norm_gct)
  colnames(ptm_unnorm_ptm_norm_corr_in_sample) <- add_prefix_to_colnames("ptm_unnorm_ptm_norm", ptm_unnorm_ptm_norm_corr_in_sample, except = "id.y")

  comb_in_corr <- merge(prot_ptm_unnorm_corr_in_sample, prot_ptm_norm_corr_in_sample, by = "id.y", all = TRUE)
  comb_in_corr <- merge(comb_in_corr, ptm_unnorm_ptm_norm_corr_in_sample, by = "id.y", all = TRUE)
  
  in_corr_path <- file.path(save_folder, "corr_in_samples.csv")
  write.csv(comb_in_corr, in_corr_path, row.names = FALSE)
  print(paste0("Per sample, across all PTM sites correlations saved at: ", in_corr_path))
  
}

prot_ptm_corr_across_samples <- function(
  prot_gct,
  ptm_gct,
  groups_colname = NULL,  # name of column indicating grouping (e.g., "pert_time")
  subset_cond = NULL,
  min_n_values = 4,
  method = "pearson")
{
  prot <- prot_gct
  ptm <- match_ptm_to_proteome(ptm_gct, prot)
  comb_ptm_prot <- merge_ptm_prot_df(ptm, prot)
  if (!is.null(subset_cond)) {
    comb_ptm_prot <- comb_ptm_prot %>% filter(!!rlang::parse_expr(subset_cond))
  }
  corr_store <- list()
  
  # get correlations for all PTMs
  corr_store[["all"]] = list()
  for (ptm_id in unique(comb_ptm_prot$id.x)) {
    # print(ptm_id)
    acr_samples <- comb_ptm_prot[comb_ptm_prot$id.x == ptm_id, ]  # %>% filter(id.x == ptm_id)
    corr_store[["all"]][[ptm_id]] <- run_corr(acr_samples$value.prot, acr_samples$value, min_n_values, method = method)
  }
  
  # get correlations for PTMs by grouping
  if (!is.null(groups_colname)) {
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
    
    for (group in sample_groups) {
      group_name <- paste0(groups_colname, ":", group)
      corr_store[[group_name]] = list()
      
      for (ptm_id in unique(comb_ptm_prot$id.x)) {
        # print(groups_colname)
        # print(ptm_id)
        acr_samples <- comb_ptm_prot[(comb_ptm_prot$id.x == ptm_id & comb_ptm_prot[[groups_colname]] == group), ] # %>% filter(id.x == ptm_id) %>% filter(!!rlang::sym(groups_colname) == group)
        # TODO: save()
        corr_store[[group_name]][[ptm_id]] <- run_corr(acr_samples$value.prot, acr_samples$value, min_n_values, method = method)
      }
    }
  }
  
  corr_df <- corr_df_from_list(corr_store)
  return(corr_df)
}

ptm_corr_across_samples <- function(ptm_gct, ptm_normalized_gct, groups_colname = NULL, subset_cond = NULL, min_n_values = 4, method = "pearson") {
  ptm <- ptm_gct
  ptm_norm <- ptm_normalized_gct
  ptm_vals <- melt_gct(ptm)[, c("id.y", "id.x", "value")]
  ptm_norm_vals <- melt_gct(ptm_norm)[, c("id.y", "id.x", "value", ..groups_colname)]
  
  comb_ptm <- merge(ptm_vals, ptm_norm_vals, by = c("id.x", "id.y"))  # TOOD: tidy join fn
  colnames(comb_ptm) <- c("id.x", "id.y", "value", "value.norm", groups_colname)  # TODO: this hard-coded order is quite dangerous...
  
  if (!is.null(comb_ptm) & !is.null(subset_cond)) {
    comb_ptm <- comb_ptm %>% filter(!!rlang::parse_expr(subset_cond))
  }
  corr_store <- list()
  
  # get correlations for all PTMs
  corr_store[["all"]] = list()
  for (ptm_id in unique(comb_ptm$id.x)) {
    acr_samples <- comb_ptm[comb_ptm$id.x == ptm_id, ]  # %>% filter(id.x == ptm_id)
    corr_store[["all"]][[ptm_id]] <- run_corr(acr_samples$value, acr_samples$value.norm, min_n_values, method = method)
  }
  
  # get correlations for PTMs by grouping
  if (!is.null(groups_colname)) {
    sample_groups <- unique(comb_ptm[[groups_colname]])
    
    for (group in sample_groups) {
      group_name <- paste0(groups_colname, ":", group)
      corr_store[[group_name]] = list()
      # print(paste0("group: ", group_name))
      
      for (ptm_id in unique(comb_ptm$id.x)) {
        # print(paste0("ptm: ", ptm_id))
        acr_samples <- comb_ptm[(comb_ptm$id.x == ptm_id & comb_ptm[[groups_colname]] == group), ] # %>% filter(id.x == ptm_id) %>% filter(!!rlang::sym(groups_colname) == group)
        corr_store[[group_name]][[ptm_id]] <- run_corr(acr_samples$value, acr_samples$value.norm, min_n_values, method = method)
      }
    }
  }
  
  corr_df <- corr_df_from_list(corr_store)
  return(corr_df)
}

prot_ptm_corr_per_sample <- function(prot_gct, ptm_gct) {
  prot <- prot_gct
  ptm <- match_ptm_to_proteome(ptm_gct, prot)
  comb_ptm_prot <- merge_ptm_prot_df(ptm, prot)
  
  corr_store <- list()
  for (sample_name in unique(comb_ptm_prot$id.y)) {
    sample <- comb_ptm_prot %>% filter(id.y == sample_name)
    corr_store[sample_name] <- cor(sample$value.prot, sample$value, method = "pearson", use = "complete.obs")
  }
  
  corr_store <- stack(corr_store)
  colnames(corr_store) <- c("all", "id.y")
  return(corr_store)
}

ptm_corr_per_sample <- function(ptm_gct, ptm_normalized_gct) {
  ptm <- ptm_gct
  ptm_norm <- ptm_normalized_gct
  ptm_vals <- melt_gct(ptm)[, c("id.y", "id.x", "value")]
  ptm_norm_vals <- melt_gct(ptm_norm)[, c("id.y", "id.x", "value")]
  
  comb_ptm <- merge(ptm_vals, ptm_norm_vals, by = c("id.x", "id.y"))
  colnames(comb_ptm) <- c("id.x", "id.y", "value", "value.norm")
  
  corr_store <- list()
  for (sample_name in unique(comb_ptm$id.y)) {
    sample <- comb_ptm %>% filter(id.y == sample_name)
    corr_store[sample_name] <- cor(sample$value, sample$value.norm, method = "pearson", use = "complete.obs")
  }
  
  corr_store <- stack(corr_store)
  colnames(corr_store) <- c("all", "id.y")
  return(corr_store)
}

run_corr <- function(list_1, list_2, min_n_values = 4, method = "pearson") {
  if (sum(!(is.na(list_1) | is.na(list_2))) >= min_n_values) {
    corr <- cor(list_1, list_2, method = method, use = "complete.obs")
  } else {
    corr <- NA
  }
  
  return(corr)
}

ptm_log_fold_stats <- function(ptm_gct_path) {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_df <- as.data.frame(ptm@mat)
  
  mean_log_fold <- colMeans(ptm_df, na.rm = TRUE)
  log_fold_geq_1 <- colSums(ptm_df >= 1, na.rm = TRUE)
  log_fold_leq_neg1 <- colSums(ptm_df <= -1, na.rm = TRUE)
  not_na <- colSums(!is.na(ptm_df), na.rm = TRUE)
  log_fold_geq_1_frac <- log_fold_geq_1 / not_na
  log_fold_leq_neg1_frac <- log_fold_leq_neg1 / not_na
  
  res <- as.data.frame(rbind(mean_log_fold, log_fold_geq_1, log_fold_leq_neg1, not_na, log_fold_geq_1_frac, log_fold_leq_neg1_frac))
  res <- cbind(model_stat = c("mean_log_fold", "log_fold_geq_1", "log_fold_leq_neg1", "not_na", "log_fold_geq_1_frac", "log_fold_leq_neg1_frac"), res)
  rownames(res) <- NULL
  out_path <- file.path(dirname(ptm_gct_path), "ptm_log_fold_stats.csv")
  write.csv(res, out_path, row.names = FALSE)
}

ptm_log_fold_stats_per_group <- function(ptm_gct_path, groups_colname = NULL, groups_colname_2 = NULL) {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_df <- as.data.frame(ptm@mat)
  
  # for each sample, store a list PTM sites for which absolute value log fold change is greater or equal (geq) to 1 AND less or equal (leq) to -1
  log_fold_geq_1_store <- list()
  log_fold_leq_neg1_store <- list()
  not_na_store <- list()
  for (col in colnames(ptm_df)) {
    log_fold_geq_1_store[[col]] <- rownames(ptm_df)[which(ptm_df[col] >= 1)]
    log_fold_leq_neg1_store[[col]] <- rownames(ptm_df)[which(ptm_df[col] <= -1)]
    not_na_store[[col]] <- rownames(ptm_df)[which(!is.na(ptm_df[col]))]
  }
  
  mean_log_fold_store_per_group <- list()
  log_fold_geq_1_count_per_group <- list()
  log_fold_leq_neg1_count_per_group <- list()
  not_na_count_per_group <- list()
  # fractions of |log| >= 1
  log_fold_geq_1_frac_per_group <- list()
  log_fold_leq_neg1_frac_per_group <- list()
  
  mean_log_fold_store_per_group[["all"]] <- mean(colMeans(ptm_df, na.rm = TRUE))
  log_fold_geq_1_count_per_group[["all"]] <- mean(lengths(log_fold_geq_1_store))  # UNION: length(unique(flatten_chr(log_fold_geq_1_store)))
  log_fold_leq_neg1_count_per_group[["all"]] <- mean(lengths(log_fold_geq_1_store))  # UNION: length(unique(flatten_chr(log_fold_leq_neg1_store)))
  not_na_count_per_group[["all"]] <- length(unique(flatten_chr(not_na_store)))
  log_fold_geq_1_frac_per_group[["all"]] <- log_fold_geq_1_count_per_group[["all"]] / not_na_count_per_group[["all"]]
  log_fold_leq_neg1_frac_per_group[["all"]] <- log_fold_leq_neg1_count_per_group[["all"]] / not_na_count_per_group[["all"]]
  
  if (!is.null(groups_colname) & is.null(groups_colname_2)) {
    for (group in unique(ptm@cdesc[[groups_colname]])) {
      # which samples belong to the same group?
      samples <- c(rownames(ptm@cdesc %>% filter(!!rlang::sym(groups_colname) == group)))
      
      # consecutive union of lists of PTM sites within samples of the group
      group_name <- paste0(groups_colname, ":", group)
      mean_log_fold_store_per_group[[group_name]] <- mean(colMeans(ptm_df[samples], na.rm = TRUE))
      log_fold_geq_1_count_per_group[[group_name]] <- mean(lengths(log_fold_geq_1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_geq_1_store[samples])))
      log_fold_leq_neg1_count_per_group[[group_name]] <- mean(lengths(log_fold_leq_neg1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_leq_neg1_store[samples])))
      not_na_count_per_group[[group_name]] <- length(unique(flatten_chr(not_na_store[samples])))
      log_fold_geq_1_frac_per_group[[group_name]] <- log_fold_geq_1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
      log_fold_leq_neg1_frac_per_group[[group_name]] <- log_fold_leq_neg1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
    }
  } else if (!is.null(groups_colname) & !is.null(groups_colname_2)) {
    for (group in unique(ptm@cdesc[[groups_colname]])) {
      for (group_2 in unique(ptm@cdesc[[groups_colname_2]])) {
        samples <- c(rownames(ptm@cdesc %>% filter(!!rlang::sym(groups_colname) == group) %>% filter(!!rlang::sym(groups_colname_2) == group_2)))
        
        group_name <- paste0(group, ".", group_2)
        mean_log_fold_store_per_group[[group_name]] <- mean(colMeans(ptm_df[samples], na.rm = TRUE))
        log_fold_geq_1_count_per_group[[group_name]] <- mean(lengths(log_fold_geq_1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_geq_1_store[samples])))
        log_fold_leq_neg1_count_per_group[[group_name]] <- mean(lengths(log_fold_leq_neg1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_leq_neg1_store[samples])))
        not_na_count_per_group[[group_name]] <- length(unique(flatten_chr(not_na_store[samples])))
        log_fold_geq_1_frac_per_group[[group_name]] <- log_fold_geq_1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
        log_fold_leq_neg1_frac_per_group[[group_name]] <- log_fold_leq_neg1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
      }
    }
  }
  
  comb_metrics <- rbind(
    as.data.frame(mean_log_fold_store_per_group),
    as.data.frame(log_fold_geq_1_count_per_group),
    as.data.frame(log_fold_leq_neg1_count_per_group),
    as.data.frame(not_na_count_per_group),
    as.data.frame(log_fold_geq_1_frac_per_group),
    as.data.frame(log_fold_leq_neg1_frac_per_group)
  )
  comb_metrics <- cbind("model_stat" = c("mean_log_fold", "log_fold_geq_1", "log_fold_leq_neg1", "not_na", "log_fold_geq_1_frac", "log_fold_leq_neg1_frac"), comb_metrics)
  if (!is.null(groups_colname_2)) {
    out_path <- file.path(dirname(ptm_gct_path), paste0("ptm_log_fold_stats", "-", groups_colname, "-", groups_colname_2, ".csv"))
  } else {
    out_path <- file.path(dirname(ptm_gct_path), "ptm_log_fold_stats_per_group.csv")
  }
  write.csv(comb_metrics, out_path, row.names = FALSE)
}

ptm_log_fold_stats_per_group_targeted <- function(ptm_gct_path, groups_colname = NULL, groups_colname_2 = NULL, cell_ids = NULL, target_accessions = NULL) {
  ptm <- parse_gctx(ptm_gct_path)
  ptm_df <- as.data.frame(ptm@mat)  # instead make = NA where column thing is not equal to target accession

  for (col in colnames(ptm_df)) {
    accessions_filter <- target_accessions[which(sapply(cell_ids, function(x) { grepl(x, col)}))]
    not_target_sites <- !grepl(accessions_filter, rownames(ptm_df[col]))
    not_na <- !is.na(ptm_df[col])
    ptm_df[not_target_sites & not_na , col] <- NA  # to count only respective target sites
  }
  
  # for each sample, store a list PTM sites for which absolute value log fold change is greater or equal (geq) to 1 AND less or equal (leq) to -1
  log_fold_geq_1_store <- list()
  log_fold_leq_neg1_store <- list()
  not_na_store <- list()
  for (col in colnames(ptm_df)) {
    accession <- target_accessions[which(sapply(cell_ids, function(x) { grepl(x, col)}))]  # fetch the corresponding accession to search for
    targets <- ptm_df[grepl(accession, rownames(ptm_df)), col]
    
    log_fold_geq_1_store[[col]] <- rownames(ptm_df)[which(targets >= 1)]
    log_fold_leq_neg1_store[[col]] <- rownames(ptm_df)[which(targets<= -1)]
    not_na_store[[col]] <- rownames(ptm_df)[which(!is.na(targets))]
  }
  
  mean_log_fold_store_per_group <- list()
  log_fold_geq_1_count_per_group <- list()
  log_fold_leq_neg1_count_per_group <- list()
  not_na_count_per_group <- list()
  # fractions of |log| >= 1
  log_fold_geq_1_frac_per_group <- list()
  log_fold_leq_neg1_frac_per_group <- list()
  
  # mean_log_fold_store_per_group[["all"]] <- mean(colMeans(ptm_df, na.rm = TRUE))
  # log_fold_geq_1_count_per_group[["all"]] <- mean(lengths(log_fold_geq_1_store))  # UNION: length(unique(flatten_chr(log_fold_geq_1_store)))
  # log_fold_leq_neg1_count_per_group[["all"]] <- mean(lengths(log_fold_geq_1_store))  # UNION: length(unique(flatten_chr(log_fold_leq_neg1_store)))
  # not_na_count_per_group[["all"]] <- length(unique(flatten_chr(not_na_store)))
  # log_fold_geq_1_frac_per_group[["all"]] <- log_fold_geq_1_count_per_group[["all"]] / not_na_count_per_group[["all"]]
  # log_fold_leq_neg1_frac_per_group[["all"]] <- log_fold_leq_neg1_count_per_group[["all"]] / not_na_count_per_group[["all"]]
  
  if (!is.null(groups_colname) & is.null(groups_colname_2)) {
    for (group in unique(ptm@cdesc[[groups_colname]])) {
      # which samples belong to the same group?
      samples <- c(rownames(ptm@cdesc %>% filter(!!rlang::sym(groups_colname) == group)))
      
      # consecutive union of lists of PTM sites within samples of the group
      group_name <- paste0(groups_colname, ":", group)
      mean_log_fold_store_per_group[[group_name]] <- mean(colMeans(ptm_df[samples], na.rm = TRUE))
      log_fold_geq_1_count_per_group[[group_name]] <- mean(lengths(log_fold_geq_1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_geq_1_store[samples])))
      log_fold_leq_neg1_count_per_group[[group_name]] <- mean(lengths(log_fold_leq_neg1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_leq_neg1_store[samples])))
      not_na_count_per_group[[group_name]] <- length(unique(flatten_chr(not_na_store[samples])))
      log_fold_geq_1_frac_per_group[[group_name]] <- log_fold_geq_1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
      log_fold_leq_neg1_frac_per_group[[group_name]] <- log_fold_leq_neg1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
    }
  } else if (!is.null(groups_colname) & !is.null(groups_colname_2)) {
    for (group in unique(ptm@cdesc[[groups_colname]])) {
      for (group_2 in unique(ptm@cdesc[[groups_colname_2]])) {
        samples <- c(rownames(ptm@cdesc %>% filter(!!rlang::sym(groups_colname) == group) %>% filter(!!rlang::sym(groups_colname_2) == group_2)))
        
        group_name <- paste0(group, ".", group_2)
        mean_log_fold_store_per_group[[group_name]] <- mean(colMeans(ptm_df[samples], na.rm = TRUE))
        log_fold_geq_1_count_per_group[[group_name]] <- mean(lengths(log_fold_geq_1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_geq_1_store[samples])))
        log_fold_leq_neg1_count_per_group[[group_name]] <- mean(lengths(log_fold_leq_neg1_store[samples]))  # UNION: length(unique(flatten_chr(log_fold_leq_neg1_store[samples])))
        not_na_count_per_group[[group_name]] <- length(unique(flatten_chr(not_na_store[samples])))
        log_fold_geq_1_frac_per_group[[group_name]] <- log_fold_geq_1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
        log_fold_leq_neg1_frac_per_group[[group_name]] <- log_fold_leq_neg1_count_per_group[[group_name]] / not_na_count_per_group[[group_name]]
      }
    }
  }
  
  comb_metrics <- rbind(
    as.data.frame(mean_log_fold_store_per_group),
    as.data.frame(log_fold_geq_1_count_per_group),
    as.data.frame(log_fold_leq_neg1_count_per_group),
    as.data.frame(not_na_count_per_group),
    as.data.frame(log_fold_geq_1_frac_per_group),
    as.data.frame(log_fold_leq_neg1_frac_per_group)
  )
  comb_metrics <- cbind("model_stat" = c("mean_log_fold", "log_fold_geq_1", "log_fold_leq_neg1", "not_na", "log_fold_geq_1_frac", "log_fold_leq_neg1_frac"), comb_metrics)
  if (!is.null(groups_colname_2)) {
    out_path <- file.path(dirname(ptm_gct_path), paste0("ptm_log_fold_stats", "-", groups_colname, "-", groups_colname_2, "-targeted.csv"))
  } else {
    out_path <- file.path(dirname(ptm_gct_path), "ptm_log_fold_stats_per_group-targeted.csv")
  }
  write.csv(comb_metrics, out_path, row.names = FALSE)
}
