#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Normalizes ptm to protein level by fitting a global linear model and returning residuals
# Harry Kane, Pierre Beltran, D. R. Mani, Karsten Krug, Zachery Gillette, Surya Mani, Khoi Pham Munchic
#
# Fits ptm = beta_0 + beta_1*protein to all matched points in a dataset, and returns residuals as
# protein-corrected ptm values.
#
#########################################################################################################


library("pacman")
p_load(cmapR)
p_load(reshape)
p_load(yaml)
p_load(tidyr)
p_load(limma)
p_load(lme4)
p_load(stats)

source("/prot/proteomics/Projects/PGDAC/src/helper.R")
source("/prot/proteomics/Projects/PGDAC/src/metrics.R")
source("/prot/proteomics/Projects/PGDAC/src/analyze.R")


#' @title Normalize PTM levels with respect to protein levels
#'
#' @param proteome.gct Path to proteome GCT inside of "`data_dir`/input"
#' @param ptm.gct Path to PTM GCT inside of "`data_dir`/input"
#' @param data_dir Path to directory with "input/" subdirectory containing input GCTs
#' @param output.prefix Output file name (default: NULL)
#'
#' @return Long format PTM site level and cognate protein level
normalize_ptm <- function(proteome.gct, ptm.gct, data_dir, output.prefix = NULL,
                          try.all.accession.numbers = TRUE,              # try hard to find a match (using accession_numbers)
                          accession_number = "accession_number",         # column with protein/ptm accession number
                          accession_numbers = "accession_numbers",       # accession_numbers for protein/ptm group
                          accession_sep = "|",                           # separator for each accession number in accession_numbers
                          score = "scoreUnique",                         # column with protein scores
                          ndigits = 5,                                   # precision level

                          center_data = FALSE,                           # whether to center the data with median
                          norm_method = "global",                        # normalization method (options: "global", "per_ptm_site", "subtract", "mixed_eff_protein", "mixed_eff_ptm_site", "mixed_eff_nested")
                          lm_method = "ls",                              # method argument of `lmFit`, "ls" for least squares, "robust" for M-estimation
                          lm_formula = value ~ value.prot,               # formula for linear regression (MUST include value as predicted variable and value.prot as explanatory)
                          loess = FALSE,                                 # whether to use loess regression (applicable to global, per protein, and per PTM models)
                          model_name = "proteome_relative_norm")         # subfolder name in which to store results)

                          groups_colname = NULL,                         # NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
                          subset_cond = NULL,                            # string: logical condition by which to subset the data (e.g., "as.character(pert_time) %in% c('6', '24')")
                          min_n_values = 4,                              # per_ptm_site: what least number of samples must contain both non-NA PTM and protein values

{
  # import GCT files
  proteome <- parse_gctx(file.path(data_dir, "input", proteome.gct))
  temp_ptm <- parse_gctx(file.path(data_dir, "input", ptm.gct))
  
  if (center_data) {
    proteome@mat <- sweep(proteome@mat, 2, apply(proteome@mat, 2, function(col) { median(col, na.rm = TRUE) }), "-")
    temp_ptm@mat <- sweep(temp_ptm@mat, 2, apply(temp_ptm@mat, 2, function(col) { median(col, na.rm = TRUE) }), "-")
  }
  ptm <- match_ptm_to_proteome(temp_ptm, proteome, accession_number, accession_numbers, try.all.accession.numbers)
  
  # warn if not all samples are matched
  matched.samples <- intersect(ptm@cid, proteome@cid)
  if (length(matched.samples) != length(ptm@cid)) {
    warning ('WARNING: not all samples in ptm file have matches in proteome ... unmatched samples removed.')
  }
  
  comb_ptm_prot <- merge_ptm_prot_df(ptm, proteome, accession_number)
  if (!is.null(groups_colname)) {
    comb_ptm_prot[[groups_colname]] <- as.factor(comb_ptm_prot[[groups_colname]])
  }
  if (!is.null(subset_cond)) {
    comb_ptm_prot <- comb_ptm_prot %>% filter(!!rlang::parse_expr(subset_cond))
    print(paste0("Based on subset condition, ", nrow(comb_ptm_prot), " unique sites are included."))
  }
  
  # create directories to store normalized data
  dir.create(file.path(data_dir, "out"), showWarnings = FALSE)
  model_out_dir <- file.path(data_dir, "out", model_name)
  dir.create(model_out_dir, showWarnings = FALSE)
  
  file_name <- ifelse(!is.null (output.prefix), output.prefix,
                      unlist(strsplit(basename(ptm.gct), split = ".gct", fixed = TRUE))[1])
  out_path <- file.path(model_out_dir, file_name)
  
  # fits linear model and returns updated GCT
  if (norm_method == "global") {
    norm_vals <- normalize_global(comb_ptm_prot, lm_method, lm_formula, loess, groups_colname, min_n_values, out_path)
  } else if (norm_method == "per_ptm_site" | norm_method == "per_ptm_site") {
    norm_vals <- normalize_per_ptm_site(comb_ptm_prot, lm_method, lm_formula, loess, groups_colname, min_n_values, out_path)
  } else if (norm_method == "per_protein") {
    norm_vals <- normalize_per_protein(comb_ptm_prot, lm_method, lm_formula, loess, groups_colname, min_n_values, out_path)
  } else if (norm_method == "mixed_eff_protein") {
    norm_vals <- normalize_mixed_eff_prot(comb_ptm_prot, groups_colname, min_n_values, out_path)
  } else if (norm_method == "mixed_eff_ptm_site") {
    norm_vals <- normalize_mixed_eff_ptm(comb_ptm_prot, groups_colname, min_n_values, out_path)
  } else if (norm_method == "mixed_eff_nested") {
    norm_vals <- normalize_mixed_eff_nested(comb_ptm_prot, groups_colname, min_n_values, out_path)
  } else if (norm_method == "subtract") {
    norm_vals <- normalize_subtract(comb_ptm_prot)
  }
  ptm.norm <- update_gct(ptm, norm_vals)
  print ("Success.")
  
  # writes and returns updated GCT
  write_gct(ptm.norm, paste0(out_path, ".gct"), 
            appenddim = FALSE, precision = ndigits)
  invisible (ptm.norm)
  
  return(comb_ptm_prot)
}


#' @title Global normalization
#' 
#' @description Linear regression to correct PTM levels for underlying protein levels
#' Uses all pairs of PTM level and cognate protein level across samples to build regression
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#' @param lm_method Method argument of `lmFit`, "ls" for least squares, "robust" for M-estimation
#' @param lm_formula Formula for linear regression (MUST include value as predicted variable and value.prot as explanatory)
#' @param loess Use local polynomial regression for non-linear fit (piece-wise)
#' @param groups_colname NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
#' @param min_n_values Smallest number of samples must contain both non-NA PTM and protein values to build regression
#' @param out_path Path to save the normalized PTM levels GCT 
#'
#' @return Long format site quant and cognate protein quant
normalize_global <- function(comb_ptm_prot, lm_method, lm_formula, loess, groups_colname, min_n_values, out_path) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  # fit global model
  print ("Fitting model...") 
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  regr_store <- list()
  count_fail_ptm <- 0  # count PTM sites and samples that can't be normalized
  for (group in sample_groups) {
    if (is.na(group)) {
      samples <- comb_ptm_prot
      group <- "all"
    } else {
      samples <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    }
    
    if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
      result <- samples[ , c("id.x", "id.y")]
      if (!loess) {
        lm_design <- model.matrix(lm_formula, samples)
        model <- lmFit(samples$value, lm_design, method = lm_method)  # lm(samples$value ~ samples$value.prot, na.action = na.exclude)
        fitted_ptm_levels <- unlist(as.list(fitted(model)))
      } else {
        model <- lowess(samples$value.prot, samples$value)
        fitted_ptm_levels <- spline(model$x, model$y, xout = samples$value.prot)$y
      }
      result$residuals <- samples$value - fitted_ptm_levels  # result$residuals <- residuals(model)
      
      all_results <- rbind(all_results, result)
      count_fail_ptm <- count_fail_ptm + sum(is.na(result$residuals))
      regr_store[[group]] <- model
    }
  }
  
  print(paste0("Number of sites normalized: ", nrow(comb_ptm_prot) - count_fail_ptm, "/", nrow(comb_ptm_prot)))
  saveRDS(regr_store, paste0(out_path, ".RDS"))
  return(all_results)
}


#' @title Per PTM site normalization
#' 
#' @description Linear regression to correct PTM levels for underlying protein levels for each PTM site
#' For each unique PTM site, get PTM levels and cognate protein levels across samples to build regression
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#' @param lm_method Method argument of `lmFit`, "ls" for least squares, "robust" for M-estimation
#' @param lm_formula Formula for linear regression (MUST include value as predicted variable and value.prot as explanatory)
#' @param loess Use local polynomial regression for non-linear fit (piece-wise)
#' @param groups_colname NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
#' @param min_n_values Smallest number of samples must contain both non-NA PTM and protein values to build regression
#' @param out_path Path to save the normalized PTM levels GCT 
#'
#' @return Long format site quant and cognate protein quant
normalize_per_ptm_site <- function(comb_ptm_prot, lm_method, lm_formula, loess, groups_colname, min_n_values, out_path) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  all_ptms <- unique(comb_ptm_prot$id.x)
  tot_num_regr <- length(all_ptms) * length(sample_groups)
  print(paste("Found", length(all_ptms), "PTM-protein pairs and", length(sample_groups),
              "sample group(s);", tot_num_regr, "regressions to be built", sep = " "))
  
  # create a dataframe placeholder to store results in
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  count <- 1
  count_norm_ptm <- 0  # count PTM sites and samples that are normalized
  count_fail_regr <- 0  # count number of regressions that can't be built
  # loop through all PTMs, build regression for each
  regr_store <- list()
  for (group in sample_groups) {
    if (!is.na(group)) {
      samples_grouped <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    } else {
      group <- "all"
      samples_grouped <- comb_ptm_prot
    }
    regr_store[[group]] = list()
    
    for (ptm_site in all_ptms) {
      if (count == 1 | count %% 500 == 0) {
        print(paste0("Building regressions... (", count, "/", tot_num_regr, ")"))
      }
      samples <- samples_grouped %>% filter(id.x == ptm_site)
      result <- samples[ , c("id.x", "id.y")]
      
      if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
        
        if (!loess) {
          lm_design <- model.matrix(lm_formula, samples)
          model <- lmFit(samples$value, lm_design, method = lm_method)  # lm(samples$value ~ samples$value.prot, na.action = na.exclude)
          fitted_ptm_levels <- unlist(as.list(fitted(model)))
          result$residuals <- samples$value - fitted_ptm_levels  # result$residuals <- residuals(model)
        } else {
          model <- loess(lm_formula, data = samples)
          result$residuals <- residuals(model)
        }
        regr_store[[group]][[ptm_site]] <- model
        count_norm_ptm <- count_norm_ptm + sum(!is.na(result$residuals))
        
      } else {
        result$residuals <- NA
        count_fail_regr <- count_fail_regr + 1
      }
      
      all_results <- rbind(all_results, result)
      count <- count + 1
    }
  }
  
  print(paste0("Number of sites normalized: ", count_norm_ptm, "/", nrow(comb_ptm_prot)))
  print(paste0("Regressions failed (values present < ", min_n_values, "): ", count_fail_regr))
  saveRDS(regr_store, paste0(out_path, ".RDS"))
  return(all_results)
}


#' @title Per protein normalization
#' 
#' @description Linear regression to correct PTM levels for underlying protein levels for each PTM site
#' For each unique protein, get PTM levels on that protein and protein levels across samples to build regression
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#' @param groups_colname NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
#' @param min_n_values Smallest number of samples must contain both non-NA PTM and protein values to build regression
#' @param out_path Path to save the normalized PTM levels GCT 
#'
#' @return Long format dataframe | id.x (PTM) | id.y (protein) | residual |
normalize_per_protein <- function(comb_ptm_prot, lm_method, lm_formula, loess, groups_colname, min_n_values, out_path) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  all_prot <- unique(comb_ptm_prot[, accession_number])  # avoid hard-coding "accession_number" column
  tot_num_regr <- length(all_prot) * length(sample_groups)
  print(paste("Found", length(all_prot), "unique proteins and", length(sample_groups),
              "sample group(s);", tot_num_regr, "regressions to be built", sep = " "))
  
  # create a dataframe placeholder to store results in
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  count <- 1
  count_norm_ptm <- 0  # count proteins and samples that are normalized
  count_fail_regr <- 0  # count number of regressions that can't be built
  # loop through all proteins, build regression for each
  regr_store <- list()
  for (group in sample_groups) {
    if (!is.na(group)) {
      samples_grouped <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    } else {
      group <- "all"
      samples_grouped <- comb_ptm_prot
    }
    regr_store[[group]] <- list()
    
    for (prot in all_prot) {
      if (count == 1 | count %% 500 == 0) {
        print(paste0("Building regressions... (", count, "/", tot_num_regr, ")"))
      }
      samples <- samples_grouped %>% filter(accession_number == prot)
      result <- samples[ , c("id.x", "id.y")]
      
      if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
        
        if (!loess) {
          lm_design <- model.matrix(lm_formula, samples)
          model <- lmFit(samples$value, lm_design, method = lm_method)  # lm(samples$value ~ samples$value.prot, na.action = na.exclude)
          fitted_ptm_levels <- unlist(as.list(fitted(model)))
          result$residuals <- samples$value - fitted_ptm_levels  # result$residuals <- residuals(model)
        } else {
          model <- loess(lm_formula, data = samples)
          result$residuals <- residuals(model)
        }
        regr_store[[group]][[prot]] <- model
        count_norm_ptm <- count_norm_ptm + sum(!is.na(result$residuals))
        
      } else {
        result$residuals <- NA
        count_fail_regr <- count_fail_regr + 1
      }
      
      all_results <- rbind(all_results, result)
      count <- count + 1
    }
  }
  
  print(paste0("Number of PTMs normalized: ", count_norm_ptm, "/", nrow(comb_ptm_prot)))
  print(paste0("Regressions failed (values present < ", min_n_values, "): ", count_fail_regr))
  saveRDS(regr_store, paste0(out_path, ".RDS"))
  return(all_results)
}

#' @title Subtract normalization
#' 
#' @description Subtraction to correct PTM levels for underlying protein levels for each PTM site
#' For each PTM site, subtract cognate protein level from site level
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#'
#' @return Long format dataframe | id.x (PTM) | id.y (protein) | residual |
normalize_subtract <- function(comb_ptm_prot) {
  all_results <- comb_ptm_prot[ , c("id.x", "id.y")]
  all_results$residuals <- comb_ptm_prot$value - comb_ptm_prot$value.prot
  
  all_not_na_count <- sum(!is.na(comb_ptm_prot$value) | !is.na(comb_ptm_prot$value.prot))
  norm_not_na_count <- sum(!is.na(all_results$residuals))
  print(paste0("Successfully normalized ", norm_not_na_count, "/", all_not_na_count, " sites."))
  
  return(all_results)
}

#' @title Mixed-effects per protein normalization
#' 
#' @description Mixed effects linear regression to correct PTM levels for underlying protein levels for each PTM site
#' Fixed effect: For each protein, get PTM levels on protein and protein levels across samples to build regression
#' Random effect: Each PTM site has its own random effect
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#' @param groups_colname NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
#' @param min_n_values Smallest number of samples must contain both non-NA PTM and protein values to build regression
#' @param out_path Path to save the normalized PTM levels GCT 
#'
#' @return #' @return Long format dataframe | id.x (PTM) | id.y (protein) | residual |
normalize_mixed_eff_prot <- function(comb_ptm_prot, groups_colname, min_n_values, out_path) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  print ("Fitting model...") 
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  regr_store <- list()
  count_fail_ptm <- 0  # count PTM sites and samples that can't be normalized
  for (group in sample_groups) {
    if (is.na(group)) {
      samples <- comb_ptm_prot
      group <- "all"
    } else {
      samples <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    }
    
    if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
      model <- lmer(formula = value ~ value.prot + (value.prot | accession_number),
                    data = samples,
                    control = lmerControl(optimizer = 'Nelder_Mead'),
                    REML = T)
  
      result <- samples[ , c("id.x", "id.y")]
      result$residuals <- samples$value - fitted(model)
      
      all_results <- rbind(all_results, result)
      count_fail_ptm <- count_fail_ptm + sum(is.na(result$residuals))
      regr_store[[group]] <- model
    }
  }
  
  print(paste0("Number of sites normalized: ", nrow(comb_ptm_prot) - count_fail_ptm, "/", nrow(comb_ptm_prot)))
  saveRDS(regr_store, paste0(out_path, ".RDS"))
  return(all_results)
}


#' @title Mixed-effects per protein normalization
#' 
#' @description Mixed effects linear regression to correct PTM levels for underlying protein levels for each PTM site
#' Fixed effect: For each PTM site, get PTM levels and cognate protein levels across samples to build regression
#' Random effect: Protein ID serves random effect
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#' @param groups_colname NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
#' @param min_n_values Smallest number of samples must contain both non-NA PTM and protein values to build regression
#' @param out_path Path to save the normalized PTM levels GCT 
#'
#' @return Long format dataframe | id.x (PTM) | id.y (protein) | residual |
normalize_mixed_eff_ptm <- function(comb_ptm_prot, groups_colname, min_n_values, out_path) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  print ("Fitting model...") 
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  regr_store <- list()
  count_fail_ptm <- 0  # count PTM sites and samples that can't be normalized
  for (group in sample_groups) {
    if (is.na(group)) {
      samples <- comb_ptm_prot
      group <- "all"
    } else {
      samples <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    }
    
    if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
      model <- lmer(formula = value ~ value.prot + (value.prot | id.x),  # check if it was per ptm level
                    data = samples,
                    control = lmerControl(optimizer = 'bobyqa'),  # try different optimizer 
                    REML = T)
      
      result <- samples[ , c("id.x", "id.y")]
      result$residuals <- samples$value - fitted(model)
      
      all_results <- rbind(all_results, result)
      count_fail_ptm <- count_fail_ptm + sum(is.na(result$residuals))
      regr_store[[group]] <- model
    }
  }
  
  print(paste0("Number of sites normalized: ", nrow(comb_ptm_prot) - count_fail_ptm, "/", nrow(comb_ptm_prot)))
  saveRDS(regr_store, paste0(out_path, ".RDS"))
  return(all_results)
}


#' @title Nested mixed-effects per PTM site
#' 
#' @description Mixed effects linear regression to correct PTM levels for underlying protein levels for each PTM site
#' Fixed effect: For each PTM site, get PTM levels and cognate protein levels across samples to build regression
#' Random effect: Protein ID serves random effect, and PTM site ID serves as random effect nested within protein random effect
#'
#' @param comb_ptm_prot Long format PTM site level and cognate protein level
#' @param groups_colname NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
#' @param min_n_values Smallest number of samples must contain both non-NA PTM and protein values to build regression
#' @param out_path Path to save the normalized PTM levels GCT 
#'
#' @return Long format dataframe | id.x (PTM) | id.y (protein) | residual |
normalize_mixed_eff_nested <- function(comb_ptm_prot, groups_colname, min_n_values, out_path) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- unique(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  print ("Fitting model...") 
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  regr_store <- list()
  count_fail_ptm <- 0  # count PTM sites and samples that can't be normalized
  for (group in sample_groups) {
    if (is.na(group)) {
      samples <- comb_ptm_prot
      group <- "all"
    } else {
      samples <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    }
    
    if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
      # id.x is unique to each accession_number so there is already implicit nesting
      model <- lmer(formula = value ~ value.prot + (value.prot | accession_number) + (value.prot | id.x),
                    data = samples,
                    control = lmerControl(optimizer = 'Nelder_Mead'),
                    REML = T)
      
      result <- samples[ , c("id.x", "id.y")]
      result$residuals <- samples$value - fitted(model)
      
      all_results <- rbind(all_results, result)
      count_fail_ptm <- count_fail_ptm + sum(is.na(result$residuals))
      regr_store[[group]] <- model
    }
  }
  
  print(paste0("Number of sites normalized: ", nrow(comb_ptm_prot) - count_fail_ptm, "/", nrow(comb_ptm_prot)))
  saveRDS(regr_store, paste0(out_path, ".RDS"))
  return(all_results)
}

if (!interactive()) {
  ## call via command line
  ## usage: Rscript normalize-ptm.R <proteome.gct> <ptm.gct> <output.prefix> <yaml.file>
  
  # process command line args
  args <- commandArgs (TRUE)
  proteome_gct <- as.character(args[1])
  ptm_gct <- as.character(args[2])
  output_prefix <- as.character (args[3])
  if (output_prefix == "NULL") output_prefix <- NULL
  yaml_file <- args[4]
  
  # read yaml parameters
  yaml_params <- read_yaml(yaml_file)
  accession_number_col <- yaml_params$panoply_ptm_normalization$accession_number_colname
  accession_numbers_col <- yaml_params$panoply_ptm_normalization$accession_numbers_colname
  accession_numbers_sep <- yaml_params$panoply_ptm_normalization$accession_numbers_separator
  score_col <- yaml_params$panoply_ptm_normalization$score_colname
  ndigits <- yaml_params$global_parameters$output_precision$ndigits
  
  
  # call ptm normalization
  normalize_ptm (proteome.gct=proteome_gct, ptm.gct=ptm_gct, output.prefix=output_prefix,
                 try.all.accession.numbers=ifelse (accession_numbers_col=="NULL", FALSE, TRUE), 
                 accession_number=accession_number_col,
                 accession_numbers=ifelse (accession_numbers_col=="NULL", NULL, accession_numbers_col),
                 accession_sep=accession_numbers_sep,
                 score=ifelse (score_col=="NULL", NULL, score_col),
                 ndigits=ndigits)
}
