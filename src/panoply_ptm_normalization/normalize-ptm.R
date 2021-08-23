#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Normalizes ptm to protein level by fitting a global linear model and returning residuals
# Harry Kane, Pierre Beltran, D. R. Mani, Karsten Krug, Zachery Gillette, Surya Mani, Khoi Pham (Munchic)
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

source("panoply_ptm_normalization/helper.R")
source("panoply_ptm_normalization/metrics.R")
source("panoply_ptm_normalization/analyze.R")


normalize_ptm <- function(proteome.gct, ptm.gct, output.prefix = NULL, 
                          try.all.accession.numbers = TRUE,              # try hard to find a match (using accession_numbers)
                          accession_number = "accession_number",         # column with protein/ptm accession number
                          accession_numbers = "accession_numbers",       # accession_numbers for protein/ptm group
                          accession_sep = "|",                           # separator for each accession number in accession_numbers
                          score = "scoreUnique",                         # column with protein scores
                          ndigits = 5,                                   # precision level
                          norm_method = "global",                        # normalization method (options: global, pairwise)
                          lm_method = "ls",                              # method argument of `lmFit`, "ls" for least squares, "robust" for M-estimation
                          lm_formula = value ~ value.prot,               # formula for linear regression (MUST include value as predicted variable and value.prot as explanatory)
                          groups_colname = NULL,                         # NULL if want to use all samples, otherwise name of column indicating grouping (e.g., "pert_time")
                          subset_cond = NULL,                            # string: logical condition by which to subset the data (e.g., "as.character(pert_time) %in% c('6', '24')")
                          min_n_values = 4,                              # pairwise: what least number of samples must contain both non-NA PTM and protein values
                          data_dir = "data",
                          model_name = "proteome_relative_norm")   # postfix to append to the file name of the GCT to be normalized
{
  # import GCT files
  proteome <- parse_gctx(file.path(data_dir, "input", proteome.gct))
  temp_ptm <- parse_gctx(file.path(data_dir, "input", ptm.gct))
  ptm <- match_ptm_to_proteome(temp_ptm, proteome, accession_number, accession_numbers, try.all.accession.numbers)
  
  # warn if not all samples are matched
  matched.samples <- intersect (ptm@cid, proteome@cid)
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
  
  # fits linear model and returns updated GCT
  if (norm_method == "global") {
    norm_vals <- normalize_global(comb_ptm_prot, lm_method, lm_formula, groups_colname, min_n_values)
  } else if (norm_method == "pairwise") {
    norm_vals <- normalize_pairwise(comb_ptm_prot, lm_method, lm_formula, groups_colname, min_n_values)
  }
  ptm.norm <- update_gct(ptm, norm_vals)
  print ("Success.")
  
  # writes and returns updated GCT
  file_name <- ifelse(!is.null (output.prefix), output.prefix,
                        unlist(strsplit(basename(ptm.gct), split = ".gct", fixed = TRUE))[1])
  # TODO: wow, such a hardcoding of directory :(
  dir.create(file.path(data_dir, "out"), showWarnings = FALSE)
  dir.create(file.path(data_dir, "out", model_name), showWarnings = FALSE)
  out_path <- file.path(data_dir, "out", model_name, file_name)
  
  write_gct(ptm.norm, paste0(out_path, ".gct"), 
            appenddim = FALSE, precision = ndigits)
  invisible (ptm.norm)
}


#' Applies linear regression to correct PTM levels for underlying protein levels
normalize_global <- function (comb_ptm_prot, lm_method, lm_formula, groups_colname, min_n_values) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- levels(comb_ptm_prot[[groups_colname]])
  } else {
    sample_groups <- c(NA)  # to use all samples
  }
  
  # fit global model
  print ("Fitting model...") 
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  count_fail_ptm <- 0  # count PTM sites and samples that can't be normalized
  for (group in sample_groups) {
    if (is.na(group)) {
      samples <- comb_ptm_prot
    } else {
      samples <- comb_ptm_prot %>% filter(!!rlang::sym(groups_colname) == group)
    }
    
    if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
      lm_design <- model.matrix(lm_formula, samples)
      model <- lmFit(samples$value, lm_design, method = lm_method)
        
      result <- samples[ , c("id.x", "id.y")]
      fitted_ptm_levels <- unlist(as.list(fitted(model)))
      result$residuals <- samples$value - fitted_ptm_levels
      all_results <- rbind(all_results, result)
      count_fail_ptm <- count_fail_ptm + sum(is.na(result$residuals))
    }
  }

  print(paste0("Number of sites normalized: ", nrow(comb_ptm_prot) - count_fail_ptm, "/", nrow(comb_ptm_prot)))
  return(all_results)
}


#' Applies linear regression to correct PTM levels for underlying protein levels
#' for each PTM-protein pair across selected or all samples
normalize_pairwise <- function(comb_ptm_prot, lm_method, lm_formula, groups_colname, min_n_values) {
  if (!is.null(groups_colname)) {  # if use all samples
    sample_groups <- levels(comb_ptm_prot[[groups_colname]])
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
  for (ptm_site in all_ptms) {
    for (group in sample_groups) {
      if (count == 1 | count %% 500 == 0) {
        print(paste0("Building regressions... (", count, "/", tot_num_regr, ")"))
      }
      
      samples <- comb_ptm_prot %>% filter(id.x == ptm_site)
      if (!is.na(group)) {
        samples <- samples %>% filter(!!rlang::sym(groups_colname) == group)
      }

      result <- samples[ , c("id.x", "id.y")]
      if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_n_values) {
        lm_design <- model.matrix(lm_formula, samples)
        model <- lmFit(samples$value, lm_design, method = lm_method)
        
        fitted_ptm_levels <- unlist(as.list(fitted(model)))
        result$residuals <- samples$value - fitted_ptm_levels
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

