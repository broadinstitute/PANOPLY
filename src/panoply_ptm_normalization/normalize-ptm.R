#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Normalizes ptm to protein level by fitting a global linear model and returning residuals
# Harry Kane, Pierre Beltran, D. R. Mani, Karsten Krug, Zachery Gillette and Surya Mani
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
p_load(MASS)  # rlm
p_load(limma)

source("panoply_ptm_normalization/helper.R")
source("panoply_ptm_normalization/metrics.R")


normalize_ptm <- function(proteome.gct, ptm.gct, output.prefix = NULL, 
                          try.all.accession.numbers = TRUE,        # try hard to find a match (using accession_numbers)
                          accession_number = "accession_number",   # column with protein/ptm accession number
                          accession_numbers = "accession_numbers", # accession_numbers for protein/ptm group
                          accession_sep = "|",                     # separator for each accession number in accession_numbers
                          score = "scoreUnique",                   # column with protein scores
                          ndigits = 5,
                          norm_method = "global",                       # normalization method (options: global, pairwise)
                          lm_method = FALSE,
                          use_all_samples = TRUE,                  # pairwise: use all samples to build linear model for each PTM-protein pair
                          sample_groups = c("06h", "24h", "96h"),  # pairwise: if not use all samples, what subsamples to build lm for? 
                          sample_groups_col = "pert_time",
                          groups_as_covariate = TRUE,              # pairwise: use sample groups as a covariate and build one regression for samples across groups
                          min_values_present = 4)                  # pairwise: what least number of samples must contain both non-NA PTM and protein values
{
  # import GCT files
  proteome <- parse_gctx(proteome.gct)
  temp_ptm <- parse_gctx(ptm.gct)
  ptm <- match_ptm_to_proteome(temp_ptm, proteome, accession_number, accession_numbers, try.all.accession.numbers)
  
  # fits linear model and returns updated GCT
  if (norm_method == "global") {
    ptm.norm <- normalize_global(ptm, proteome, accession_number, lm_method)
  } else if (norm_method == "pairwise") {
    ptm.norm <- normalize_pairwise(ptm, proteome, accession_number, lm_method, use_all_samples, sample_groups, sample_groups_col, groups_as_covariate, min_values_present)
  }
  
  # writes and returns updated GCT
  file.prefix <- ifelse (!is.null (output.prefix), output.prefix,
                         unlist(strsplit(ptm.gct, split = '.gct', fixed = TRUE))[1])
  
  write_gct (ptm.norm, paste(file.prefix, '-proteome-relative-norm.gct', sep = ''), 
             appenddim = FALSE, precision=ndigits)
  
  invisible (ptm.norm)
}


#' Applies linear regression to correct PTM levels for underlying protein levels
normalize_global <- function (ptm, proteome, accession_number, lm_method) {
  # warn if not all samples are matched
  matched.samples <- intersect (ptm@cid, proteome@cid)
  if (length(matched.samples) != length(ptm@cid)) {
    warning ('WARNING: not all samples in ptm file have matches in proteome ... unmatched samples removed.')
  }
  
  data <- merge_ptm_prot_df(ptm, proteome, accession_number)
  
  # fit global model
  print ("Fitting model...") 
  if (lm_method == "robust") {
    model <- rlm(value ~ value.prot, data = data, na.action = na.exclude)
  } else {
    model <- lm(value ~ value.prot, data = data, na.action = na.exclude)
  }
  residuals <- residuals(model)
  results <- data.frame (data$id.x, data$id.y, residuals)
  colnames(results) <- c('id.x', 'id.y', 'residuals')
  print ("Success.")
  print (summary(model))
  
  ptm <- update_gct(ptm, results)
  return(ptm)
}


#' Applies linear regression to correct PTM levels for underlying protein levels
#' for each PTM-protein pair across selected or all samples
normalize_pairwise <- function(ptm, proteome, accession_number, lm_method, use_all_samples, sample_groups, sample_groups_col, groups_as_covariate, min_values_present) {
  # warn if not all samples are matched
  matched.samples <- intersect (ptm@cid, proteome@cid)
  if (length(matched.samples) != length(ptm@cid)) {
    warning ('WARNING: not all samples in ptm file have matches in proteome ... unmatched samples removed.')
  }
  
  combined_data <- merge_ptm_prot_df(ptm, proteome, accession_number)
  combined_data <- combined_data %>% filter(grepl(paste(sample_groups, collapse="|"), id.y))  # filter only groups of interest
  all_ptms <- unique(combined_data$id.x)
  sample_names <- combined_data$id.y
  if (use_all_samples | groups_as_covariate) {
    sample_groups = c("")  # just so that it matches to everything
  }
  combined_data$treatment_group <- combined_data[[sample_groups_col]]  # create a variable indicating group
  combined_data$treatment_group <- as.factor(combined_data$treatment_group)
  
  tot_num_regr <- length(all_ptms) * length(sample_groups)
  print(paste("Found", length(all_ptms), "PTM-protein pairs and", length(sample_groups),
              "sample group(s);", tot_num_regr, "regressions to be built", sep = " "))
  
  # create a dataframe placeholder to store results in
  all_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all_results) <- c("id.x", "id.y", "residuals")
  
  # loop through all PTMs, build regression for each
  count <- 1
  count_fail <- 0
  for (ptm_site in all_ptms) {
    for (group in sample_groups) {
      if (count == 1 | count %% 500 == 0) {
        print(paste0("Building regressions... (", count, "/", tot_num_regr, ")"))
      }
      
      group_data <- combined_data[grep(group, sample_names)]
      samples <- group_data[group_data$id.x == ptm_site]
      
      if (sum(!(is.na(samples$value) | is.na(samples$value.prot))) >= min_values_present) {
        if (groups_as_covariate) {
          lm_design <- model.matrix(~ value.prot + treatment_group, data = samples)
        } else {
          lm_design <- model.matrix(~ value.prot, data = samples)
        }
        model <- lmFit(samples$value, lm_design, method = lm_method)
        
        result <- samples[ , c("id.x", "id.y")]
        fitted_ptm_levels <- unlist(as.list(fitted(model)))
        result$residuals <- samples$value - fitted_ptm_levels
        all_results <- rbind(all_results, result)
      } else {
        count_fail <- count_fail + 1
      }
      count <- count + 1
    }
  }
  
  print(paste0("Regressions failed (values present < ", min_values_present, "): ", count_fail))
  print("Done.")
  
  ptm <- update_gct(ptm, all_results)
  return(ptm)
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

