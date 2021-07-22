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


normalize_ptm <- function(proteome.gct, ptm.gct, output.prefix = NULL, 
                          try.all.accession.numbers = TRUE,        # try hard to find a match (using accession_numbers)
                          accession_number = "accession_number",   # column with protein/ptm accession number
                          accession_numbers = "accession_numbers", # accession_numbers for protein/ptm group
                          accession_sep = "|",                     # separator for each accession number in accession_numbers
                          score = "scoreUnique",                   # column with protein scores
                          ndigits = 5,
                          method = "global",                       # normalization method (options: global, pairwise)
                          robust = FALSE,
                          use_all_samples = TRUE,                  # pairwise: use all samples to build linear model for each PTM-protein pair
                          sample_groups = c("06h", "24h", "96h"),  # pairwise: if not use all samples, what subsamples to build lm for? 
                          min_values_present = 4)                  # pairwise: what least number of samples must contain both non-NA PTM and protein values
{
  # import GCT files
  proteome <- parse_gctx(proteome.gct)
  ptm <- parse_gctx(ptm.gct)

  # Spectrum Mill specific (extensible to other search engines if function arguments are set appropriately):
  # if PTM accession number does not have match in the proteome, will try to match all accession numbers
  # for that site in the proteome. if there are multiple protein matches, it will pick
  # the highest scoring one (if score column present) or the first one (if score column absent).
  # if this does not result in a match, search is extended to include all proteome accession numbers.
  # replaces value in the 'accession_number' column of the ptm GCT with the new
  # match (duplicates original accession number column so those values are saved).
  if (try.all.accession.numbers && !is.null (accession_numbers)) {
    original_accession_number <- ptm@rdesc[,accession_number]
    ptm@rdesc <- data.frame (ptm@rdesc, original_accession_number, stringsAsFactors = FALSE)
    ptm@rdesc <- swap.accession.numbers (ptm@rdesc, proteome@rdesc, accession_number,
                                         accession_numbers, accession_sep, score)
  }
  
  
  # fits linear model and returns updated GCT
  if (method == "global") {
    ptm.norm <- normalize_global(ptm, proteome, accession_number, robust)
  } else if (method == "pairwise") {
    ptm.norm <- normalize_pairwise(ptm, proteome, accession_number, robust, use_all_samples, sample_groups, min_values_present)
  }
  
  
  # writes and returns updated GCT
  file.prefix <- ifelse (!is.null (output.prefix), output.prefix,
                         unlist(strsplit(ptm.gct, split = '.gct', fixed = TRUE))[1])
  
  write_gct (ptm.norm, paste(file.prefix, '-proteome-relative-norm.gct', sep = ''), 
             appenddim = FALSE, precision=ndigits)
  
  invisible (ptm.norm)
}


#' Applies linear regression to correct PTM levels for underlying protein levels
normalize_global <- function (ptm, proteome, accession_number) {
  # warn if not all samples are matched
  matched.samples <- intersect (ptm@cid, proteome@cid)
  if (length(matched.samples) != length(ptm@cid)) {
    warning ('WARNING: not all samples in ptm file have matches in proteome ... unmatched samples removed.')
  }
  
  data <- merge_ptm_prot_df(ptm, proteome, accession_number)
  
  # fit global model
  print ("Fitting model...") 
  model <- lm(value ~ value.prot, data = data)
  residuals <- residuals (model)
  results <- data.frame (data$id.x, data$id.y, residuals)
  colnames(results) <- c('id.x', 'id.y', 'residuals')
  print ("Success.")
  print (summary(model))
  
  ptm <- update_gct(ptm, results)
  return(ptm)
}


#' Applies linear regression to correct PTM levels for underlying protein levels
#' for each PTM-protein pair across selected or all samples
normalize_pairwise <- function(ptm, proteome, accession_number, robust, use_all_samples, sample_groups, min_values_present) {
  # warn if not all samples are matched
  matched.samples <- intersect (ptm@cid, proteome@cid)
  if (length(matched.samples) != length(ptm@cid)) {
    warning ('WARNING: not all samples in ptm file have matches in proteome ... unmatched samples removed.')
  }
  
  combined_data <- merge_ptm_prot_df(ptm, proteome, accession_number)
  all_ptms <- unique(combined_data$id.x)
  sample_names <- combined_data$id.y
  if (use_all_samples) {
    sample_groups = c("")  # just so that it matches to everything
  }
  
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
        if (robust) {
          model <- rlm(value ~ value.prot, data = samples, na.action = na.exclude)
        } else {
          model <- lm(value ~ value.prot, data = samples, na.action = na.exclude)
        }
        
        result <- samples[ , c("id.x", "id.y")]
        result$residuals <- model$residuals
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

swap.accession.numbers <- function(ptm.rdesc, prot.rdesc, accession_number,
                                   accession_numbers, accession_sep, score)
{
  # Matches ptm and proteome by accession number, and tries all possible ptm accession numbers if primary
  # accession number doesn't have proteome match
  prot.rdesc.expanded <- prot.rdesc %>% separate_rows(!!as.name (accession_numbers), sep=accession_sep)
  
  for (row in c(1:nrow(ptm.rdesc))) {
    match <- which (prot.rdesc[,accession_number] == ptm.rdesc[row, accession_number])
    if(length(match) == 0){
      # first try to match ptm accession_numbers to protein primary ID
      matches <- prot.rdesc[which(unlist (prot.rdesc[,accession_number]) %in% accession.numbers),]
      
      if (nrow(matches) >= 1) {
        best_index <- ifelse (!is.null(score), which.max(matches[,score]), 1)
        ptm.rdesc[row, accession_number] <- matches[best_index, accession_number]
      } else {
        # if the above fails, try to match ptm accession_numbers to all proteome accession_numbers
        matches <- prot.rdesc.expanded[which(unlist(prot.rdesc.expanded[,accession_numbers]) %in% accession.numbers),]
        if (nrow(matches) >= 1) {
          best_index <- ifelse (!is.null(score), which.max(unlist (matches[,score])), 1)
          ptm.rdesc[row, accession_number] <- matches[best_index,accession_number]
        } 
      }
    }
  }
  
  return(ptm.rdesc)
}

merge_ptm_prot_df <- function(ptm, proteome, accession_number) {
  ptm.melt <- melt_gct(ptm)
  prot.melt <- melt_gct(proteome)
  prot.melt.data.only <- data.frame(prot.melt$id.y, prot.melt[ , accession_number], prot.melt$value)
  colnames(prot.melt.data.only) <- c("id.y", accession_number, "value.prot")
  data <- merge (ptm.melt, prot.melt.data.only, by = c("id.y", accession_number))
  
  # print metrics
  percent <- round(100 * nrow(data) / nrow(ptm.melt), digits = 1)
  print(paste(nrow(data), " out of ", nrow(ptm.melt), " PTM peptides to be normalized", " (", percent, "%).", sep = ""))
  
  return(data)
}

update_gct <- function(ptm, results) {
  results.mat <- result_unmelt(results)

  ptm@rdesc <- ptm@rdesc[match(rownames(results.mat), ptm@rdesc$id), ]
  ptm@cdesc <- ptm@cdesc[match(colnames(results.mat), ptm@cdesc$id), ]
  ptm@rid <- rownames(results.mat)
  ptm@cid <- colnames(results.mat)
  ptm@mat <- results.mat

  return(ptm)
}

result_unmelt <- function(res_df.melt) {
  # compile results into matrix
  results.df <- data.frame(reshape::cast(res_df.melt, id.x ~ id.y, value.var = residuals))
  results.mat <- data.matrix(results.df[, c(2:ncol(results.df))])
  rownames(results.mat) <- as.character (results.df$id.x)
  
  return(results.mat)
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
