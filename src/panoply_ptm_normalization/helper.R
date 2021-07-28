#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Normalizes ptm to protein level by fitting a global linear model and returning residuals
# Harry Kane, Pierre Beltran, D. R. Mani, Karsten Krug, Zachery Gillette and Surya Mani
#
# Helper functions for processing PTM and protein data for normalization
#
#########################################################################################################


library("pacman")
p_load(cmapR)
p_load(reshape)
p_load(yaml)
p_load(tidyr)


match_ptm_to_proteome <- function(ptm,
                                  proteome,
                                  accession_number = "accession_number",
                                  accession_numbers = "accession_numbers",
                                  try.all.accession.numbers = FALSE) {
  # Spectrum Mill specific (extensible to other search engines if function arguments are set appropriately):
  # if PTM accession number does not have match in the proteome, will try to match all accession numbers
  # for that site in the proteome. if there are multiple protein matches, it will pick
  # the highest scoring one (if score column present) or the first one (if score column absent).
  # if this does not result in a match, search is extended to include all proteome accession numbers.
  # replaces value in the 'accession_number' column of the ptm GCT with the new
  # match (duplicates original accession number column so those values are saved).
  if (try.all.accession.numbers && !is.null(accession_numbers)) {
    original_accession_number <- ptm@rdesc[,accession_number]
    ptm@rdesc <- data.frame(ptm@rdesc, original_accession_number, stringsAsFactors = FALSE)
    ptm@rdesc <- swap_accession_numbers(ptm@rdesc, proteome@rdesc, accession_number,
                                        accession_numbers, accession_sep, score)
  }
  
  return(ptm)
}

swap_accession_numbers <- function(ptm.rdesc, prot.rdesc, accession_number,
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

merge_ptm_prot_df <- function(ptm, proteome, accession_number = "accession_number") {
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

assign_treatment_label <- function(sample_name, sample_groups = c("06h", "24h")) {
  count <- 1
  for (group in sample_groups) {
    if (grepl(group, sample_name) == 1) {
      return(as.factor(count))
    }
    count <- count + 1
  }
}
