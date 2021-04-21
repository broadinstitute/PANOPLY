#
# Copyright (c) 2021 The Broad Institute, Inc. All rights reserved.
#

#########################################################################################################
#
# Normalizes PTM to protein level by fitting a global linear model and returning residuals
# Harry Kane, Pierre Beltran, D. R. Mani, Karsten Krug, Zachery Gillette and Surya Mani
#
# Fits PTM = beta_0 + beta_1*protein to all matched points in a dataset, and returns residuals as
# protein-corrected PTM values.
#
#########################################################################################################


library('pacman')
p_load(cmapR)
p_load(reshape)
p_load(yaml)
p_load(tidyr)


normalize_ptm <- function (proteome.gct, ptm.gct, output.prefix=NULL, 
                           try.all.accession.numbers=TRUE,        # try hard to find a match (using accession_numbers)
                           accession_number='accession_number',   # column with protein/PTM accession number
                           accession_numbers='accession_numbers', # accession_numbers for protein/PTM group
                           accession_sep='|',                     # spearator for each accession number in accession_numbers
                           score='scoreUnique',                   # column with protein scores
                           ndigits=5)
{
  # import GCT files
  proteome <- parse.gctx (proteome.gct)
  PTM <- parse.gctx (ptm.gct)

  # Spectrum Mill specific (extensible to other search engines if function arguments are set appropriately):
  # if PTM accession number does not have match in the proteome, will try to match all accession numbers
  # for that site in the proteome. if there are multiple protein matches, it will pick
  # the highest scoring one (if score column present) or the first one (if score column absent).
  # if this does not result in a match, search is extended to include all proteome accession numbers.
  # replaces value in the 'accession_number' column of the PTM GCT with the new
  # match (duplicates original accession number column so those values are saved).
  if (try.all.accession.numbers && !is.null (accession_numbers)) {
    original_accession_number <- PTM@rdesc[,accession_number]
    PTM@rdesc <- data.frame (PTM@rdesc, original_accession_number, stringsAsFactors = FALSE)
    PTM@rdesc <- swap.accession.numbers (PTM@rdesc, proteome@rdesc, accession_number,
                                         accession_numbers, accession_sep, score)
  }
  
  
  # fits linear model and returns updated GCT
  PTM.norm <- normalize(PTM, proteome, accession_number)
  
  
  # writes and returns updated GCT
  file.prefix <- ifelse (!is.null (output.prefix), output.prefix,
                         unlist(strsplit(PTM.gct, split = '.gct', fixed = TRUE))[1])
  
  write.gct (PTM.norm, paste(file.prefix, '-proteome-relative-norm.gct', sep = ''), 
             appenddim = FALSE, precision=ndigits)
  
  invisible (PTM.norm)
}



swap.accession.numbers <- function (PTM.rdesc, prot.rdesc, accession_number,
                                    accession_numbers, accession_sep, score)
{
  # Matches PTM and proteome by accession number, and tries all possible PTM accession numbers if primary
  # accession number doesn't have proteome match
  prot.rdesc.expanded <- prot.rdesc %>% separate_rows(!!as.name (accession_numbers), sep=accession_sep)
  
  for (row in c(1:nrow(PTM.rdesc))) {
    match <- which (prot.rdesc[,accession_number] == PTM.rdesc[row, accession_number])
    if(length(match) == 0){
      # first try to match PTM accession_numbers to protein primary ID
      accession.numbers <- unlist(strsplit(PTM.rdesc[row,accession_numbers], accession_sep))
      matches <- prot.rdesc[which(unlist (prot.rdesc[,accession_number]) %in% accession.numbers),]
      
      if (nrow(matches) >= 1) {
        best_index <- ifelse (!is.null(score), which.max(matches[,score]), 1)
        PTM.rdesc[row, accession_number] <- matches[best_index, accession_number]
      } else {
        # if the above fails, try to match PTM accession_numbers to all proteome accession_numbers
        matches <- prot.rdesc.expanded[which(unlist(prot.rdesc.expanded[,accession_numbers]) %in% accession.numbers),]
        if (nrow(matches) >= 1) {
          best_index <- ifelse (!is.null(score), which.max(unlist (matches[,score])), 1)
          PTM.rdesc[row, accession_number] <- matches[best_index,accession_number]
        } 
      }
    }
  }
  
  return(PTM.rdesc)
}



normalize <- function (PTM, proteome, accession_number) {
  # Applies linear regression to correct PTM levels for underlying protein levels
  
  # warn if not all samples are matched
  matched.samples <- intersect (PTM@cid, proteome@cid)
  if (length(matched.samples) != length(PTM@cid)) {
    warning ('WARNING: not all samples in PTM file have matches in proteome ... unmatched samples removed.')
  }
  
  # create merged data table
  PTM.melt <- melt.gct (PTM)
  prot.melt <- melt.gct (proteome)
  prot.melt.data.only <- data.frame (prot.melt$id.y, prot.melt[,accession_number], prot.melt$value)
  colnames (prot.melt.data.only) <- c ('id.y', accession_number, 'value.prot')
  data <- merge (PTM.melt, prot.melt.data.only, by = c('id.y', accession_number))
  
  # print metrics
  percent <-round (100*nrow(data)/nrow(PTM.melt), digits = 1)
  print (paste (nrow(data), ' out of ', nrow(PTM.melt), ' PTM peptides normalized',
                ' (', percent, '%).', sep = ''))
  
  # fit global model
  print ("Fitting model...") 
  model <- lm (value ~ value.prot, data = data)
  residuals <- residuals (model)
  results <- data.frame (data$id.x, data$id.y, residuals)
  colnames(results) <- c('id.x', 'id.y', 'residuals')
  print ("Success.")
  print (summary(model))
  
  # compile results into matrix
  results.df <- data.frame (reshape::cast(results, id.x ~ id.y, value.var = residuals))
  results.mat <- data.matrix (results.df[, c(2:ncol(results.df))])
  rownames(results.mat) <- as.character (results.df$id.x)
  
  # reset GCT
  PTM@rdesc <- PTM@rdesc[match(rownames(results.mat), PTM@rdesc$id),]
  PTM@cdesc <- PTM@cdesc[match(colnames(results.mat), PTM@cdesc$id),]
  PTM@rid <- rownames(results.mat)
  PTM@cid <- colnames(results.mat)
  PTM@mat <- results.mat
  
  return(PTM)
}


if (!interactive()) {
## call via command line
## usage: Rscript normalize-ptm.R <proteome.gct> <ptm.gct> <output.prefix> <yaml.file>

# process command line args
args <- commandArgs (TRUE)
proteome_gct <- as.character (args[1])
ptm_gct <- as.character (args[2])
output_prefix <- as.character (args[3])
if (output_prefix == "NULL") output_prefix <- NULL
yaml_file <- args[4]

# read yaml parameters
yaml_params <- read_yaml (yaml_file)
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
