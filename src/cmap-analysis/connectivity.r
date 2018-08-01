
## Path for R-utilites (for I/O and other misc functions)
Rutil.path <- switch (Sys.info()[['sysname']],
                      Windows = {'//flynn-cifs/prot_proteomics/Projects/R-utilities'},
                      Darwin = {'/Volumes/prot_proteomics/Projects/R-utilities'},
                      Linux = {'/prot/proteomics/Projects/R-utilities'})

Source <- function (f) {
  # adaptive source file location -- from current directory or from R-utiliites
  if (file.exists (f)) {source (f)}
  else {
    f.os <- file.path (Rutil.path, f)
    if (file.exists (f.os)) source (f.os)
    else stop (paste ("Can't find file", f.os, '-- R-utilities missing?'))
  }
}


Source ('gct-io.r')

pacman::p_load (dplyr)


cmap.connectivity <- function (subset.scores.dir, results.prefix, tau.threshold=90, threshold.count=4) {
  ## calculate CMAP connectivity scores after assembling subset ssGSEA (un-normalized ES) scores into a single table
  ## determine significant genes: 
  ##   a gene has a significant effect if it is CIS-enriched
  ##   enrichment is defined as at least threshold.count profiles having tau > tau.threshold
  #     (defaults for tau.threshold and threshold.count are set to mimic the mean_rankpt_4 metric)
  
  ES <- NULL
  for (f in dir (subset.scores.dir, pattern='-scores\\.gct$')) {
    d <- parse.gctx (file.path (subset.scores.dir, f))
    if (is.null (ES)) ES <- d
    else ES <- add.cols.gct (ES, d@mat, d@cdesc)
  }
  es <- data.frame (t (ES@mat))
  
  
  write.score <- function (d, f) {
    # utility function
    # write out d to file f filling in annotations using ES
    D <- ES
    D@mat <- as.matrix (d)
    write.gct (D, f)
  }
  
  
  # See Subramanian et. al., Cell 171, 1437-1452. 2017.
  # (for normalized scores and connectivity scores calculated below)
  
  # normalized ES (p,c,s) = ES (p,c,s) / mean ( ES(p,c,s) ) where
  #  p=perturbation type
  #  c=cell line
  #  s=sign of ES
  # i.e., scores are normalized within each perturbation type and cell line,  
  # separately for +ve and -ve scores
  # NB: in dplyr, a formula is converted to a function; for a single argument to the function, use .
  nes <- es %>% group_by (ES@cdesc$pert_type, ES@cdesc$cell_id)  %>% 
    mutate_all (~ ifelse(. > 0, . / mean (.[.>0], na.rm=T), . / mean (abs (.[.<=0]), na.rm=T)))
  nes <- nes [, 1:ncol(es)]  # remove grouping columns at the end (added by dplyr)
  NES <- data.frame (t (nes))
  colnames (NES) <- colnames (ES@mat)
  write.score (NES, sprintf ("%s-nes.gct", results.prefix))
  
  # connectivity map score (tau)
  # standarized NES for each signature (ie., sample or id), using current set of queries as compendium
  tau <- nes %>% mutate_all (~ sign (.) * 100 / length(.) * sapply (., function (x) sum (abs(.) < abs(x))) )
  tau.t <- t (tau)
  colnames (tau.t) <- colnames (ES@mat)
  write.score (tau.t, sprintf ("%s-connectivity.gct", results.prefix))
  
  # summarize TAU scores for each gene by counting the number of tau's > tau.threshold
  # take into account +ve and -ve extremes, and penalize for mixed extremes in a given gene
  sig.counts <- tau %>% group_by (gene=ES@cdesc$pert_iname) %>% summarise_all (~ sum (. > tau.threshold) - sum (. < -tau.threshold))
  write.csv (sig.counts, sprintf ("%s-sigcounts.csv", results.prefix), row.names=FALSE)
  # a gene is significant if it is CIS-enriched in both CNAdel (+ve) and CNAamp (-ve)
  # genes uses for the analysis are listed in the *-gene-list.csv file output by cmap-input.r
  genes <- read.csv ( sprintf ("%s-gene-list.csv", results.prefix), as.is=TRUE )[,1]
  sig.amp <- sig.del <- sig.both <- NULL
  sig.cols <- colnames (sig.counts)
  for (g in genes) {
    print (g)
    g.row <- which (sig.counts [,'gene'] == g)
    if (length (g.row) != 1) next   # gene not found in sig.counts table
    amp.col <- sprintf ("%s_CNAamp", g)
    amp <- ifelse (amp.col %in% sig.cols, sig.counts [g.row, amp.col], 0)
    del.col <- sprintf ("%s_CNAdel", g)
    del <- ifelse (del.col %in% sig.cols, sig.counts [sig.counts [,'gene'] == g, del.col], 0)
    if (amp <= -threshold.count && del >= threshold.count) { sig.both <- c (sig.both, g) }
    else {
      if (amp <= -threshold.count) sig.amp <- c (sig.amp, g)
      if (del >= threshold.count) sig.del <- c (sig.del, g)
    }
  }
  write.table (sig.both, sprintf ("%s-sig-genes-bidirectional.txt", results.prefix))
  write.table (c (sig.amp, sig.del), sprintf ("%s-sig-genes-unidirectional.txt", results.prefix))
}



## read argument options from command line
args <- commandArgs (TRUE)

scores.dir <- toString (args[1])
gr <- toString (args[2])
typ <- toString (args[3])


# calculate connectivity scores
cmap.connectivity (scores.dir, paste (gr, 'cmap', typ, sep='-'))
