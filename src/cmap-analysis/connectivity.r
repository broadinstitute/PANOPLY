
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


cmap.connectivity <- function (subset.scores.dir, results.prefix, rankpt.n=4, mean.rankpt.threshold=90) {
  ## calculate CMAP connectivity scores after assembling subset ssGSEA (un-normalized ES) scores into a single table
  ## determine significant genes: 
  ##   a gene has a significant effect if it is CIS-enriched
  ##   enrichment is defined as mean.rankpt<n> > mean.rankpt.threshold

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
  # (for normalized scores and connectivity scores calculated below; rankpt's are similar to
  #  perturbagen centric measure of connectivity)
  
  # normalized ES (p,c,s) = ES (p,c,s) / mean ( ES(p,c,s) ) where
  #  p=perturbation type
  #  c=cell line
  #  s=sign of ES
  # i.e., scores are normalized within each perturbation type and cell line,  
  # separately for +ve and -ve scores
  # NB: in dplyr, a formula is converted to a function; for a single argument to the function, use .
  #     group_by groups rows, and mutate_all transforms each column (by row groups)
  nes <- es %>% group_by (ES@cdesc$pert_type, ES@cdesc$cell_id)  %>% 
    mutate_all (~ ifelse(. > 0, . / mean (.[.>0], na.rm=T), . / mean (abs (.[.<=0]), na.rm=T)))
  nes <- nes [, 1:ncol(es)]  # remove grouping columns at the end (added by dplyr)
  NES <- data.frame (t (nes))
  colnames (NES) <- colnames (ES@mat)
  write.score (NES, sprintf ("%s-nes.gct", results.prefix))
  
  # connectivity map score (tau) -- not used in calculating rankpt's (ala old scoring method)
  # standarized NES for each signature (ie., sample or id), using current set of queries as compendium
  ## tau <- NES %>% mutate_all (~ sign (.) * 100 / length(.) * sapply (., function (x) sum (abs(.) < abs(x))) )

  # connectivity map rank-point scores -- score percentile NES over all perturbations/cell-lines for each query (== gene set)
  rankpts <- nes %>% mutate_all (~ sign (.) * 100 / length(.) * sapply (., function (x) sum (abs(.) < abs(x))) )
  rankpts.t <- t (rankpts)
  colnames (rankpts.t) <- colnames (ES@mat)
  write.score (rankpts.t, sprintf ("%s-connectivity.gct", results.prefix))
  
  # determine mean rankpt_n -- mean of (best of +ve / -ve) n rank points
  mean.rankpt <- function (x, n=rankpt.n) {
    sx <- sort (x)
    m.lo <- mean ( sx[1:n], na.rm=TRUE )
    m.hi <- mean ( rev(sx)[1:n], na.rm=TRUE )
    ifelse (abs(m.hi) >= abs(m.lo), m.hi, m.lo)
  }
  mean_rankpt.n <- rankpts %>% group_by (gene=ES@cdesc$pert_iname) %>% summarise_all (~ mean.rankpt (.))
  write.csv (mean_rankpt.n, sprintf ("%s-meanrankpt%d.csv", results.prefix, rankpt.n), row.names=FALSE)
  # a gene is significant if it is CIS-enriched in both CNAdel (+ve) and CNAamp (-ve)
  # genes used for the analysis are listed in the *-gene-list.csv file output by cmap-input.r
  genes <- read.csv ( sprintf ("%s-gene-list.csv", results.prefix), as.is=TRUE )[,1]
  sig.amp <- sig.del <- sig.both <- NULL
  sig.cols <- colnames (mean_rankpt.n)
  for (g in genes) {
    print (g)
    g.row <- which (mean_rankpt.n [,'gene'] == g)
    if (length (g.row) != 1) next   # gene not found in mean_rankpt.n table
    amp.col <- sprintf ("%s_CNAamp", g)
    amp <- ifelse (amp.col %in% sig.cols, mean_rankpt.n [g.row, amp.col], 0)
    del.col <- sprintf ("%s_CNAdel", g)
    del <- ifelse (del.col %in% sig.cols, mean_rankpt.n [mean_rankpt.n [,'gene'] == g, del.col], 0)
    if (amp <= -mean.rankpt.threshold && del >= mean.rankpt.threshold) { sig.both <- c (sig.both, g) }
    else {
      if (amp <= -mean.rankpt.threshold) sig.amp <- c (sig.amp, g)
      if (del >= mean.rankpt.threshold) sig.del <- c (sig.del, g)
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
