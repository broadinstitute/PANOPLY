

summarize.cmap.results <- function (subset.scores.dir, results.prefix, group, dtype, nperm, 
                                    permuted.scores.dir, parallel=TRUE, ...) {
  # determine connectivity scores for (i) actual genesets; and (ii) permuted genesets
  # output results for actual genesets, and determine FDR using permuted genesets
  
  ## All paths, libraries/source files and support functions are in this function to enable easier use of parallelism using %dopar%
  ## Path for R-utilites (for I/O and other misc functions)

  
  ###
  ### Main function for determining connectivity
  ###
  cmap.connectivity <- function (subset.scores.dir, results.prefix, group, dtype, permutation=NULL,
                                 rankpt.n=4, mean.rankpt.threshold=85, cis.fdr=0.05) {
    ## calculate CMAP connectivity scores after assembling subset ssGSEA (un-normalized ES) scores into a single table
    ## determine significant genes: 
    ##   a gene has a significant effect if it is CIS-enriched
    ##   enrichment is defined as mean.rankpt<n> > mean.rankpt.threshold
    
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
    
    
    ES <- NULL
    filename.pattern <- '-scores\\.gct$'
    if (!is.null (permutation)) filename.pattern <- sprintf ("-%d%s", permutation, filename.pattern)
    for (f in dir (subset.scores.dir, pattern=filename.pattern)) {
      d <- parse.gctx (file.path (subset.scores.dir, f))
      if (is.null (ES)) ES <- d
      else ES <- add.cols.gct (ES, d@mat, d@cdesc)
    }
    es <- data.frame (t (ES@mat))
    
    
    write.score <- function (d, f) {
      # utility function
      # write out d to file f filling in annotations using ES
      if (!is.null (permutation)) return ()    # do not write for permutation scores
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
    get.rankpt <- function (rx, cx) {
      # get rankpt value from row rx, column cx; return 0 is NA or not present
      ret.val <- ifelse (cx %in% sig.cols, 
                         {val <- mean_rankpt.n [rx, cx]; ifelse (is.na (val), 0, val)}, 
                         0)
      return (ret.val)
    }
    mean_rankpt.n <- rankpts %>% group_by (gene=ES@cdesc$pert_iname) %>% summarise_all (~ mean.rankpt (.))
    write.csv (mean_rankpt.n, sprintf ("%s-meanrankpt%d.csv", results.prefix, rankpt.n), row.names=FALSE)
    # a gene is significant if it is CIS-enriched in both CNAdel (+ve) and CNAamp (-ve),
    # ... and gene has positive, significant cis-correlation (cna vs dtype)
    # genes used for the analysis are listed in the *-gene-list.csv file output by cmap-input.r
    genes <- read.csv ( sprintf ("%s-cmap-%s-gene-list.csv", group, dtype), as.is=TRUE )[,1]
    cis.sigtable <- read.csv ( sprintf ("%s-%s-cis-correlation.csv", group, dtype), row.names=1 )
    sig.amp <- sig.del <- sig.both <- NULL
    sig.cols <- colnames (mean_rankpt.n)
    for (g in genes) {
      if (cis.sigtable [g, 'adj.pvalue'] > cis.fdr || cis.sigtable [g, 'correlation'] < 0) next
      print (g)
      g.row <- which (mean_rankpt.n [,'gene'] == g)
      if (length (g.row) != 1) next   # gene not found in mean_rankpt.n table
      amp <- get.rankpt (g.row, sprintf ("%s_CNAamp", g))
      del <- get.rankpt (g.row, sprintf ("%s_CNAdel", g))
      if (amp <= -mean.rankpt.threshold && del >= mean.rankpt.threshold) { sig.both <- c (sig.both, g) }
      else {
        if (amp <= -mean.rankpt.threshold) sig.amp <- c (sig.amp, g)
        if (del >= mean.rankpt.threshold) sig.del <- c (sig.del, g)
      }
    }
    write.table (sig.both, sprintf ("%s-sig-genes-bidirectional.txt", results.prefix))
    write.table (c (sig.amp, sig.del), sprintf ("%s-sig-genes-unidirectional.txt", results.prefix))
    
    return (list (amp.and.del=sig.both, amp.or.del=c(sig.amp, sig.del)))
  }
  
  ###
  ### support functions
  ###
  poisson.score.ci <- function (x, alpha=0.05) {
    # the random permuations are treated as poisson trials with false positive rate lambda
    # n * lambda is small -- hence Wald CI is not appropriate
    # hence use score CI
    # see "A Comparison of Nine Confidence Intervals for a Poisson Parameter When the Expected Number
    #      of Events Is <=5", Lawrence Barker, The American Statistician, Vol. 56, No. 2, pp. 85-89, 2002.
    x.bar <- mean (x)
    n <- length(x)
    z <- qnorm (alpha/2, lower.tail=FALSE)
    x.est <- x.bar + z^2 / (2*n)
    ci.factor <- z * (4*x.bar + z^2/n)^0.5 / (4*n)^0.5
    return (list (estimate=x.est, ci=c (x.est-ci.factor, x.est+ci.factor)))
  }
  
  
  print.FDR.summary <- function (act, perm, type, output.file) {
    ## write out FDR summary
    str <- '\n   '
    act.n <- ifelse (type=='and', length (act$amp.and.del), length (act$amp.or.del))
    cat ('Genes CIS-enriched in CNAAMP', type, 'CNADEL:', act.n, '\n', file=output.file, append=TRUE)
    cat (' Gene list:', str, file=output.file, append=TRUE)
    if (type=='and') cat ( paste (act$amp.and.del, collapse=str), '\n', file=output.file, append=TRUE)
    else cat ( paste (act$amp.or.del, collapse=str), '\n', file=output.file, append=TRUE)
    cat ('FDR:', min (round (perm$estimate / act.n, digits=4), 1), '[',
         min (round (perm$ci[1]/act.n, digits=4), 1), ',',
         min (round (perm$ci[2]/act.n, digits=4), 1), ']\n\n',
         file=output.file, append=TRUE)
  }
  
  
  ### Main body of summarize.cmap.results
  actual <- cmap.connectivity (scores.dir, paste (group, 'cmap', dtype, sep='-'), group, dtype, ...)
  
  if (parallel) {
    # process permutations in parallel
    pacman::p_load(doParallel)
    pacman::p_load(foreach)
    
    cl <- makeCluster (detectCores() - 1)
    registerDoParallel (cl)
    permuted <- foreach (i = 0:(nperm-1), .combine='c') %dopar% {
      p <- cmap.connectivity (permuted.scores.dir, paste (group, 'cmap', dtype, 'permutation', i, sep='-'), group, dtype, permutation=i, ...)   
      list (list (amp.and.del=p$amp.and.del, amp.or.del=p$amp.or.del,
                  N.and=length(p$amp.and.del), N.or=length(p$amp.or.del)))
    }
    on.exit (stopCluster (cl))
  } else {
    permuted <- NULL
    for (i in 0:(nperm-1)) {
      p <- cmap.connectivity (permuted.scores.dir, paste (group, 'cmap', dtype, 'permutation', i, sep='-'), group, dtype, permutation=i, ...)   
      permuted <- c (permuted, list (list (amp.and.del=p$amp.and.del, amp.or.del=p$amp.or.del,
                                           N.and=length(p$amp.and.del), N.or=length(p$amp.or.del))))
    }
  }
  # FDR
  perm.and <- poisson.score.ci (sapply (permuted, function (x) x$N.and))
  perm.or <- poisson.score.ci (sapply (permuted, function (x) x$N.or))
  
  ## write out FDR calc summary
  out.file <- sprintf ("%s-sig-genes-with-fdr.txt", results.prefix)
  if (file.exists (out.file)) file.remove (out.file)  # since all results are appended to file
  print.FDR.summary (actual, perm.and, 'and', out.file)
  print.FDR.summary (actual, perm.or, 'or', out.file)
}



## read argument options from command line
args <- commandArgs (TRUE)

scores.dir <- toString (args[1])
gr <- toString (args[2])
typ <- toString (args[3])
n <- as.integer (args[4])
perm.dir <- ifelse (n > 0, toString (args[5]), "tmp")
rpt <- as.integer (args[6])
rpt.thresh <- as.integer (args[7])


# calculate connectivity scores
summarize.cmap.results (scores.dir, paste (gr, 'cmap', typ, sep='-'), gr, typ, n, perm.dir, 
                        rankpt.n=rpt, mean.rankpt.threshold=rpt.thresh)
