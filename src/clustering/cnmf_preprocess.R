# Filename: CNMF_Preprocess.R
#  Authors: Pablo Tamayo (tamayo@genome.wi.mit.edu), R. Zupko II (rzupko@broadinstitute.org)
#
#  Purpose: This file is used by GDAC_CNMF.R to preprocess the datasets that are provided.
options(warn = -1)

GSEA.NormalizeCols <- function(V) {
  # Standardize columns of a gene expression matrix
  col.mean <- apply(V, MARGIN = 2, FUN = mean)
  col.sd   <- apply(V, MARGIN = 2, FUN = sd)
  for (i in 1:ncol(V)) {
    if (col.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[,i] <- (V[,i] - col.mean[i]) / col.sd[i]
    }
  }
  return(V)
}

GSEA.Threshold <- function(V, thres, ceil) {
  # Threshold and ceiling pre-processing for gene expression matrix
  V[V < thres] <- thres
  V[V > ceil]  <- ceil
  return(V)
}


MSIG.Gct2Frame <- function(filename = "NULL") {
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  ds <- read.delim(filename, header = T, sep = "\t", skip = 2, blank.lines.skip = T, comment.char = "", as.is = T, na.strings = "")
  t  <- which(is.na(ds[,1]))
  if (length(t) > 0) {
    ds <- ds[-t,]
  }
  rownames(ds) <- ds[,1]
  descs <- ds[,1]
  ds    <- ds[,-c(1:2)]
  row.names <- row.names(ds)
  names <- names(ds)
  return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

MSIG.NormalizeCols.Rank <- function(V) {
  # Column rank normalization
  cols <- length(V[1,])
  rows <- length(V[,1])
  for (j in 1:cols) {  
    V[,j] <- rank(V[,j], ties.method = "average")
  }
  return(V)
}

MSIG.NormalizeCols.Rescale <- function(V) {
  # Column rank normalization
  epsilon <- 0.00001
  for (j in 1:ncol(V)) {
    max.v <- max(V[,j])
    min.v <- min(V[,j])
    V[,j] <- (V[,j] - min.v + epsilon) / (max.v - min.v)
  }
  return(V)
}

MSIG.Preprocess.Dataset <- function(input.ds, output.ds, thres = NULL, ceil = NULL, shift = NULL,
                                    fold = NULL, delta = NULL, normalization = NULL, cntrl.genes = NULL) {
  # Read dataset
  dataset <- MSIG.Gct2Frame(filename = input.ds)
  m <- data.matrix(dataset$ds)
  gs.names <- dataset$row.names
  gs.descs <- dataset$descs
  sample.names <- dataset$names
  # threshold, ceiling and shift
  if (!is.null(thres)) {
    m[m < thres] <- thres
  }
  if (!is.null(ceil)) {
    m[m > ceil] <- ceil
  }
  if (!is.null(shift)) {
    m <- m + shift
  }
  # identify and save control genes
  if (!is.null(cntrl.genes)) {
    gene.names2 <- intersect(cntrl.genes, gs.names)
    locs        <- match(gene.names2, gs.names, nomatch=0)
    msig.cntrl  <- m[locs, ]
    msig.cntrl.genes <- gs.names[locs]
    msig.cntrl.descs <- gs.descs[locs]
    m <- m[-locs, ]
    gs.names <- gs.names[-locs]
    gs.descs <- gs.descs[-locs]
  }
  # variation filter
  if ((!is.null(fold)) && (!is.null(delta))) {
    temp <- MSIG.VarFilter(V = m, fold = fold, delta = delta, gene.names = gs.names, gene.descs = gs.descs)
    m <- temp$V
    gs.names <- temp$new.gene.names
    gs.descs <- temp$new.gene.descs
    dim(m)
  }
  # restore control genes
  if (!is.null(cntrl.genes)) {
    m <- rbind(m, msig.cntrl)
    gs.names <- c(gs.names, msig.cntrl.genes)
    gs.descs <- c(gs.descs, msig.cntrl.descs)
  }
  # normalization
  if (!is.null(normalization)) {
    if (normalization == 1) {
      m <- MSIG.NormalizeCols.Rank(m)
    } else if (normalization == 2) {
      m <- MSIG.NormalizeCols.Rank(m)/length(m[,1])
    } else if (normalization == 3) {
      m <- GSEA.NormalizeCols(m) + 3
      m <- GSEA.Threshold(m, 0.001, 100000)
    } else if (normalization == 4) {
      m <- MSIG.NormalizeCols.Rank(m)/length(m[,1])
    } else if (normalization == 5) {
      m <- MSIG.NormalizeCols.Rescale(m)
    } else if (normalization == 6) {
      cols <- length(m[1,])
      # Column rank normalization from 0 to N - 1
      for (j in 1:cols) {  
        m[,j] <- rank(m[,j], ties.method = "average") - 1
      }
      m <- 10000*m/(length(m[,1]) - 1)
      } else if (normalization == 7) {
        m <- ((100*MSIG.NormalizeCols.Rank(m))%/%length(m[,1]) + 1)
      } else if (normalization == 8) {
        row.mean <- apply(m, MARGIN=1, FUN=mean)
        for (i in 1:length(m[,1])) {
          m[i,] <- m[i,] / row.mean[i]
        }
      }
    }
    V <- data.frame(m)
    names(V) <- sample.names
    row.names(V) <- gs.names
    write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
}

MSIG.VarFilter <- function(V, fold, delta, gene.names = "", gene.descs = "") {
  # Variation filter pre-processing for gene expression matrix
  cols <- length(V[1,])
  rows <- length(V[,1])

  row.max <- apply(V, MARGIN=1, FUN=max)
  row.min <- apply(V, MARGIN=1, FUN=min)

  flag <- array(dim=rows)
  flag <- (row.max /row.min >= fold) & (row.max - row.min >= delta)
  size <- sum(flag)
  B <- matrix(0, nrow = size, ncol = cols)
  j <- 1
  if (length(gene.names) == 1) {
    for (i in 1:rows) {
      if (flag[i]) {
        B[j,] <- V[i,]
        j <- j + 1
      }
    }
    return(B)
  }
  new.gene.names <- vector(mode = "character", length = size)
  new.gene.descs <- vector(mode = "character", length = size)
  for (i in 1:rows) {
    if (flag[i]) {
      B[j,] <- V[i,]
      new.gene.names[j] <- gene.names[i]
      new.gene.descs[j] <- gene.descs[i]
      j <- j + 1
    }
  }
  return(list(V = B, new.gene.names = new.gene.names, new.gene.descs = new.gene.descs, locations = flag))
}


write.gct <- function(gct.data.frame, descs = "", filename) {
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  names <- names(gct.data.frame)
  cat("\t", names[1], file = f, append = TRUE, sep = "")
  if (length(names) > 1) {
    for (j in 2:length(names)) {
      cat("\t", names[j], file = f, append = TRUE, sep = "")
    }
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
  m[, 1] <- row.names(gct.data.frame)
  if (length(descs) > 1) {
    m[, 2] <- descs
  } else {
    m[, 2] <- row.names(gct.data.frame)
  }
  index <- 3
  for (i in 1:dim(gct.data.frame)[2]) {
    m[, index] <- gct.data.frame[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
}
