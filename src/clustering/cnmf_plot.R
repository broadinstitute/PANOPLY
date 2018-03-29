# Filename: CNMF_Plot.R
#  Authors: Pablo Tamayo (tamayo@genome.wi.mit.edu), R. Zupko II (rzupko@broadinstitute.org)
#
#  Purpose: This file is used by CNMF_GDAC.R to plot the results of the processing.
options(warn = -1)

CNMF.ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {
  # Plots a heatmap plot of a consensus matrix
  cols <- length(V[1,])
  B <- matrix(0, nrow=cols, ncol=cols)
  max.val <- max(V)
  min.val <- min(V)
  for (i in 1:cols) {
      for (j in 1:cols) {
        k <- cols - i + 1
        B[k, j] <-  max.val - V[i, j] + min.val
      }
  }
  col.names2 <- rev(col.names)
  col.labels2 <- rev(col.labels)
  D <- matrix(0, nrow=(cols + 1), ncol=(cols + 1))
  col.tag <- vector(length=cols, mode="numeric")
  current.tag <- 0
  col.tag[1] <- current.tag
  for (i in 2:cols) {
     if (col.labels[i] != col.labels[i - 1]) {
       current.tag <- 1 - current.tag
     }
     col.tag[i] <- current.tag
  }
  col.tag2 <- rev(col.tag)
  D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
  D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
  D[(cols + 1), 1] <- 1.03
  D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]
  col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), "#BBBBBB", "#333333", "#FFFFFF") # , gamma = 1.5
  image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)
  for (i in 1:cols) {
      col.names[i]  <- paste("      ", substr(col.names[i], 1, 12), sep="")
      col.names2[i] <- paste(substr(col.names2[i], 1, 12), "     ", sep="")
  }
  axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.50, font.axis=1, line=-1)
  axis(2, at=1:cols, labels=col.labels2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)
  axis(3, at=2:(cols + 1), labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=0.50, font.axis=1, line=-1)
  axis(3, at=2:(cols + 1), labels=as.character(col.labels), adj = 1, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)
  return()
}

CNMF.matrix.abs.plot <- function(V, axes = FALSE, log = FALSE, norm = TRUE, transpose = TRUE,
                                 matrix.order = TRUE, max.v = 1, min.v = 0, main = " ", sub = " ",
                                 xlab = " ", ylab = "  ") {
  # Plot CNMF matrix as abs plot
  rows <- nrow(V)
  cols <- ncol(V)
  if (log == T) {
    V <- log(V)
  }
  B <- matrix(0, nrow=rows, ncol=cols)
  for (i in 1:rows) {
    for (j in 1:cols) {
      if (matrix.order == T) {
        k <- rows - i + 1
      } else {
        k <- i
      }
      if (norm == T) {
        if ((max.v == 1) && (min.v == 0)) {
          max.val <- max(V)
          min.val <- min(V)
        } else {
          max.val = max.v
          min.val = min.v
        }
      }
      B[k, j] <-  max.val - V[i, j] + min.val
    }
  }
  if (transpose == T) {
    B <- t(B)
  }
  if (norm == T) {
    image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), main = main, sub = sub, xlab = xlab, ylab = ylab)
  } else {
    image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9), main = main, sub = sub, xlab = xlab, ylab = ylab)
  }
  return(list(B, max.val, min.val))
}

CNMF.metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
  index <- 1:ncol(H)
  plot(index, H[1,], xlim = c(1, ncol(H)), ylim = c(min(H), max(H)), main = main, sub = sub, ylab = ylab, xlab = xlab, type = "n")
  for (i in 1:nrow(H)) {
    lines(index, H[i,], type = "l", col = i, lwd = 2)
  }
}
