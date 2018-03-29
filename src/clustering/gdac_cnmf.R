library(optparse)
library(doParallel)
library(Cairo)
options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-m", "--expfile"), action="store", type="character", help=""),
  make_option(c("-u", "--k_int"), action="store", type="character", help=""),
  make_option(c("-v", "--k_final"), action="store", type="character", help=""),
  make_option(c("-o", "--output_prefix"), action="store", type="character", help=""),
  make_option(c("-l", "--libdir"), action="store", type="character", default=getwd(), help=""),
  make_option(c("-p", "--cpu"), action="store", type="integer", default=1, help="")
)

source_files <- function(libdir) {
  # Source the files that we need
  source(file.path(libdir, "cnmf_preprocess.R"))
  source(file.path(libdir, "cnmf_data_processing.R"))
  source(file.path(libdir, "cnmf_plot.R"))
}

NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
  # Calculate the Non-negative matrix factorization of V, using the parameters and
  # constraints specified.

  # @param V - matrix
  # @param k - the current rank
  # @param maxniter - maximum number of iterations
  # @param seed - rand seed for selection of W and H
  # @param stopconv - stop convergence
  # @param stopfreq - stop frequency

  # Precalcuate the correct size
  M    <- ncol(V)
  size <- nrow(V) * M
  # Randomize and initialize W and H with random numbers
  set.seed(seed)
  W <- matrix(runif(nrow(V) * k), nrow = nrow(V), ncol = k)
  H <- matrix(runif(k * M), nrow = k, ncol = M)
  # Prepare the rest of the variables
  new.membership  <- vector(mode = "numeric", length = M)
  old.membership  <- vector(mode = "numeric", length = M)
  no.change.count <- 0
  # Run until the max iterator or there has been convergence
  for (t in 1:maxniter) {
    # Calculate the updated W and H   
    H    <- H * (t(W) %*% (V / (W %*% H))) + .Machine$double.eps
    norm <- apply(W, MARGIN = 2, FUN = sum)
    for (i in 1:k) {
       H[i,] <- H[i,] / norm[i]
    }
    W    <- W * ((V / (W %*% H)) %*% t(H)) + .Machine$double.eps
    norm <- apply(H, MARGIN = 1, FUN = sum)
    for (i in 1:k) {
      W[,i] <- W[,i] / norm[i]
    }
    # Check to see if there has been convergence
    if (t %% stopfreq == 0) {
      for (j in 1:M) {
        new.membership[j] <- order(H[,j], decreasing=T)[1]
      }
      new.membership <- as.vector(new.membership)
      if (sum(new.membership == old.membership) == M) {
        no.change.count <- no.change.count + 1
        if (no.change.count == stopconv) { 
          break
        }
      } else {
        no.change.count <- 0
      }
      old.membership <- new.membership
    }
  }
  return(list(H = H))
}

# Write the parameters for the module to the file indicated.
writeParameters <- function(filename, input.ds, k.init, k.final, num.clusterings, maxniter, error.function, 
                            rseed, directory, stopconv, stopfreq, non.interactive.run, doc.string) {
  # @param filename name of output file
  # @param input.dx

  time.string <- as.character(as.POSIXlt(Sys.time(),"GMT"))
  write(paste("Run of NMF on ", time.string), file=filename)
  write(paste("input.ds =", input.ds, sep=" "), file=filename, append=T)
  write(paste("k.init = ", k.init, sep=" "), file=filename, append=T)
  write(paste("k.final =", k.final, sep=" "), file=filename, append=T)
  write(paste("num.clusterings =", num.clusterings, sep=" "), file=filename, append=T)
  write(paste("maxniter =", maxniter, sep=" "), file=filename, append=T)
  write(paste("error.function =", error.function, sep=" "), file=filename, append=T)
  write(paste("rseed =", rseed, sep=" "), file=filename, append=T)
  write(paste("directory =", directory, sep=" "), file=filename, append=T)
  write(paste("stopconv =", stopconv, sep=" "), file=filename, append=T)
  write(paste("stopfreq =", stopfreq, sep=" "), file=filename, append=T)
  write(paste("non.interctive.run =", non.interactive.run, sep=" "), file=filename, append=T)
  write(paste("doc.string =", doc.string, sep=" "), file=filename, append=T)
}



consensusNMF <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function,
                         rseed=123456789, stopconv=40, stopfreq=10, non.interactive.run=FALSE,
                         doc.string = "", directory="", libdir = "", ...) {
  # Check to make sure the error function is supported, the only one in this
  # module is divergence
  if (error.function != "divergence") {
    stop(paste("Un-supported error function: ", error.function, sep = ""))
  }
  # Save input parameters
  filename    <- paste0(doc.string, ".params.txt")
  writeParameters(filename, input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed,
                  directory, stopconv, stopfreq, non.interactive.run, doc.string)

  # Prepare the parameters for the run
  k.init          <- as.integer(k.init)
  k.final         <- as.integer(k.final)
  num.clusterings <- as.integer(num.clusterings)
  n.iter          <- as.integer(maxniter)

  stopifnot(!is.na(rseed))

  seed <- as.integer(rseed)

  D         <- CNMF.read.dataset(input.ds)
  col.names <- names(D)
  A         <- data.matrix(D)

  # Threshold negative values to small quantity
  A[A < 0] <- .Machine$double.eps
  cols     <- ncol(A)

  num.k    <- k.final - k.init + 1
  rho      <- vector(mode = "numeric", length = num.k)
  k.vector <- vector(mode = "numeric", length = num.k)
  k.index  <- 1

  connect.matrix.ordered <- array(0, c(num.k, cols, cols))
  all.membership         <- c()

  print(paste("INFO: Cores: ", getDoParWorkers()))
  print(paste("INFO: DoMC Version: ", getDoParVersion()))

  # Run for the required number of ranks
  for (k in k.init:k.final) {
    # Loop over the clusterings and calculate the divergence
    assign <- foreach(i = 1:num.clusterings, .combine = cbind) %dopar% {
      # Find the nonnegative matrix factorization
      NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
      # Find the membership
      members <- c()
      for (j in 1:cols) {
        members <- c(members, order(NMF.out$H[,j], decreasing = T)[1])
      }
      return(members)
    }
    assign <- matrix(assign, nrow = num.clusterings, byrow = TRUE)
    # Compute consensus matrix
    connect.matrix <- matrix(0, nrow = cols, ncol = cols)
    for (i in 1:num.clusterings) {
      for (j in 1:cols) {
        for (p in 1:cols) {
          if (j == p) {
            connect.matrix[j, p] <- connect.matrix[j, p] + 1
          } else if (assign[i, j] == assign[i, p]) {
            connect.matrix[j, p] <- connect.matrix[j, p] + 1
          }
        }
      }
    }
    connect.matrix <- connect.matrix / num.clusterings
    # Compute the distance matrix
    dist.matrix <- as.dist(1 - connect.matrix)
    HC <- hclust(dist.matrix, method="average")
    # Update the rho for the index
    dist.coph <- cophenetic(HC)
    k.vector[k.index] <- k
    rho[k.index] <- signif((cor(dist.matrix, dist.coph)), digits = 4)
    # Order the matrix
    for (i in 1:cols) {
      for (j in 1:cols) {
        connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
      }
    }
    # Compute consensus clustering membership and update the k index
    all.membership <- cbind(all.membership, cutree(HC, k = k))
    k.index <- k.index + 1
  }

  # Save consensus matrices in one file
  # Plot the all clusters as a single image
  CairoPNG(filename = paste(directory, doc.string, ".", "consensus.all.k.plot.png", sep=""))
  nf <- layout(matrix(c(1:8), 2, 4, byrow=T), c(1, 1), c(1, 1, 1, 1), TRUE)
  for (k in 1:num.k) { 
    CNMF.matrix.abs.plot(connect.matrix.ordered[k,,], log=FALSE, main=paste("k=", k.vector[k]), 
                         sub=paste("Cophenetic coef.=", rho[k]), ylab="samples", xlab="samples")
  }

  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab="k", 
       ylab="Cophenetic correlation", type="n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
  dev.off() 

  # Plot each image individually
  for (k in 1:num.k) {
    CairoPNG(filename = paste(directory, doc.string, ".", "consensus.plot.k", k.vector[k], ".png", sep=""), width = 600, height = 600)
    CNMF.matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), ylab = "samples", xlab ="samples")
    dev.off()
  }

  # Plot the cophenetic coefficient
  CairoPNG(filename = paste(directory, doc.string, ".cophenetic.coefficient.png", sep = ""), width = 600, height = 800)
  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
  dev.off()

  # Prepare the results to be saved
  results.membership           <- data.frame(all.membership)
  rownames(results.membership) <- col.names
  colnames(results.membership) <- c("membership", paste("membership.", 1:(ncol(results.membership) - 1), sep = ""))

  # Write the membership matrix
  membership.filename          <- paste(directory, doc.string, ".", "membership", ".txt", sep="", collapse="")
  # rownames(results.membership) <- gsub("\\.","-", rownames(results.membership))
  write.table(results.membership, file=membership.filename, append=FALSE, quote=FALSE, sep="\t", eol="\n", 
              col.names=TRUE, row.names=TRUE)

  # Write the cophenetic coefficient for down stream use
  cophenetic.filename          <- paste(directory, doc.string, ".cophenetic.coefficient.txt", sep = "")
  rho                          <- data.frame(rho, row.names = k.init:k.final)
  write.table(rho, file=cophenetic.filename, append=FALSE, quote=FALSE, sep="\t", eol="\n", col.names=TRUE, row.names=TRUE)

}

main <- function(opt) {
  # Load the data set
  MSIG.Preprocess.Dataset(input.ds = opt$expfile, output.ds = paste(opt$output_prefix, ".normalized.gct", sep = ""), normalization = 4)
  norfile <- list.files(getwd(), pattern = "\\.normalized.gct", full.names = TRUE)
  # Run the main algorithm
  consensusNMF(input.ds = norfile[1], k.init = opt$k_int, k.final = opt$k_final, num.clusterings = 20, 
               maxniter = 1500, error.function = "divergence", rseed = 12345678, stopconv = 40, stopfreq = 10, 
               non.interactive.run = FALSE, doc.string = opt$output_prefix, libdir = libdir)
}

opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

# Prepare for parallel code
registerDoParallel(opt$cpu)

source_files(opt$libdir)

main(opt)
