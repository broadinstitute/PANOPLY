
#
#  $Id$
#





by.collapse <- function (...) {
  # an extension of the 'by' function to return a data frame
  # instead of a list/array that by usually does

  result <- by (...)
  dimensions <- expand.grid ( dimnames (result) )
  # NB: In expand.grid, the first dimension varies fastest;
  #     the same happens in the output of by, and hence the
  #     table of dimensions and results match
  
  # if the result of the applied function is a single value ... use cbind ... else use do.call (rbind, ...)
  result.length <- length (result[[1]])
  if (result.length==1) data <- cbind ( unlist (lapply (result, unlist)) )
  else data <- do.call (rbind, list (result))
  
  final.result <- cbind (dimensions, data)
  invisible (final.result)
}





repeated <- function (data) {
  # function returns a boolean vector indicating whether each row in data is repeated;
  # this is differenent from the 'duplicated' function in that it includes the original
  # row along with the duplicates of that row

  # duplicated identifies the second and later replicates;
  # use duplicated after reversing the order of the rows to find the first occurance
  if (!is.data.frame (data)) data <- data.frame (cbind (data))
  reverse <- rev ( 1:nrow(data) )
  repeated <- duplicated (data) | duplicated (data[reverse,])[reverse]

  return (repeated)
}



ppm <- function (m) { return (m/1e6) }


ppm.var <- function (y) {
  # calculate ppm variation (range) for values in vector y
  r <-range (y)
  1e6*(r[2]-r[1])/r[2]
}

 

ppm.var.no.outliers <- function (y) {
  # calculate ppm variation (range) for values in vector y, after eliminating outliers
  y.quantiles <- boxplot (y, plot=FALSE)$stats
  r <- c (y.quantiles[1,1], y.quantiles[5,1])
  1e6*(r[2]-r[1])/r[2]
}



same.peak <- function (peak1, peak2, mz.threshold=25, rt.threshold=1e-3, intensity.threshold=1e-3) {
  # peak1 and peak2 should be lists with (mz, rt, charge, intensity)
  # determine if mz1 and mz2 are within mz.threshold ppm, provided charge is identical,
  # and the rt,intensity differenece is small (to account for rounding and other numeric errors)
  same <- ifelse (peak1$charge == peak2$charge && abs (peak1$rt-peak2$rt) <= rt.threshold &&
                  abs (peak1$intensity-peak2$intensity) <= intensity.threshold,
                  ppm.var ( c(peak1$mz, peak2$mz) ) <= mz.threshold, FALSE)
  return (same)
}


range.diff <- function (y) {
  # calculate variation (range) for values in vector y
  r <- range (y)
  r[2] - r[1]
}


range.diff.no.outliers <- function (y) {
  # calculate variation (range) for values in vector y, after eliminating outliers
  y.quantiles <- boxplot (y, plot=FALSE)$stats
  r <- c (y.quantiles[1,1], y.quantiles[5,1])
  r[2] - r[1]
}


get.subset <- function (data, start.n, ppm.tol) {
  # extracts and plots subset of 'data', to include all rows with m/z values within
  # ppm.tol of m/z for row 'start.n'
  mz.n <- data [start.n, 'mz']
  x <- data [ data[,'mz'] >= (mz.n - ppm.tol*ppm(mz.n)) & data[,'mz'] <= (mz.n + ppm.tol*ppm(mz.n)), ]
  plot (x[,c('rt','mz')], col=x[,'charge'], pch='+', cex=2, main=t)

  invisible (x)
}



calculate.nonzero.summary <- function (x, combine.fn=mean) {
  # calculate mean/median of vector x, after eliminating 0's
  x <- as.numeric(x)
  y <- x [ x!=0 ] 
  ifelse (length(y)==0, 0, combine.fn (y, na.rm=TRUE))
}



calculate.nonzero.cv <- function (x) {
  # calculate cv = sd/mean using only non-zero values in x
  mu <- calculate.nonzero.summary (x, mean)
  stdev <- calculate.nonzero.summary (x, sd)
  return ( stdev/mu )
}


consistent <- function (x) {
  # check if a vector representing peptide intensities in adjacent fractions is
  # "consistent" -- i.e., all non-zero intensities are adjacent to each other

  state <- 0        # start with 0 intensity
  n.changes <- 0    # number of state changes

  for ( i in 1:length(x) ) {
    if ( xor (state, x[i]>0) ) {
      # state and x[i] are different -- process state change
      state <- ! state
      n.changes <- n.changes + 1
    }
  }

  return (n.changes <= 2)
}




tolerant.merge <- function (data1, data2, tolerant.dim, by.dims,
                            tolerant.tolerance, by.tolerance, all=TRUE) {
  # merge data1 and data2, subject to the following conditions:
  # i.  the tolerant.dim (usually mz) should be within tolerant.tolerance ppm of each other
  # ii. all the other by.dims (usually rt, charge, intensity) should be within by.tolerance
  # if all==TRUE, all rows of data2 are included; if all=FALSE, only matched rows are retained

  all.dims <- c (tolerant.dim, by.dims)
  if ( ! (all.dims %in% colnames (data1)) ||
       ! (all.dims %in% colnames (data2)) )
    stop ('Required dimensions not present in data1 or data2.')
  add.dims <- setdiff ( colnames(data2), colnames(data1) )
  
  if (all) {
    data <- data1
    data.na <- data.frame (matrix (nrow=nrow(data), ncol=length(add.dims)))
    colnames (data.na) <- add.dims
    data <- cbind (data, data.na)
  } else {
    data <- NULL
  }

  temp <- lapply (1:nrow(data2),
                  function (i) {
                    row <- data2 [i, ]
                    # find closest item in tolerant.dim
                    diff <- abs (data1[,tolerant.dim] - row[1,tolerant.dim])
                    diff.index <- diff == min (diff)
                    # find items satisfying specified by.dims tolerances
                    for (d in by.dims) {
                      diff.d <- abs (data1[,d] - row[1,d])
                      diff.index <- diff.index & ( diff.d <= by.tolerance )
                    }
                    # closest item -- should be unique
                    closest <- which (diff.index)

                    # ensure closest item (if it exists) satisfies tolerance for tolerant.dim
                    if ( length (closest) > 0 ) {
                      same <- ppm.var ( c(row[1,tolerant.dim], data1[closest,tolerant.dim]) ) <= tolerant.tolerance
                      if (same) {
                        if (all) {
                          for (d in add.dims) 
                            data [closest, d] <<- ifelse (is.factor (data2[,d]), toString (row[1,d]), row[1,d])
                        } else {
                          new.row <- data1 [closest,]
                          for (d in add.dims) 
                            new.row <- c (new.row, ifelse (is.factor (data2[,d]), toString (row[1,d]), row[1,d]))
                          data <<- rbind (data, new.row)
                        }
                      }
                    }
                  })
  return (data)
}

  
    

rms <- function (x, y) {
  # calculates RMS (root-mean-square) value of x and y, after eliminating missing value
  index <- !is.na (x) & !is.na (y)
  x <- x [index]
  y <- y [index]
  sqrt ( sum ((x-y)^2) / length(x) )
}



mz.merge <- function (data1, data2, mz.tolerance, rt.tolerance, suffixes=c('.1','.2'),
                      outfile=NULL, parallel=FALSE, debug=0) {
  # merge two datasets by mz, rt and charge
  # mz and rt should be within specified tolerance
  # charge should be identical
  # only common rows are retained
  # columns having the same name in data1 and data2 are disambiguated using suffixes
  
  # make sure required dimensions are present
  dims <- c ('mz', 'rt', 'charge')
  if ( ! all ( dims %in% colnames (data1) ) ||
       ! all ( dims %in% colnames (data2) ) )
    stop ('Required dimensions not present in data1 or data2.')

  # disambiguate and order dimensions
  dims1 <- setdiff (colnames (data1), dims)
  dims2 <- setdiff (colnames (data2), dims)
  # for data1
  dims1.common <- dims1 %in% dims2
  dims1.orig <- c (dims, dims1 [! dims1.common], dims1[dims1.common])
  dims1.new <- c (dims, dims1 [! dims1.common],
                  paste (dims1[dims1.common], suffixes[1], sep=''))
  # for data2
  dims2.common <- dims2 %in% dims1
  dims2.orig <- c (dims, dims2 [! dims2.common], dims2[dims2.common])
  dims2.new <- c (paste (dims, suffixes[2], sep=''), dims2 [! dims2.common],
                  paste (dims2[dims2.common], suffixes[2], sep=''))

  d1 <- data1
  d2 <- data2
  if (! parallel) {
    # sort the data by mz
    # (will already be sorted for a parallel call)
    d1 <- data1 [ sort (data1[,'mz'], index=T)$ix, ]
    d2 <- data2 [ sort (data2[,'mz'], index=T)$ix, ]
  }

  # convert factors to strings so that names, etc. (if present) are preserved
  for (l in dims1) if (is.factor (d1[,l])) d1[,l] <- unlist (lapply (d1[,l], toString))
  for (l in dims2) if (is.factor (d2[,l])) d2[,l] <- unlist (lapply (d2[,l], toString))


  i.1 <- i.2 <- 1
  data <- NULL
  while ( i.1 <= nrow(d1) && i.2 <= nrow(d2) ) {
    if (debug) print ( paste (' ... matching rows', i.1, 'and', i.2), quote=FALSE )
    mz.diff <- d1[i.1,'mz'] - d2[i.2,'mz']
    ppm.diff <- ppm.var ( c (d1[i.1,'mz'], d2[i.2,'mz']) )
    if (ppm.diff > mz.tolerance) {
      # current rows don't match
      if (mz.diff < 0) i.1 <- i.1 + 1
      else i.2 <- i.2 + 1
    } else {
      # current mz's match
      # find set of mz's in d2 within tolerance
      j.end <- j.start <- i.2
      ppm.diff.list <- ppm.diff
      while ( i.1 <= nrow (d1) && (j.end+1) <= nrow (d2) &&
              (mz.ppm <- ppm.var ( c (d1[i.1,'mz'], d2[j.end+1,'mz']) )) <= mz.tolerance ) {
        j.end <- j.end + 1
        ppm.diff.list <- c (ppm.diff.list, mz.ppm)
      }
      # within rows that match mz, check for rt and charge match
      match.all <- abs ( d2 [j.start:j.end,'rt'] - d1[i.1,'rt'] ) <= rt.tolerance &
                   d2 [j.start:j.end,'charge'] == d1[i.1,'charge']
      match.final <- (j.start:j.end)[match.all]

      if (length (match.final) > 1) {
        # too many matches -- pick the closed mz
        ppm.matches <- ppm.diff.list [ match.all ]
        match.final <- (j.start:j.end) [which ( ppm.diff.list == min (ppm.matches) )]
      }

      if (length (match.final) == 1) {
        # found the best match
        matched.row <- c ( d1[i.1, dims1.orig], d2[match.final, dims2.orig] )
        data <- rbind (data, matched.row)
      }

      i.1 <- i.1 + 1   # i.2 was part of a block, so don't increment it
    }
  }
  if (! is.null (dim (data)) )
    # data has some rows
    colnames (data) <- c (dims1.new, dims2.new)
  
  if (! is.null(outfile)) write.table (data, outfile, sep=',', row.names=F, quote=F)
  invisible (data)
}



parallel.mz.merge <- function (file1, file2, outfile, mz.tolerance, rt.tolerance, suffixes=c('.1','.2'),
                               nproc=100, scratch.dir='LSF', parallelization.limit=1000, debug=0) {
  # splits data in to multiple strips and performs tolerant merge in parallel
  # merge two datasets by mz, rt and charge
  # mz and rt should be within specified tolerance
  # charge should be identical
  # only common rows are retained
  # columns having the same name in data1 and data2 are disambiguated using suffixes



  split.for.mz.merge <- function (data1, data2, mz.tolerance, output.prefix,
                                  data1.prefix="d1", data2.prefix="d2", nproc=100) {
    # split datasets data1 and data2 into nproc pieces so that the merge
    # can be run in parallel


    # make sure required dimensions are present
    dims <- c ('mz', 'rt', 'charge')
    if ( ! all ( dims %in% colnames (data1) ) ||
         ! all ( dims %in% colnames (data2) ) )
      stop ('Required dimensions not present in data1 or data2.')

    # sort the data by mz
    # this is required for correct partitioning
    d1 <- data1 [ sort (data1[,'mz'], index=T)$ix, ]
    d2 <- data2 [ sort (data2[,'mz'], index=T)$ix, ]

    # size of d1 chunks
    k <- ceiling ( nrow (d1)/ nproc )
    
    for (i in 1:nproc) {
      # find d1 chunk
      start.1 <- k * (i-1) + 1
      end.1 <- k * i
      if (end.1 > nrow (d1)) end.1 <- nrow (d1)
      
      # find d2 chunk to correspond with d1 chunk
      # note that d2 chunks can overlap to account for multiple rows with mz's close to
      # start.1/end.1 mz's
      start.2.ppm <- abs (1e6 * ( d1[start.1, 'mz'] - d2[,'mz'] ) / d1[start.1, 'mz'])
      start.2.list <- which (start.2.ppm <= mz.tolerance)
      start.2 <- ifelse (length(start.2.list) > 0,
                         min (start.2.list),
                         which (start.2.ppm == min (start.2.ppm)))   # if no mz was within tolerance,
                                                                     # pick the closest mz
      
      end.2.ppm <- abs (1e6 * ( d1[end.1, 'mz'] - d2[,'mz'] ) / d1[end.1, 'mz'])
      end.2.list <- which (end.2.ppm <= mz.tolerance)
      end.2 <- ifelse (length(end.2.list) > 0,
                       max (end.2.list),
                       which (end.2.ppm == min (end.2.ppm)))         # if no mz was within tolerance,
                                                                     # pick the closest mz

      if (debug) {
        print ( paste ('    split:', i), quote=FALSE )
        print ( paste ('      start.1:end.1 = ', start.1, ':', end.1, sep=''), quote=FALSE )
        print ( paste ('      start.2:end.2 = ', start.2, ':', end.2, sep=''), quote=FALSE )
      }
      
      # write output files
      write.table ( d1[ start.1:end.1, ], paste (output.prefix, data1.prefix, i, 'csv', sep='.'),
                    sep=',', row.names=FALSE, quote=FALSE )
      write.table ( d2[ start.2:end.2, ], paste (output.prefix, data2.prefix, i, 'csv', sep='.'),
                    sep=',', row.names=FALSE, quote=FALSE )
    }
  }
    

  print ('Reading data files ...', quote=FALSE)
  data1 <- read.csv (file1, comment.char='')
  data2 <- read.csv (file2, comment.char='')

  
  # make sure required dimensions are present
  dims <- c ('mz', 'rt', 'charge')
  if ( ! all ( dims %in% colnames (data1) ) ||
       ! all ( dims %in% colnames (data2) ) )
    stop ('Required dimensions not present in data1 or data2.')


  if (nrow (data1) < parallelization.limit) {
    # no need to parallelize small datasets
    mz.merge (data1, data2, mz.tolerance, rt.tolerance, suffixes=suffixes, outfile=outfile,
              parallel=FALSE, debug=debug)
  } else {
    # parallelize merge
    print ('Parallelizing merge ...', quote=FALSE)
    data1.prefix <- 'd1'
    data2.prefix <- 'd2'
    
    pid <- Sys.getpid()
    temp.dir <- paste (scratch.dir, '.', pid, sep='')
    dir.create (temp.dir)
    setwd (temp.dir)

    tmp.data.prefix <- paste ('temp', pid, sep='.')
    r.file.prefix <- paste ('temp-r', pid, 'r', sep='.')
    outfile.prefix <- paste (tmp.data.prefix, 'merged', sep='.')

    print ('  splitting data for parallel processes', quote=FALSE)
    split.for.mz.merge (data1, data2, mz.tolerance, tmp.data.prefix,
                        data1.prefix=data1.prefix, data2.prefix=data2.prefix, nproc=nproc)

    print ('  starting parallel jobs', quote=FALSE)
    for (i in 1:nproc) {
      if (debug) print ( paste ('    process:', i), quote=FALSE )
      r.file <- paste (r.file.prefix, i, sep='.')
      script <- file (r.file, 'w')
      # assmemble R call, and dispatch to LSF
      cat ("code.dir <- Sys.getenv ('PROTEOMICS_CODEBASE')\n", file=script)
      cat ("source  ( paste (code.dir, 'lc-ms.r', sep='/') )\n\n", file=script)
      cat ("data1 <- read.csv ('", paste (tmp.data.prefix, data1.prefix, i, 'csv', sep='.'),
           "', comment.char='')\n", file=script, sep='')
      cat ("data2 <- read.csv ('", paste (tmp.data.prefix, data2.prefix, i, 'csv', sep='.'),
           "', comment.char='')\n", file=script, sep='')
      merge.cmd <- paste ('mz.merge (data1, data2, ', mz.tolerance, ', ', rt.tolerance,
                          ',  suffixes=c("', suffixes[1], '","', suffixes[2], '")',
                          ', outfile="', paste (outfile.prefix, i, 'csv', sep='.'),
                          '", parallel=TRUE, debug=', debug, ')', sep='')
      cat (merge.cmd, file=script)
      close (script)

      cmd <- paste ('bsub -J mz_merge_', pid, ' R CMD BATCH --vanilla ', r.file, sep='')
      system (cmd)
    }

    # wait for parallel jobs to finish and assemble results
    print ('Waiting for parallel jobs to finish', quote=FALSE)
    while ( as.numeric (system (paste ('bjobs -J mz_merge_', pid, ' | wc -l', sep=''),
                                intern=TRUE)) > 1 ) Sys.sleep (120)

    print ('Assembing final merged table', quote=FALSE)
    data.merge <- NULL
    for (i in 1:nproc) {
      result <- read.csv ( paste (outfile.prefix, i, 'csv', sep='.'), comment.char='' )
      data.merge <- rbind (data.merge, result)
    }
    
    setwd ('..')

    print ('Writing output table', quote=FALSE)
    write.table (data.merge, outfile, sep=',', row.names=FALSE, quote=FALSE)
  }
}
