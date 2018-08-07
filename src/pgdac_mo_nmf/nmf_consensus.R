

nmf2class <- function(cdesc, 
                      class.variable){
  
  ## #######################################################################################################
  ## add NMF to color list for pheatmap
  ## - map NMF classes to levels in 'class.variable'
  mb.nmf.map <- data.frame(mb=cdesc[, class.variable], nmf_basis=NMF.basis, nmf_cons=NMF.consensus)
  
  ## map consensus to 'class.variable'
  cons.map <- tapply(mb.nmf.map$mb, mb.nmf.map$nmf_cons, table)
  
  ## relative frequency
  cons.map.rel <- lapply(cons.map, function(x) x/sum(x))
  
  ## pick the maximum
  cons.map.max <- sapply(cons.map.rel, function(x) {m=max(x); ifelse( sum(x %in% m) > 1, paste(names(x)[which(x %in% m)] , collapse='|'), names(x)[which.max(x)])})
  
  
  ## assign colors 
  NMF.consensus.col <- cdesc.color[[class.variable]][ cons.map.max ]
  names(NMF.consensus.col) <- cons.map.max
  
  ######################################################  
  ## check whether all NMF basis have color
  if( sum(is.na(NMF.consensus.col)) > 0 ){
    
    #idx.tmp <- which( !(levels(NMF.consensus) %in% names(NMF.consensus.col)) )
    idx.tmp <- which( is.na(NMF.consensus.col) )
    
    col.tmp <-  rev( palette())[ 1:length(idx.tmp) ]
    #names(col.tmp) <- idx.tmp
    #NMF.consensus.col <- c(NMF.consensus.col, col.tmp)
    NMF.consensus.col[idx.tmp] <- col.tmp
  }
  cdesc.color$NMF.consensus <- NMF.consensus.col
  names(cdesc.color$NMF.consensus) <- 1:length(NMF.consensus.col)
  
}


## #############################################
## make a matrix non-negative
## separate up/down
make.nn <- function(m){
  

  ## up
  m.up <-  t( apply(m, 1, function(x){
    up=rep(0, length(x))
    up[which(x > 0)] <- x[which(x>0)]
    up
  }))
  colnames(m.up) <- colnames(m)
  rownames(m.up) <- paste(rownames(m), 'up', sep='_')
  ## down
  m.down <-  t( apply(m, 1, function(x){
    down=rep(0, length(x))
    down[which(x < 0)] <- abs(x[which(x<0)])
    down
  }))
  colnames(m.down) <- colnames(m)
  rownames(m.down) <- paste(rownames(m), 'down', sep='_')
  
  ## combine
  m.comb <- rbind(m.up, m.down)
  m.comb <- m.comb[order(rownames(m.comb)),]
  
  ## remove rows only containing zero values
  keep.idx <- which(apply(m.comb, 1, function(x) sum(x != 0) ) > 0)
  m.comb <- m.comb[keep.idx, ]
  
  return(m.comb)
}
  
  
  
## ###########################################################################################
## consensus clustering
consensus.clust <- function(m,                 ## data matrix
                          method=c('hclust', 'nmf', 'kmeans'), ## cluster method  
                          nmf.nrun=10,         ## number of random restarts for NMF
                          bs.nrun=5,           ## number of bootstrap runs for consensus
                          nmf.opts=NULL,       ## additional options for function NMF
                          nmf.method='brunet', ## nmf method
                          seed='random',       ## nmf seed method
                          rank=3,              ## nmf rank (number of clusters)
                          make.nn=FALSE,       ## make matrix non-negative?
                          ...
){
  
  if(!require(pacman)) install.packages('pacman')
  p_load(NMF)
  
  
  ## parse cluster method
  method <- match.arg(method)
  
  ## number of sample columns
  ns <- ncol(m)
  ## sample names
  samp <- colnames(m)
 
  ## vectors to count occurences of samp across bootstrap runs
  samp.count <- vector('numeric', length(samp))
  names(samp.count) <- samp
  
  # initialize consensus and indicator matrices
  M <- I <- matrix(0, nrow=ns, ncol=ns, dimnames = list(samp, samp))
  
  ## ################################
  ## bootstrap iterations
  for(bs in 1:bs.nrun){
    
    ## bootstrap
    samp.bs <- sample(samp, ns, replace=T)
    m.bs <- m[, samp.bs]
    keep.idx <- which(apply(m.bs, 1, function(x) sum(x != 0) ) > 0)
    m.bs <- m.bs[keep.idx, ]
    
    ## count bootstrap samples
    I[samp.bs, samp.bs] <- I[samp.bs, samp.bs] + 1
    
    ## nmf clustering
    if(method == 'nmf'){
        
      if(make.nn | min(m.bs, na.rm=T) < 0)
        m.bs <- make.nn(m.bs)
        
        nmf.bs <- nmf(m.bs, rank=rank, method=nmf.method, seed=seed, nrun=nmf.nrun, .options = nmf.opts)
        
        ## cluster assignment
        if(nmf.nrun > 1)
          bs.cons <- predict(nmf.bs, 'consensus')
        else
          bs.cons <- predict(nmf.bs)
    }
    if(method == 'hclust'){
      bs.cons <- cutree( hclust(dist(t(m.bs))), rank)
    }
    if(method == 'kmeans'){
      bs.cons <-  kmeans(t(m.bs), rank)$cluster
    }
    
    ## connectivity
    for(i in 1:rank){
      samp.i <- names(bs.cons[ bs.cons == i ])
      M[samp.i, samp.i] <-  M[samp.i, samp.i] + 1
    }

  } ## end boostrap iterations

  ## #################################
  ##  consensus matrix
  M <- M/I
  
  return(M)
  
}




# 
# ###########################################
# ## simulated data, 3 clusters
# data <- matrix( c(rnorm(10000, 100, 1),
#                   rnorm(10000, 101, 1),
#                   rnorm(10000, 99, 0.5)),
#                 nrow=1000, dimnames=list(paste('feature', 1:1000, sep='-'), paste('sample', 1:30,sep='-'))
#                 )
# 
# data2 <- syntheticNMF(1000, r = 3, p = 30)
# rownames(data2) <- paste('feature', 1:1000, sep='-')
# colnames(data2) <- paste('sample', 1:30,sep='-')
# 
# # 
#  data3 <- matrix( c(rnorm(10000, 0, 1),
#                    rnorm(10000, -.3, .1),
#                    rnorm(10000, .4, .5)),
#                  nrow=1000, dimnames=list(paste('feature', 1:1000, sep='-'), paste('sample', 1:30,sep='-'))
#  )
# 


#####################################################
#library(pheatmap)
#debug(nmf.consensus)

#opts <-  paste('vp3t', sep='')

#tmp<-consensus.clust(data, bnmf.opts=opts)
#tmp3<-consensus.clust(data2, method='nmf', bs.nrun = 100, nmf.nrun = 1)
#tmp3<-consensus.clust(data3, method='kmeans', bs.nrun = 100, nmf.nrun = 1)
#pheatmap(tmp3, symm=T)



