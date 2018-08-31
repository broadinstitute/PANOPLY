

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
## extract tarbar and read configuration file
prepare.data.sets <- function( tar.file, tmp.dir){
  
  ## #################################
  ## extract tar ball
  if(!dir.exists(tmp.dir))
    dir.create(tmp.dir)
  untar(tar.file, exdir=tmp.dir)
  
  ##  import config file
  conf <- read.delim(paste( tmp.dir, 'nmf.conf', sep='/'), row.names = NULL, stringsAsFactors = F, header=F)
  data.str <- paste('..', tmp.dir, conf[, 2], sep='/')
  names(data.str) <- conf[,1]
  
  return(data.str)
}

## #############################################
## import and merge data tables
#import.data.sets <- function(data.str, zscore.cnv=F){
import.data.sets <- function(tar.file, tmp.dir, zscore.cnv=F){
  
  ## ################################
  ## extract tar file
  data.str <- prepare.data.sets(tar.file, tmp.dir)
  
  ## ################################
  ## loop over file path
  for(i in 1:length(data.str)){
    
    ##gct
    gct <-  parse.gctx2(data.str[i])
    gct@cid <- make.names(gct@cid)
    ## extract data matrix
    data.tmp <- gct@mat
    ## append data type to row names
    rownames(data.tmp) <- paste(names(data.str)[i], gct@rid, sep='-')
    colnames(data.tmp) <- gct@cid
    
    ## remove columns with no data values
    keep.idx <- which(apply(data.tmp, 2, function(x) ifelse(sum(is.na(x))/length(x) == 1, F, T)))
    data.tmp <- data.tmp[, keep.idx]
    
    ## row annotations    
    rdesc.tmp <- gct@rdesc
    if(nrow(rdesc.tmp) == 0){
      rdesc.tmp <- data.frame(id2=gct@rid)
    }
    rownames(rdesc.tmp) <- paste(names(data.str)[i], gct@rid, sep='-')
    rdesc.tmp <- data.frame( Data.Type.ID=rownames(rdesc.tmp), rdesc.tmp)
    
    ## column annotations
    cdesc.tmp <- gct@cdesc %>% data.frame
    #if(nrow(cdesc.tmp) == 0){
    #  cdesc.tmp <- data.frame(id2=gct@cid)
    #}
    
    #cdesc.tmp[] <- lapply(cdesc.tmp, as.character) %>% as.data.frame()
    if(nrow(cdesc.tmp) > 0){
      cdesc.tmp <- data.frame(ID=gct@cid, cdesc.tmp)
      cdesc.tmp <- cdesc.tmp[ keep.idx,  ]
      rownames(cdesc.tmp) <- gct@cid
    }
    
    ## Z-score CNV?
    if(names(data.str)[i] == 'CNV' & zscore.cnv)
      data.tmp <- apply(data.tmp, 2, function(x) (x-median(x, na.rm=T))/mad(x, na.rm=T))
    
    ## merge data tables, row and column annotations
    if(i == 1){
      data <- data.tmp
      samp <- colnames(data.tmp)
      rdesc <- rdesc.tmp
      cdesc <- cdesc.tmp
      
    } else {
      ## only sample columns present in all tables will be carried over
      samp <- intersect(samp, colnames(data.tmp))
      data <- data[, samp]
      cdesc <- cdesc[samp, ]
      cdesc[] <- lapply(cdesc, function(x){if(sum(is.na(as.numeric(x))) == 0){as.numeric(x);} else {x} })
      
      ## append to data matrix
      data <- rbind(data, data.tmp[, colnames(data)])
      
      ## row annotations
      rdesc <- full_join(rdesc, rdesc.tmp)
      ## add back rownames
      rownames(rdesc) <- rdesc$Data.Type.ID
      ## ensure same order as data
      rdesc <- rdesc[rownames(data), ]
      
      ## column annotations
      if(nrow(cdesc.tmp) > 0){
        cdesc.tmp[] <- lapply(cdesc.tmp, function(x){if(sum(is.na(as.numeric(x))) == 0){as.numeric(x);} else {x} })
        #cdesc.tmp[] <- lapply(cdesc.tmp, as.character) %>% as.data.frame()
        cdesc <- left_join(cdesc, cdesc.tmp)
        rownames(cdesc) <- cdesc$ID
        #cdesc2 <- full_join(data.frame(t(cdesc)), data.frame(t(cdesc.tmp)))
      }
    }
  }
  
  ## expression data
  expr <- data.matrix(data)
  
  return(list(expr=expr, cdesc=cdesc, rdesc=rdesc))

}


## #############################################
## make a matrix non-negative
## separate up/down
make.non.negative <- function(m){
  
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
  
nmf.tmp <- function(data, k){
  
  #data <- make.non.negative(data)
  clust <- cluster(data, rank=k,  method='nmf', nmf.method='brunet', seed='random', nmf.nrun=3, make.nn = T, nmf.opts=NULL)
  return(clust)
} 
nmf.dist <- function(data) return(data)

## ###########################################################################################
## consensus clustering
consensus.clust <- function(m,                 ## data matrix p x n
                          method=c('hclust', 'nmf', 'kmeans'), ## cluster method  
                          nmf.nrun=10,         ## number of random restarts for NMF
                          bs.nrun=5,           ## number of bootstrap runs for consensus
                          nmf.opts=NULL,       ## additional options for function NMF
                          nmf.method='brunet', ## nmf method
                          seed='random',       ## nmf seed method
                          rank=3,              ## nmf rank (number of clusters)
                          make.nn=FALSE,       ## make matrix non-negative?
                          ncore=2,             ## number of cores
                          plot=F,              ## if TRUE, a heatmap of the consensus matrix will be plotted
                          ...
){
  
  if(!require(pacman)) install.packages('pacman')
  p_load(NMF)
  p_load(doParallel)
  p_load(pheatmap)
  p_load( RColorBrewer)
  #p_load(doFuture)
  
  ## parse cluster method
  method <- match.arg(method)
  
  
  ## ############################
  ## helper function
  cluster <- function(data,  ## p x n matrix, p-features, n-samples
                      method='nmf'
                      , rank, seed, nmf.method, nmf.nrun, nmf.opts, make.nn){
    #cat('\njep\n')
    ## #############################################
    ## helper function
    ## make a matrix non-negative
    
    ## separate up/down
    make.non.negative <- function(m){
      #cat('\ntest..\n')
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
    #debug(make.non.negative)
    ## nmf clustering
    if(method == 'nmf'){
      
      if(make.nn | min(data, na.rm=T) < 0)
        data <- make.non.negative(data)
      
      nmf.bs <- nmf(data, rank=rank, method=nmf.method, seed=seed, nrun=nmf.nrun, .options = nmf.opts)
      
      ## cluster assignment
      if(nmf.nrun > 1)
        bs.cons <- predict(nmf.bs, 'consensus')
      else
        bs.cons <- predict(nmf.bs)
    }
    if(method == 'hclust'){
      bs.cons <- cutree( hclust(dist(t(data))), rank)
    }
    if(method == 'kmeans'){
      bs.cons <-  kmeans(t(data), rank)$cluster
    }
    
    return(bs.cons)
  } ## end function cluster
  
  ## ###############################
  
  
  ## number of sample columns
  ns <- ncol(m)
  
  ## sample names
  samp <- colnames(m)
 
  ## vectors to count occurences of samp across bootstrap runs
  samp.count <- vector('numeric', length(samp))
  names(samp.count) <- samp
  
  
  ## ################################
  ## parallelize
  ncore.tot <- detectCores()
  ncore <- min(ncore.tot, ncore)
  
  cat('\nDetected', ncore.tot+1, 'cores.\nUsing', ncore, 'cores for bootstrapping.\n')
  
  if(method=='nmf'){
      ncore.nmf <- max( ncore.tot - ncore, 1)
      if(is.null(nmf.opts)) nmf.opts <- paste('vp', ncore.nmf, 't', sep='')
      cat('Using', ncore.nmf, 'cores for NMF.\n')
      
      }
  
  
  
  ## ################################
  ## bootstrap iterations
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  #cl <- makeCluster(ncore)
  
  
  #for(bs in 1:bs.nrun){
  M2 <- foreach(bs = 1:bs.nrun) %dopar% {
  #M2 <- list()
  #for(bs in 1:bs.nrun){
    library(NMF)
    
    # initialize consensus and indicator matrices
    M <- I <- matrix(0, nrow=ns, ncol=ns, dimnames = list(samp, samp))
    
    ## bootstrap
    samp.bs <- sample(samp, ns, replace=T)
    m.bs <- m[, samp.bs]
    keep.idx <- which(apply(m.bs, 1, function(x) sum(x != 0) ) > 0)
    m.bs <- m.bs[keep.idx, ]
   
    ## indicator matrix
    #I[samp.bs, samp.bs] <- I[samp.bs, samp.bs] + 1
    I[samp.bs, samp.bs] <- 1
    
    ## cluster
    bs.clust <- cluster(data=m.bs, method=method, 
                       rank=rank, seed=seed, nmf.method=nmf.method, 
                       nmf.nrun=nmf.nrun, nmf.opts=nmf.opts, make.nn=make.nn)

    ## connectivity matrix
    for(i in 1:rank){
      samp.i <- names(bs.clust[ bs.clust == i ])
      #M[samp.i, samp.i] <-  M[samp.i, samp.i] + 1
      M[samp.i, samp.i] <- 1
      
    }
    
    list(M=M, I=I)
  } ## end boostrap iterations
  
  on.exit(stopCluster(cl))
  
  ## #################################
  ##  consensus matrix
  M.list <- lapply(M2, function(x)x[['M']])
  I.list <- lapply(M2, function(x)x[['I']])
  M <- M.list[[1]]
  I <- I.list[[1]]
  
  for(i in 2:length(M.list)){
    M <- M + M.list[[i]]
    I <- I + I.list[[i]]
  }
  M.abs <- M
  M <- M/I
  
  ## #################################
  ## clustering of actual data
  cat('\nPerforming clustering on unperturbed data\n')
  clust <- cluster(data=m, method=method, 
                   rank=rank, seed=seed, nmf.method=nmf.method,
                   nmf.nrun=nmf.nrun, nmf.opts=nmf.opts, make.nn=make.nn)
  
  ## ##################################
  ## plot
  if(plot){
    anno <- data.frame(cluster=as.factor(clust))
    clust.col <- RColorBrewer::brewer.pal(rank, "Dark2")
    names(clust.col) <- unique(anno$cluster)
    anno.col <- list(cluster=clust.col)
    
    hm.col <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)
    breaks <- seq(0,1, length.out = length(hm.col)+1)
    
    pheatmap(M, symm=T, 
             col=hm.col,
             breaks=breaks,
             annotation_col = anno,
             annotation_row = anno,
             annotation_colors = anno.col)
  }
  
  ## assemble output
  out <-c()
  out$M <- M
  out$M.abs <- M.abs
  out$I <- I
  out$clust <- clust
  return(out)
}


## #################################################################
## consensus clustering:
## calculate cluster consensus m_k
##
clust_cons <- function(M,      ## consensus matrix
                       clust   ## cluster assignment for rows/columns in M
){
  ## number of clusters
  K <- max(clust)
  
  m_k <- sapply(1:K, function(x) {
    m=M[ clust==x, clust==x ]
    m=m[upper.tri(m, diag = F)]
    mean(m)
    } )
  return(m_k)
}

## #################################################################
## consensus clustering:
## calculate item's consensus
item_cons <- function(M,
                      clust){
  
}




###########################################
## simulated data, 3 clusters
#data.sim <- matrix( c(rnorm(1000, 101, 1),
#                  rnorm(1000, 103, 1),
#                  rnorm(1000, 99, 0.5)),
#               nrow=100, dimnames=list(paste('feature', 1:100, sep='-'), paste('sample', 1:30,sep='-'))
#                )
data.sim <- matrix( c(rnorm(1000, 0, .6),
                  rnorm(1000, -0.5, .7),
                  rnorm(1000, 0.6, 0.5)),
               nrow=100, dimnames=list(paste('feature', 1:100, sep='-'), paste('sample', 1:30,sep='-'))
                )
pca <- prcomp(t(data.sim))

plot(pca$x[,1], pca$x[,2], col=p$clust, pch=20, cex=2)

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



