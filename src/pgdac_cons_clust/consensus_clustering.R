
if(!require('pacman')) install.packages('pacman')
p_load(NMF)
p_load(NbClust)
p_load(factoextra)
p_load(glue)
p_load(magrittr)
p_load(cmapR)
p_load(doParallel)


## ###########################################################################################
## consensus clustering
consensus_clustering_single_k <- function(m,                 ## data matrix p x n
                          method=c('hclust', 'nmf', 'kmeans'), ## cluster method  
                          nmf.nrun=10,         ## number of random restarts for NMF
                          bs.nrun=5,           ## number of bootstrap runs for consensus
                          nmf.opts=NULL,       ## additional options for function NMF
                          nmf.method='brunet', ## nmf method
                          seed='random',       ## nmf seed method
                          k=3,                 ## number of clusters
                          make.nn=FALSE,       ## make matrix non-negative?
                          ncore=2,             ## number of cores
                          plot=F,              ## if TRUE, a heatmap of the consensus matrix will be plotted
                          logfile="",
                          ...
){
  
  if(!require(pacman)) install.packages('pacman')
  p_load(NMF)
  p_load(doParallel)
  p_load(pheatmap)
  p_load( RColorBrewer)
  p_load(factoextra)
  #p_load(doFuture)
  
  ## parse cluster method
  method <- match.arg(method)
  
  ## ############################
  ## helper function
  cluster <- function(data,                    ## p x n data matrix, p-features, n-samples
                      d=NULL,                  ## distance matrix, if specified hierarchical clustering will be performed
                      method='nmf',            ## character: hclust, kmeans, nmf
                      k,                       ## numeric, number of clusters 
                      seed='random',           ## nmf: seed to initialize matrices
                      nmf.method='brunet',     ## nmf: factorization method, default is 'brunet' 
                      nmf.nrun,                ## nmf: number of random restarts
                      nmf.opts=NULL,           ## nmf: further options passed to function 'nmf' 
                      make.nn=T                ## nmf: logical, if TRUE 'data' will be made non-negative 
                      ){
    
    ## #############################################
    ## helper function
    ## make a matrix non-negative
    ## - separate up/down
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
    
    ## if a distance matrix was specified
    if(!is.null(d)){
        method <- 'hclust'
    }
    
    ## ###########################
    ## nmf clustering
    if(method == 'nmf'){
      
      if(make.nn | min(data, na.rm=T) < 0)
        data <- make.non.negative(data)
      
      clust <- nmf(data, rank=k, method=nmf.method, seed=seed, nrun=nmf.nrun, .options = nmf.opts)
      
      ## cluster assignment
      if(nmf.nrun > 1)
        membership <- predict(clust, 'consensus')
      else
        membership <- predict(clust)
      
      clust$cluster <- membership
    }
    ## ##########################
    ## hierarchical clustering
    if(method == 'hclust'){
      
      if(is.null(d))
        d <- dist(t(data), method = 'euclidean')
      
      clust <- hclust(d, method = 'ward.D2')#'average')
      membership <- cutree( clust, k)
      clust$cluster <- membership
    }
    
    ## ###########################
    ## k-means clustering
    if(method == 'kmeans'){
      clust <- kmeans(t(data), k)
      membership <-  clust$cluster
    }
    
    ## assemble output
    out <- c()
    out$membership <- membership
    out$cluster.object <- clust
    out$cluster.object$dist <- d
    
    return(out)
  } ## end function cluster
  ## ###############################
  cat(glue('\nk={k} \n'))
  cat(glue('\n\n{paste(rep("#", 20), collapse="")}\nPartitioning data in k={k} clusters using {method}-clustering.\n'), file=logfile, append=T)
  
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
  ncore <- min(ncore.tot-1, ncore)
  
  cat('\nDetected', ncore.tot, 'cores.\nUsing', ncore, 'cores for bootstrapping.\n', file=logfile, append=T)
  
  if(method=='nmf'){
      ncore.nmf <- max( ncore.tot - ncore, 1)
      if(is.null(nmf.opts)) nmf.opts <- paste('vp', ncore.nmf, 't', sep='')
      cat('Using', ncore.nmf, 'cores for NMF.\n')
  }
  
  ## ################################
  ## bootstrap iterations
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  
  cat('Bootstrapping (', bs.nrun, 'iterations).\n', file=logfile, append=T)
  M2 <- foreach(bs = 1:bs.nrun) %dopar% {
    library(NMF)
    
    # initialize consensus and indicator matrices
    M <- I <- matrix(0, nrow=ns, ncol=ns, dimnames = list(samp, samp))
    
    ## bootstrap
    samp.bs <- sample(samp, ns, replace=T)
    m.bs <- m[, samp.bs]
    keep.idx <- which(apply(m.bs, 1, function(x) sum(x != 0) ) > 0)
    m.bs <- m.bs[keep.idx, ]
   
    ## indicator matrix
    I[samp.bs, samp.bs] <- 1
    
    ## cluster
    clust.bs <- cluster(data=m.bs, method=method, 
                       k=k, seed=seed, nmf.method=nmf.method, 
                       nmf.nrun=nmf.nrun, nmf.opts=nmf.opts, make.nn=make.nn)
    membership.bs <- clust.bs$membership
    
    ## connectivity matrix
    for(i in 1:k){
      samp.i <- names(membership.bs[ membership.bs == i ])
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
  ## final clustering on consensus matrix
  cat('\nPerforming clustering on consenus matrix\n', file=logfile, append=T)
  #clust <- cluster(data=m, method=method, 
  clust <- cluster(data=NULL, 
                   d=as.dist(1 - M),
                   method='hclust', 
                   k=k, seed=seed, nmf.method=nmf.method,
                   nmf.nrun=nmf.nrun, nmf.opts=nmf.opts, make.nn=make.nn)
  cat('done\n', file=logfile, append=T)
  ## ##################################
  ## plot
  if(plot){
    
    hm.consensus(M, clust$membership )

    ## cluster plot
    fviz_cluster_plot( m,  clust$membership )
  }
  
  ## assemble output
  out <-c()
  out$M <- M
  out$M.abs <- M.abs
  out$I <- I
  out$cluster.results <- clust
  
  return(out)
}

## ######################################################################
## performs the actual consensus clustering
## - loop over cluster numbers k.min to k.max
## - perform bootstraping-based consensus clustering
## - determine optimal number of clusters
##
consensus_clustering <- function(m,                 ## data matrix p x n
                                 method=c('hclust', 'nmf', 'kmeans'), ## cluster method  
                                 nmf.nrun=10,         ## number of random restarts for NMF
                                 bs.nrun=5,           ## number of bootstrap runs for consensus
                                 nmf.opts=NULL,       ## additional options for function NMF
                                 nmf.method='brunet', ## nmf method
                                 seed='random',       ## nmf seed method
                                 k.min=2,             ## min. number of clusters
                                 k.max=4,             ## max. number of clusters
                                 make.nn=FALSE,       ## make matrix non-negative? NMF only
                                 ncore=2,             ## number of cores
                                 plot=T,              ## if TRUE, a heatmap of the consensus matrix will be plotted
                                 prefix='out',
                                 cdesc=NULL,
                                 ...
){
  
  
  method <- match.arg(method)
  
  
  ## ################################
  ## loop over cluster numbers  
  cons.res <- lapply( k.min:k.max, function(k) 
    consensus_clustering_single_k(m, method=method, nmf.nrun=nmf.nrun, bs.nrun=bs.nrun, nmf.opts=nmf.opts, seed=seed, k=k,
                                  make.nn=make.nn, ncore=ncore, plot=F, ...) 
  )
  names(cons.res) <- k.min:k.max
  
  
  ## #################################
  ## select best number of clusters  
  k_opt <- select_best_k(cons.res, data=m, plot=plot, prefix=prefix)
  K <- k_opt$k.opt[1]
  
  ## #################################
  ## extract results for optimal K
  res.opt <- cons.res[[as.character(K)]]
  write.table(k_opt$all.metrics, file=glue('{prefix}_cluster_metrics.csv'), row.names = T, quote = F, col.names = NA, sep=',')
  write.table(k_opt$best.metrics, file=glue('{prefix}_best_K.csv'), row.names = T, quote=F, col.names = NA, sep=',')
  
  

  ## #################################
  ## export final clustering
  cluster.final <- res.opt$cluster.results$membership
  if(plot)
    pdf(glue('{prefix}.bestclus.silhouette.pdf'))
  silhoutte.final <- silhouette_scores(res.opt$cluster.results$cluster.object)
  sil.val <- silhoutte.final[, 'sil_width']
  if(plot)
    dev.off()
  
  ## assemble final table
  bestclus <- data.frame(
    SampleName=colnames(m),
    cluster=cluster.final,
    SilhouetteValue=sil.val
  )
  fn <- glue ("{prefix}.bestclus.txt")
  write.table(bestclus, file=fn, row.names=F, sep='\t', quote=F)
  
  ## ##################################
  ## plot
  if(plot){
    
    ## ##############################
    ## best K
    hm.consensus(M = res.opt$M, membership = res.opt$cluster.results$membership,
                 cdesc=cdesc,
                 show_rownames=F,
                 show_colnames=F,
                 filename=glue('{prefix}_consensus_matrix_K{K}.pdf'))  
    
    pdf(glue('{prefix}_cluster_pca_K{K}.pdf'))
    fviz_cluster_plot(m, res.opt$cluster.results$membership)
    dev.off()
    
    ## ###############################
    ## all K
    dummy <- sapply(names(cons.res), function(k){
      cl=cons.res[[k]]
      hm.consensus(cl$M, cl$cluster.results$membership, cdesc=cdesc, 
                   show_rownames=F, show_colnames=F, 
                   filename=glue('{prefix}_consensus_matrix_k{k}.png'))
      })
  }
  
  
  ## assemble output
  
  return(cons.res)
  
}


## ###############################################################
## select the best number of clusters
select_best_k <- function(cons.res, 
                          data,              ## the data
                          delta.auc.sig=0.1, ## threshold for delta auc of consensus cdf  
                          plot=T,
                          prefix='out'){
  
  k.max <- names(cons.res) %>% as.numeric %>% max
  k.min <- names(cons.res) %>% as.numeric %>% min
  n.clust <- length(cons.res)
  k.all <- names(cons.res)
  
  ## ##########################################
  ## cluster metrics
  
  ## delta auc
  clust.auc <- sapply(cons.res, function(k) auc_cdf_cons(k$M) )
  clust.auc.delta <-  sapply(1:(n.clust), function(k, auc) {
    if(k == 1) return(auc[k])
    return( (auc[k] - auc[k-1])/auc[k-1])
  }, clust.auc)
  clust.auc.delta.diff <- sapply(1:(n.clust-1), function(i, val){ val[i+1] - val[i]}, clust.auc.delta ) 
  clust.auc.delta.diff <- c(0, clust.auc.delta.diff) %>% abs
  
  ## calculate silhouette scores
  clust.sil <- lapply( cons.res, function(k) silhouette_scores(k$cluster.results$cluster.object, data=data, plot=F))
  clust.sil.avg <- sapply(clust.sil, function(k) mean( k[, 'sil_width'], na.rm=T ))
  
  ## cophenetic correlation
  clust.coph <- sapply( cons.res, function(k) cophenetic_correlation(k$M))
  clust.coph.diff <- sapply(1:(n.clust-1), function(i, val){ val[i+1] - val[i]}, clust.coph ) 
  clust.coph.diff <- c(0, clust.coph.diff)
  

  ## ###############################################################
  ## assemble all metrics into a data frame
  cm.all <- data.frame(#cdf.auc=clust.auc,
    delta.auc=clust.auc.delta,
    delta.auc.diff=clust.auc.delta.diff,
    silhouette.score=clust.sil.avg,
    cophenetic.correlation=clust.coph,
    cophenetic.correlation.diff=clust.coph.diff,
    stringsAsFactors = F
  )
  
  ## ###########################################
  ## determine best k for each metric
  cm.best <- matrix(NA, nrow=ncol(cm.all), ncol=3, dimnames=list(colnames(cm.all), c('best.k', 'best.k.idx', 'best.k.score')))
  for(i in rownames(cm.best)){
    
    best.k <- best.idx <- best.val <- NA

    ## delta area  
    if(i == 'delta.auc.diff'){ # choose K before biggest drop in delta AUC
      best.idx <- max(which.max(cm.all[, i])-1, 1)
      #best.idx <- max(which( cm.all[, i] > delta.auc.sig))
    }
    ## cluster consensus / silhouette score
    if(i == 'cluster.consensus' |  i == 'silhouette.score'){ # choose K with maximal silhouette score
      best.idx <- which.max(cm.all[, i])
    }
    ## cophenetic correlation
    if(i == 'cophenetic.correlation'){ # choose K with maximal cophentic correlation
      best.idx <- which.max( cm.all[, i])
    }
      
    ## fill matrix
    if(!is.na(best.idx)){  
      best.val <- cm.all[best.idx, i]
      best.k <- k.all[best.idx]
      cm.best[i, ] <- c(best.k, best.idx, best.val) %>% as.numeric
    }
  }
  
  ## ############################################
  ## plot
  if(plot){
    
    ## number of rows in plot
    nr <- ceiling((nrow(cm.best)-1)/2)
    
    pdf(glue('{prefix}_cluster_metrics_K{k.min}-K{k.max}.pdf'), 6, 3*nr)
    par(mfrow=c(nr ,2))
    for(i in colnames(cm.all)){

      best.idx <- cm.best[i, 'best.k.idx' ]
      metric <- cm.all[, i]
      
      if(i == 'delta.auc.diff'){
        metric <- cm.all[, 'delta.auc']
      }
      if( !(i %in% c('delta.auc', 'cophenetic.correlation.diff') ) ) {
        plot(1:n.clust,  metric, col='black', pch=20, cex=2, main=i, xaxt='n', xlab='Number of clusters', ylab='score', type='b')
        points(c(1:n.clust)[best.idx], metric[best.idx], col='red', pch=20, cex=2)
        axis(1, at=1:n.clust, labels=k.all)
        legend('topright', legend=glue('k={cm.best[i, "best.k"]}'), bty='n')
      }
    }
    dev.off()
  }
  
  ## #############################################
  ## determine best k
  ## -  pick K according to majority of three metrics
  ## - if all three netrix don't agree, use delta AUC
  ##   as metric
  k.opt <- table(cm.best[ , 'best.k'])
  if(length(k.opt) < 3){
    k.opt <- names(k.opt)[ which.max(k.opt) ] %>% as.numeric %>% sort
  } else {
    k.opt <- cm.best['delta.auc.diff', 'best.k']
  }
  ## ############################################
  ## assemble output
  out <- c()
  out$all.metrics <- cm.all
  out$best.metrics <- cm.best
  out$k.opt <- k.opt
  
  return(out) 
}



## #################################################################
## heatmap of consensus matrix
hm.consensus <- function(M, 
                         membership,
                         cdesc=NULL,
                         method='ward.D2',
                         ...){
  k <- membership %>% as.character %>% as.numeric %>% max

  
  anno <- data.frame(cluster=as.factor(membership))
  if(!is.null(cdesc))
    anno <- data.frame(anno, cdesc)
  
    
  clust.col <- RColorBrewer::brewer.pal(max(k, 3), "Dark2")[1:k]
  names(clust.col) <- unique(anno$cluster)
  
  anno.col <- list(cluster=clust.col)
  if(!is.null(cdesc)){
    n.class <- length(unique(cdesc[, 1]))
    col.cdesc <- unlist( brewer.pal( max(n.class, 3), "Set1") )[1:n.class]
    names(col.cdesc) <- unique( cdesc[, 1] )
      
    anno.col[[2]] <- col.cdesc
    names(anno.col)[2] <- colnames(cdesc)[1]
  }
  
  hm.col <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)
  breaks <- seq(0, 1, length.out = length(hm.col) + 1)
  
  pheatmap(M, symm=T, 
           col=hm.col,
           breaks=breaks,
           annotation_col = anno,
           annotation_row = anno,
           annotation_colors = anno.col,
           clustering_method=method,
           clustering_distance_row=as.dist(1 - M),
           clustering_distance_col=as.dist(1 - M),
                      ...)

}
## ##################################################################
## pca cluster plot
fviz_cluster_plot <- function(data, membership){
  p <- fviz_cluster(list(cluster=membership, data=t(na.omit(data))), palette='Dark2', repel=T, labelsize = 8)
  plot(p)
}

## #################################################################
## consensus clustering:
## calculate cluster consensus m_k
##
clust_cons <- function(M,           ## consensus matrix
                       membership   ## cluster assignment for rows/columns in M
){
  ## number of clusters
  K <- membership %>% as.character %>% as.numeric %>% max
  
  m_k <- sapply(1:K, function(x) {
    m=M[ membership==x, membership==x ]
    m=m[upper.tri(m, diag = F)]
    mean(m)
    } )
  return(m_k)
}


## #################################################################
## consensus clustering
## - area under cdf
auc_cdf_cons <- function(M){

  ## function to calculate CDF
  cdf_cons <- function(c){
    m <- ncol(M)*(ncol(M)-1)/2 
    M <- M[upper.tri(M,diag = F)]
    return( sum(M < c)/m )
  }
  
  m <- ncol(M)*(ncol(M)-1)/2 
  x <- M[upper.tri(M,diag = F)] %>%  sort
    
  auc <- sum(sapply(2:m, function(i) (x[i] - x[i-1])*cdf_cons(x[i]) ))
    
  return(auc)
}

## #################################################################
## 

## #################################################################
## consensus clustering:
## calculate item's consensus
item_cons <- function(M,
                      clust){

  ## TODO  
}


## ##################################################################
## silhoutte scores:
silhouette_scores <- function(cluster.object,
                              data=NULL,
                              plot=T,
                              distance=c('euclidean')
                      ){
  require(glue)
  distance <- match.arg( distance )
  
  method <- attributes(cluster.object)$class
  
  #if(method == 'NMFfitX1')
  if( pmatch('NMF', method, nomatch=0 ) == 1){
    #sil <- silhouette(cluster.object$cluster)
    d <- dist(t(data), method = distance)
  }
  if(method == 'kmeans'){
    d <- dist(t(data), method = distance)
    #d <- daisy(t(data), metric = 'euclidean')
  }
  
  if(method == 'hclust'){
    d <- cluster.object$dist
    #sil <- silhouette(cluster.object$cluster, d)
  }
  
  #sil <- silhouette(cluster.object$cluster, d) %>% sortSilhouette()
  sil <- silhouette(cluster.object$cluster, d)
  
  #if(!is.null(data))#{
  #  rownames(sil) <- colnames(data)[ rownames(sil) %>% as.numeric ]
  #} else {
  #  rownames(sil) <- colnames(d)[ rownames(sil) %>% as.numeric ]
  #}
  if(plot)
    plot(sil, main=glue('Silhouette plot, cluster method: {method}'))
  
  return(sil)
}

## ###################################################################
## cophenetic correlation
## - code snippets copied from gdac cnmf cluster module
cophenetic_correlation <- function(connect.matrix){
  
  # Compute the distance matrix
  dist.matrix <- as.dist(1 - connect.matrix)
  HC <- hclust(dist.matrix, method="average")
  # Update the rho for the index
  dist.coph <- cophenetic(HC)
  
  rho <- signif((cor(dist.matrix, dist.coph)), digits = 4)

  return(rho)
}

## #############################################################
## simulate data by drawing form normal distribution
## export as gct for testing cluster module in PGDAC pipeline
simulate.data <- function(p=100,       # number of features
                          n=40,        # number of samples 
                          k=4          # number of classes, i.e. Gaussian distr with different mean/sd
){
  
  N <- p*n
  Nk <- N/k
  nk <- n/k
  
  comp <- lapply(1:k, function(i) rnorm(Nk, rnorm(1, 0, 3), rnorm(1, 1, 0.1) ) )
  #comp <- lapply(comp, function(x) matrix(x, nrow=p, dimnames=list(paste('feature', 1:p, sep='-'), paste('sample', 1:(n/k), sep=''))) )
  comp <- lapply(comp, function(x) matrix(x, nrow=p))
  
  mat <- Reduce(f = cbind, comp)
  dimnames(mat) <- list(paste('feature', 1:p, sep='-'), paste('sample', 1:n, sep='')) 
  
  rid <- rownames(mat)
  cid <- colnames(mat)
  rdesc <- data.frame(dummy=rep('ABC', p))
  cdesc <- data.frame(class=lapply(1:k, function(x) rep(glue('c{x}'), nk)) %>% unlist)
  rownames(cdesc) <- colnames(mat)
  
  gct <- new('GCT')
  gct@mat <- mat
  gct@cdesc <- cdesc
  gct@rdesc <- rdesc
  gct@rid <- rid
  gct@cid <- cid
  
  write.gct(gct, ofile=glue('proteome-ratio-norm-NArm'), appenddim = F)
  
  pheatmap(mat, annotation_col=cdesc)
}



## ##########################################################################################
## function to determine a range of 'best' number of clusters using the NbClust package
## - does not always work and will fail with:
## 
## "Error in NbClust(data, method = "kmeans", index = "all", min.nc = min.k,  : 
## The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated."
##
## - this happens on matrices without ANY missing vsalues, too.
## - function is not usable in this state.
# best_k <- function(data,
#                    method='kmeans',
#                    distance='euclidean',
#                    min.k=2,
#                    max.k=10,
#                    index='all',
#                    
#                    plot=T){
# 
#   require(NbClust)
# 
#   ## calculate ~30 indices to determine best k
#   index <- NbClust(data, method='kmeans', index='all', min.nc=min.k, max.nc=max.k)
# 
#   ## results
#   index.nc <-table( index$Best.nc )
#   best.nc <- index.nc[ which(index.nc == max(index.nc, na.rm=T)) ]
# 
#   return(best.nc)
# }
# # 
