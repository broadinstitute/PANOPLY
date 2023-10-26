#####################################################################################
kmeans_clustering<-function(clustering_data, kmid, kalt) {
  
  set.seed(2020)
  #center and scale expression matrices
  normmatrix <- t(scale(t(clustering_data), scale=TRUE, center=TRUE)) 
  
  #test different clustering configurations, and use the silhouette index to pick the best one
  #cat("Calculating distance matrix\n")
  dist.mat <- dist(normmatrix)
  threshold <- NULL
  s.index <- NULL
  row <- 1
  for (k in seq(kmid-kalt,kmid+kalt,by=floor(kalt/2))){
    #cat("Testing k=",k,"\n",sep="")
    #k-means clustering
    km.all <- suppressWarnings(kmeans(normmatrix, k,  nstart=5, iter.max=20))
    s <- silhouette(km.all$cluster,dist.mat)
    s.sum <- summary(s)
    threshold[row] <- k
    s.index[row] <- s.sum$avg.width
    row = row+1
  }
  
  #save the best configuration, which is the largest silhouette index 
  my.k <- threshold[s.index==max(s.index)]
  #cat(paste("Choose k=",my.k,"\n",sep=""))
  km.all <- kmeans(normmatrix, my.k, nstart=50, iter.max=20)
  results <- as.data.frame(km.all$cluster)
  colnames(results) <- "clusters"
  
  return(results)
}
