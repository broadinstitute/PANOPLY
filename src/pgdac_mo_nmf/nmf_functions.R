
## ##################################################################
## function to map NMF cluster to levels in 'class.variable'
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
  #data.str <- paste('..', tmp.dir, conf[, 2], sep='/')
  data.str <- paste(tmp.dir, conf[, 2], sep='/')
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
  
  return(list(expr=expr, cdesc=cdesc, rdesc=rdesc, data.str=data.str))

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



## #############################################
## function ...
nmf.plots <- function(ws){
  
  ## ############################
  ## import workspace
  cat('loading ws...')
  load(ws)
  cat('done.\n')
  
  ## ###################################################
  ## class variable/annotation tracks of interest
  class.variable <- opt$class.variable
  variable.other <- strsplit( opt$variable.other, ';' ) %>% unlist
  class.variable <- gsub('\\\'', '', class.variable) 
  variable.other <- gsub('\\\'', '', variable.other) 
  
  ## pheatmap parameters
  cw <- 15
  ch <- 15
  max.val <- 10
  
  
  ## ############################################
  ## save a copy
  #expr.full.org <- expr.full
  expr.org <- expr
  cdesc.org <- cdesc
  rdesc.org <- rdesc
  variable.other.org <- variable.other
  class.variable.org <- class.variable
  
  ## determine 'omics data types
  data.omes <- sub('(^.*?)-.*', '\\1', rownames(expr.org)) %>% unique 
  
  
  ## #################################################
  ##                  colors
  ## #################################################
  class.colors <- opt$class.colors
  color.all <- class.colors %>% strsplit(., '\\|') %>% unlist
  
  ## put into a format suitable for 'pheatmap'
  cdesc.color <- vector('list', length(color.all))
  names(cdesc.color) <- sub('^(.*?)=.*', '\\1', color.all)
  
  color.all <- sub('.*?=','', color.all)
  
  for(i in 1:length(cdesc.color)){
    color.tmp <- color.all[i] %>% strsplit(., ';') %>% unlist %>% strsplit(. , ':')
    color.names <- sapply(color.tmp, function(x)x[1])
    color.tmp <- sapply(color.tmp, function(x)x[2])
    names(color.tmp) <- color.names
    color.tmp <- gsub('\\\'', '', color.tmp) 
    names(color.tmp) <- gsub('\\\'', '', names(color.tmp)) 
    cdesc.color[[i]] <- color.tmp
  }
  cdesc.color <- lapply(cdesc.color, function(x){names(x)=sub('\\n.*', '',names(x));x})
  cdesc.color.org <- cdesc.color
  
  
  ## ################################
  ## parse genes of interest
  gene.column <- opt$gene.column
  genes.of.interest <- opt$genes.of.interest
  if(nchar(genes.of.interest) < 3){
    genes.of.interest <- NULL
  } else{
    genes.of.interest <- strsplit( genes.of.interest, ';' ) %>% unlist %>% unique
  }
  
  
  ## ####################################################
  ##              loop over ranks
  for(rank in rank.top) {
    
    #for(rank in as.character(ranks)){
    #for(rank in as.character(c(4,5,6))){
    dir.create(paste('K', rank, sep='_'))
    setwd(paste('K', rank, sep='_'))
    
    
    ## extract NMF results
    res <- res.rank[[as.character( rank )]]
    
    
    ## ##########################################
    ## silhoutte plots
    pdf(paste('1.0_silhouette_K_', rank, '.pdf', sep=''), 10, 6)
    #par(mfrow=c(1,2))
    plot(rank.sil[[rank]], main=paste('K=', rank, sep=''), col=palette()[1:as.numeric(rank)+1] )
    #plot(rank.sil.random[[rank]], main=paste('Randomized data', sep=''), col=palette()[1:as.numeric(rank)+1])
    dev.off()
    
    
    ## colors
    cdesc.color <- cdesc.color.org
    cdesc <- cdesc.org
    variable.other <- variable.other.org
    class.variable <- class.variable.org
    
    ## data matrix
    expr <- expr.org
    
    
    ## ########################################
    ## add NMF classification to cdesc object
    NMF.basis <- predict(res)
    NMF.consensus <- predict(res, 'consensus')
    
    ## map consensus to basis
    NMF.consensus.map <- NMF.consensus
    for(ii in 1:as.numeric(rank))
      NMF.consensus[ NMF.consensus == ii ] <-  NMF.basis[NMF.consensus == ii]
    
    ## add to cdesc
    cdesc <- data.frame(cdesc, NMF.basis=as.factor(NMF.basis), NMF.consensus=as.factor(NMF.consensus), stringsAsFactors=F)
    
    ## export
    write.table(cdesc, file=paste('clin_anno_nmf.txt'), sep='\t', quote=F, na='', col.names=NA)
    
    ## add to heatmap annotation tracks
    variable.other <- c('NMF.consensus',  'NMF.consensus.mapped', variable.other)
    
    
    
    ## #######################################################################################################
    ## add NMF to color list for pheatmap
    ## - map NMF classes to levels in 'class.variable'
    mb.nmf.map <- data.frame(mb=cdesc[, class.variable], nmf_basis=NMF.basis, nmf_cons=NMF.consensus)
    
    ## map consensus to 'class.variable'
    cons.map <- tapply(mb.nmf.map$mb, mb.nmf.map$nmf_cons, table)
    
    ## relative frequency
    cons.map.rel <- lapply(cons.map, function(x) x/sum(x))
    ## matrix
    cons.map.rel.mat <- matrix(0, nrow=as.numeric(rank), ncol=length(unique(cdesc[, class.variable])), dimnames = list(1:as.numeric(rank), unique(cdesc[, class.variable])))
    for(j in 1:nrow(cons.map.rel.mat))
      cons.map.rel.mat[j, names(cons.map.rel[[j]])] <- cons.map.rel[[j]]
    write.table(cons.map.rel.mat, col.names = NA, sep='\t', file=paste('nmf_vs_', class.variable, '.txt', sep=''), quote=F)
    
    ## pick the maximum
    #cons.map.max <- sapply(cons.map.rel, function(x) {m=max(x); ifelse( sum(x %in% m) > 1, paste(names(x)[which(x %in% m)] , collapse='|'), names(x)[which.max(x)])})
    cons.map.max <- sapply(cons.map.rel, function(x) {
      m=max(x); 
      ifelse( sum(x %in% m) > 1, 
              paste(names(x)[which(x %in% m)] , collapse='|'),
              paste(names(x)[x > 0.3], collapse='|' ) )
    })
    
    ## assign colors 
    NMF.consensus.col <- cdesc.color[[class.variable]][ cons.map.max ]
    names(NMF.consensus.col) <- cons.map.max
    
    ######################################################  
    ## check whether all NMF basis have color
    if( sum(is.na(NMF.consensus.col)) > 0 ){
      
      #idx.tmp <- which( !(levels(NMF.consensus) %in% names(NMF.consensus.col)) )
      idx.tmp <- which( is.na(NMF.consensus.col) )
      
      col.tmp <-  rev( palette() )[ (1:length(idx.tmp) ) +  1] 
      #names(col.tmp) <- idx.tmp
      #NMF.consensus.col <- c(NMF.consensus.col, col.tmp)
      NMF.consensus.col[idx.tmp] <- col.tmp
    }
    cdesc$NMF.consensus.mapped <-  cdesc$NMF.consensus
    cdesc.color$NMF.consensus.mapped <- NMF.consensus.col
    names(cdesc.color$NMF.consensus.mapped) <- 1:length(NMF.consensus.col)
    
    
    
    ## ########################################
    ## order according to nmf clustering
    nmf.ord <-  rownames(cdesc)[cdesc$NMF.consensus %>% as.numeric %>% order ] 
    expr <- expr[, nmf.ord]
    cdesc <- cdesc[nmf.ord, ]
    
    
    
    
    ## ########################################
    ## heatmap coefficients
    #coefmap(res, annCol=cdesc[, c(class.variable, variable.other)], annColors=cdesc.color, filename='1_coefmap.pdf', tracks='consensus', main='Metagenes (consensus)')
    #coefmap(res, annCol=cdesc[, c(class.variable, variable.other)], annColors=cdesc.color, filename='2.0_coefmap.pdf')
    
    # pheatmap version
    H <- res@fit@H[ , nmf.ord ]
    H.norm <- apply(H, 2, function(x)x/max(x))
    pheatmap(H.norm, cluster_row=F, annotation_col=cdesc[, rev(c(class.variable, variable.other))],
             annotation_colors = cdesc.color,
             color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), 
             filename='2.1_coefmap_pheatmap.pdf', cellwidth = cw, cellheight = ch)
    
    ## disable clustering
    #pheatmap(H.norm[, order(cdesc[, 'NMF.consensus'])], cluster_row=F, cluster_cols = F, annotation_col=cdesc[order(cdesc[, 'NMF.consensus']), rev(c(class.variable, variable.other) )],
    pheatmap(H.norm, cluster_row=F, cluster_cols = F, annotation_col=cdesc[ , rev(c(class.variable, variable.other) )],
             annotation_colors = cdesc.color,
             color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), 
             filename='2.2_coefmap_pheatmap_sorted.pdf', cellwidth = cw, cellheight = ch)
    
    
    #pheatmap(H.norm, cluster_row=F, color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), filename='2_coefmap_pheatmap.pdf')
    
    
    ## #######################################
    ## consensus map
    cons <- res@consensus[nmf.ord, nmf.ord]
    #colnames(cons) <- rownames(cons) <- colnames(H)
    
    res2 <- res
    res2@consensus <- cons
    consensusmap(res2, filename=paste('3.0_consensusmap_nrun_', opt$nrun, '.pdf', sep=''))
    
    #cons.ord <- order(NMF.consensus)
    #cons <- cons[cons.ord, cons.ord ]
    
    pheatmap(cons, cluster_row=F, cluster_col=F, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
             annotation_colors=cdesc.color, color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), 
             cellwidth = cw, cellheight = ch,
             filename=paste('3.1_consensusmap_nrun_', opt$nrun,'_pheatmap.pdf', sep=''))
    
    
    ## ##########################################
    ## consensus matrix from bootstrapping
    #cons.bs <- res.cons[[as.character( rank )]]
    #try(
    #  pheatmap(cons.bs, cluster_row=T, cluster_col=T, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
    #         annotation_row = cdesc[ , rev(c(class.variable, variable.other))],
    #         annotation_colors=cdesc.color, color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), 
    #         cellwidth = cw-10, cellheight = ch-10, symm=T, fontsize_row = 5, fontsize_col = 5,
    #         filename=paste('3.2_consensus_matrix_bootstrap_nrun_', nrun.bs,'.pdf', sep=''))
    #    )
    
    ## ordered by NMF cluster assignment
    #cdesc.bs <- cdesc[ , rev(c(class.variable, variable.other))]
    #ord.idx <- order(cdesc.bs$NMF.consensus)
    #cons.bs.ord <- cons.bs[ ord.idx, ord.idx ]
    
    #pheatmap(cons.bs.ord, cluster_row=F, cluster_col=F, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
    #         annotation_row = cdesc[ , rev(c(class.variable, variable.other))],
    #         annotation_colors=cdesc.color, color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), 
    #         cellwidth = cw-10, cellheight = ch-10, symm=T, fontsize_row = 5, fontsize_col = 5,
    #         filename=paste('3.3_consensus_matrix_bootstrap_nrun_', nrun.bs,'_NMF_ordered.pdf', sep=''))
    
    ## ###############################################################
    ##
    ##                   extract NMF features
    ##
    ## ###############################################################
    
    # success <- F
    # pdf('4.0_histogram_feature_scores.pdf')
    # j <- 1
    # while(!success){
    # 
    #   feature.method=feature.methods[[j]]
    #   feature.method.score=feature.methods.scores[[j]]
    #   
    #   ## feature scores
    #   s.i <- featureScore(res, method=feature.method.score)
    #   
    #   hist( unlist(s.i), main=paste("Feature scores (method:", feature.method,")"), xlab='Feature score' )
    # 
    #   ## extract features
    #   s <- extractFeatures(res, method=feature.method)
    # 
    #   if( sum( sapply(s, function(x) sum(is.na(x))) ) == 0 ){
    #     success=T
    #   } else {
    #     warning(paste("\nExtracting NMF features using method '", feature.method, "' did not return features for all basis\n"))    
    #     j <- j + 1
    #     if(j > length(feature.methods)){
    #      break
    #       warning('No success!')
    #     }  
    #   }
    # }
    # dev.off()
    # 
    #success <- F
    ## ###############################################################################
    feature.methods <- list('kim', 'max', 1, 10L)
    feature.methods.scores <- list('kim', 'max', 'max', 'max')
    
    
    ## ##############################################################
    ##
    ##           extract meta feature matrix
    ##
    ## used as ranking for cluster marker selection
    ## ##############################################################
    W <- res@fit@W
    colnames(W) <- 1:ncol(W)
    ## row-normalize
    W.norm <- t(apply(W, 1,  function(x)x/sum(x)))
    
    ## combine up/down: take highest normalized score
    feat.comb <- unique( sub('^(.*-.*)_(up|down)$' , '\\1', rownames(W.norm)) )
    W.norm.comb <- matrix(NA, nrow=length(feat.comb), ncol=ncol(W))
    dimnames(W.norm.comb) <- list(feat.comb, colnames(W))
    
    W.norm.dir <- W.norm
    direction.tmp <- sub('.*_(.*)$', '\\1', rownames(W.norm.dir))
    ## sign the coefficients
    direction.tmp[direction.tmp == 'down' ] <- -1
    direction.tmp[direction.tmp == 'up' ] <- 1
    direction.tmp <- as.numeric(direction.tmp)
    W.norm.dir <- direction.tmp*W.norm.dir
    
    for(ii in feat.comb){
      idx.tmp <- grep(glue("^{ii}_(up|down)$"), rownames(W.norm))
      
      #direction.tmp <- sub('.*_(.*)$', '\\1', rownames(W.norm)[idx.tmp])
      #for(iii in 1:length(idx.tmp))
      
      W.tmp <- W.norm.dir[ idx.tmp, ]
      if(length(idx.tmp) > 1){
        W.norm.comb[ii, ] <- apply(W.tmp, 2, function(x) x[ which.max(abs(x)) ])
      } else {
        W.norm.comb[ii, ] <- W.tmp
      }
      
    }
    ## gct
    #w.comb.rdesc <- data.frame(Type=sub('^(.*?)-.*','\\1', rownames(W.norm.comb)), stringsAsFactors = F)
    if(opt$gene.column %in% colnames(rdesc)){
      w.comb.rdesc <- data.frame(
        Type=sub('^(.*?)-.*','\\1', rownames(W.norm.comb)), 
        geneSymbol=gct.comb@rdesc[rownames(W.norm.comb), opt$gene.column],
        stringsAsFactors = F)
    } else {
      w.comb.rdesc <- data.frame(
        Type=sub('^(.*?)-.*','\\1', rownames(W.norm.comb)), 
        geneSymbol=rownames(W.norm.comb), 
        stringsAsFactors = F)
      
    }
    rownames(w.comb.rdesc) <- rownames(W.norm.comb)
    w.gct <- new('GCT')
    w.gct@mat <- data.matrix(W.norm.comb)
    w.gct@rid <- rownames(W.norm.comb)
    w.gct@cid <- colnames(W.norm.comb)
    w.gct@rdesc <- w.comb.rdesc
    write.gct(w.gct, ofile=glue("matrix_W_combined_signed"))
    
    ## export
    #write.table(W.norm.comb, sep='\t', quote=F,  na='', col.names=NA, file=glue("matrix_W_combined.txt"))
    
    
    ## #############################################################
    ## 
    ##                      extract features
    ##
    ## ##############################################################
    success <- F
    pdf('4.0_histogram_feature_scores.pdf')
    j <- 1
    while(!success){
      
      feature.method=feature.methods[[j]]
      feature.method.score=feature.methods.scores[[j]]
      
      ## feature scores
      s.i <- featureScore(res, method=feature.method.score)
      
      hist( unlist(s.i), main=paste("Feature scores (method:", feature.method,")"), xlab='Feature score' )
      
      ## extract features
      s <- extractFeatures(res, method=feature.method)
      
      if( sum( sapply(s, function(x) sum(is.na(x))) ) < length(s) ){
        success=T
      } else {
        warning(paste("\nExtracting NMF features using method '", feature.method, "' did not return any features\n"))
        j <- j + 1
        if(j > length(feature.methods)){
          break
          warning('No success!')
        }
      }
    }
    dev.off()
    
    
    #W.norm <- W.norm*log2(W.norm)
    #H.norm <- apply(H, 2, function(x)x/max(x))
    
    ## feature scores
    #s.i <- featureScore(res, method=feature.methods.scores[[1]])
    
    ## extract features
    #s <- extractFeatures(res, method=feature.methods[[1]])
    
    #names(s) <- names(s.i) <- paste('features_basis', 1:length(s), sep='_')
    names(s) <- paste(1:length(s))
    
    ## cluster with no features
    rm.idx <- which( sapply(s, function(x) sum(is.na(x))) > 0)
    if(length(rm.idx) > 0){
      s <- s[-c(rm.idx)]
      # s.i <- s.i[-c(rm.idx)]
      
      cons.map <- cons.map[-c(rm.idx)]
      cons.map.max <- cons.map.max[-c(rm.idx)]
    }
    
    #pdf('4.0_histogram_feature_scores.pdf')
    #hist( unlist(s.i), main=paste("Feature scores (method:", feature.methods[[1]],")"), xlab='Feature score' )
    #dev.off()
    
    
    ## ####################################
    ##  annotate the features
    
    ## get accession numbers
    s.acc <- lapply(s, function(x) data.frame(Accession=sub('_up|_down','',rownames(expr.comb)[x]) , id=rownames(expr.comb)[x], stringsAsFactors = F) ) 
    #names(s.acc) <- paste('features_basis', 1:length(s.acc), sep='_')
    
    
    ## features plus scores
    s.acc.scores <- lapply(s.acc, function(x) data.frame(x, Score=s.i[x$id], stringsAsFactors = F))
    
    ## ###################################
    ## upset plot
    upset.mat <- matrix(0, ncol=length(s.acc), nrow=length(unique(unlist(sapply(s.acc, function(x)x$Accession )))), dimnames=list(unique(unlist(sapply(s.acc, function(x)x$Accession ))), names(s.acc)))
    for(ii in names(s.acc))
      upset.mat[s.acc[[ii]]$Accession, ii] <- 1
    
    pdf('4.1_upset_NMF_markers.pdf')
    upset(data.frame(upset.mat), point.size = 4, text.scale = 1.5)
    dev.off()
    
    ## ####################################
    ## gene names
    if(opt$gene.column %in% colnames(rdesc)){
      s.gn <- lapply(s.acc, function(x) rdesc[x$Accession, opt$gene.column] )
      s.gn.red <- s.gn
      
    } else {
      ## unique gene names
      s.gn <- lapply(s.acc, function(x) unique(sub('\\|.*','',sub('.*?\\|(.*).*', '\\1', sub('^CNV-(.*?)\\|.*', '\\1', x$Accession)))))
      ##  redundant gene names
      s.gn.red <- lapply(s.acc, function(x) sub('\\|.*','',sub('.*?\\|(.*).*', '\\1', sub('^CNV-(.*?)\\|.*', '\\1', x$Accession))))
      
    }
    
    ## keep accession
    s.acc2 <- lapply(s.acc, function(x)sub('^(.*?)\\|.*$', '\\1', x$Accession))
    
    ## annotations as data frame
    s.acc.gn <- s.acc
    for(i in 1:length(s.acc)){
      
      tmp <- data.frame(
        Type=sub('-.*','', s.acc[[i]]$Accession),
        Accession=sub('.*?-','', s.acc2[[i]]),
        Direction=sub('.*_(up|down)$','\\1', s.acc[[i]]$id),
        SYMBOL=s.gn.red[[i]],
        Score=s.acc.scores[[i]]$Score,
        stringsAsFactors=F)
      s.acc.gn[[i]] <- tmp[order(tmp$SYMBOL),]
    }
    s.acc.gn.anno <- s.acc.gn
    
    ## ################################################################
    if(opt$gene.column %in% colnames(rdesc)){
      ## add description and enzyme codes
      s.acc.gn.anno <- lapply( s.acc.gn, function(x)
        AnnotationDbi::select(org.Hs.eg.db, keys=x$SYMBOL , column=c( 'GENENAME',  'ENZYME'), keytype='SYMBOL', multiVals='first')
      )
      
      ## join in a single data frame
      for(i in 1:length(s.acc)){
        tmp <- merge(s.acc.gn[[i]], s.acc.gn.anno[[i]], 'SYMBOL' )
        tmp <- tmp[which(!duplicated(apply(tmp, 1, paste, collapse=' '))), ]
        s.acc.gn.anno[[i]] <- tmp
      }
      
      
      ## ##########################################################
      ## extract protein kinases and phosphatases
      s.acc.gn.kinase <- lapply(s.acc.gn.anno, function(x) x[ grep('2\\.7\\.1[0-2]|3\\.1\\.3\\.(16|48)', x$ENZYME), ])
      names(s.acc.gn.kinase) <- paste(names(s.acc.gn.kinase), 'KINASE', sep='_')
      
      s.acc.gn.kinase <- lapply(s.acc.gn.kinase, function(x) x[which(!duplicated(apply(x, 1, paste, collapse=' '))),] )
      
      ## sort
      s.acc.gn.kinase <- lapply(s.acc.gn.kinase, function(x)x[order(x$Score, decreasing=T), ])
      
      ## kinases
      list2env(s.acc.gn.kinase , envir=.GlobalEnv)
      WriteXLS(names(s.acc.gn.kinase), ExcelFileName=paste('NMF_KINASE_PHOSPHATASE_features_.xlsx', sep=''), FreezeRow=1, FreezeCol=1, SheetNames=make.unique(substr(names(s.acc.gn.kinase), 1, 20)), row.names=F, BoldHeaderRow=T, AutoFilter=T)
      
    }
    
    ## sort
    s.acc.gn.anno <- lapply(s.acc.gn.anno, function(x)x[order(x$Score, decreasing=T), ])
    
    
    ## ##########################################################
    ## export as txt
    ## all features
    list2env(s.acc.gn.anno , envir=.GlobalEnv)
    WriteXLS(names(s.acc.gn.anno), ExcelFileName=paste('NMF_features_N_', length(unique(unlist(s.acc2))),'.xlsx', sep=''), FreezeRow=1, FreezeCol=1, SheetNames=make.unique(substr(names(s.acc.gn.anno), 1, 20)), row.names=F, BoldHeaderRow=T, AutoFilter=T)
    
    
    
    
    ## #####################################################################################
    ##
    ##               contribution of each data type to NMF coefficients
    ##
    ## ######################################################################################
    data.types <- lapply(s.acc, function(x) table(sub('-.*','', x$Accession)))
    tmp <- lapply(data.types, names)
    tmp <- unique(unlist(tmp))
    data.types.mat <- matrix(0, nrow=length(data.types), ncol=length(tmp), dimnames=list(names(s), tmp) )
    for(i in rownames(data.types.mat))
      data.types.mat[i, names(data.types[[i]])] <- data.types[[i]]
    
    ## add mapping to class variable
    #idx.tmp <- match( cons.map.max, rownames(data.types.mat))
    
    idx.tmp <- match( names(cons.map), rownames(data.types.mat))
    
    rownames(data.types.mat)[idx.tmp] <- paste(rownames(data.types.mat)[idx.tmp], ' (', cons.map.max, ')', sep='' )
    
    
    col <- RColorBrewer::brewer.pal(length(tmp), "Set2")[1:length(tmp)]
    
    
    ## #############
    ## barplot
    names(col) <- tmp
    pdf('4.2_barplot_features_per_basis2.pdf')
    fancyBarplot(t(data.types.mat), srt=30, xlab='NMF basis', col=col, main='Proteogenomic features', ylab='No. features')
    legend('top', legend=names(col), fill=col, bty='n')
    legend('topright', legend=paste('total #:', length(unique(unlist( s.acc2 )))), bty='n')
    dev.off()
    
    ## ##############
    ## table
    m <- data.types.mat
    pheatmap(m, cluster_rows=F, cluster_cols=F, display_numbers=T, legend=F, number_format='%.0f', cellwidth=30, cellheight=30, fontsize=15, number_color='grey20', main='NMF features', scale='none', col='white', border_color='blue', filename='4.3_matrix_features_per_NMF.consensus.pdf')
    
    
    ## ###########################################################
    ## plot the actual expression values for the extracted features
    #for(i in 1:length(s.gn)){
    cc <- 1
    for(i in names(s.gn)){
      
      
      m <- expr[ s.acc[[i]]$Accession, ]
      
      if( !is.null(nrow(m))){
        
        if(is.null(max.val)){
          m.max <- ceiling(max(abs(m), na.rm=T))
        } else {
          m.max <- max.val
          m[m > m.max] <- m.max
          m[m < -m.max] <- -m.max
        }
        #pheatmap(m[, order(cdesc$NMF.consensus)], scale='none', fontsize_row=4, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
        # pheatmap(m, scale='none', fontsize_row=4, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
        #         annotation_colors=cdesc.color, filename=paste('5.',cc-1,'_heatmap_expression_basis_', i, '.pdf', sep=''),  
        #         breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), 
        #         main=names(NMF.consensus.col)[i], cluster_cols = F, height=8, width=8)
        
        ## #####################################
        ## complexheatmap
        col.hm <- colorRamp2(seq(-m.max, m.max, length.out=11), rev(brewer.pal (11, "RdBu")))
        cdesc.ha <- HeatmapAnnotation(df=cdesc[ , rev(c(class.variable, variable.other))], col=cdesc.color)
        #hm <- Heatmap(m, col=rev(brewer.pal (11, "RdBu")), 
        hm <- Heatmap(m, col=col.hm, 
                      
                      cluster_columns = F,
                      top_annotation = cdesc.ha, 
                      split = sub('-.*','',rownames(m)),
                      row_dend_side = 'right',
                      name='NMF features',
                      show_row_names = F)
        ## plot
        pdf(paste('5.',cc-1,'_ComplexHeatmap_expression_basis_', i, '.pdf', sep=''))
        draw(hm)
        for(an in rev(c(class.variable, variable.other)) ) {
          decorate_annotation(an, {
            # annotation names on the right
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
            # annotation names on the left
            #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
          })
        }
        dev.off()
        
        cc <- cc +1
        
      } # end !is.null(nrow(m))
      
    }
    
    
    ## ############################################################
    ##
    ##                plot ALL markers
    ##
    feat.all <- unique(unlist(s.acc2)) %>% sort
    rdesc.feat <- matrix(0, nrow=length(feat.all), ncol=length(s.acc), dimnames=list( feat.all, names(s.acc) ))
    for(i in 1:length(s.acc))
      rdesc.feat[ s.acc[[i]]$Accession, i] <- 1
    rdesc.feat <- data.frame(rdesc.feat)
    
    m <- expr[ feat.all, ]
    
    if(is.null(max.val)){
      
      m.max <- ceiling(max(abs(m), na.rm=T))
      
    } else {
      
      m.max <- max.val
      m[m > m.max] <- m.max
      m[m < -m.max] <- -m.max
      
    }
    #pheatmap(m, scale='none', fontsize_row=0, annotation_row = rdesc.feat, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, 
    #         clustering_distance_col='euclidean', clustering_method='ward.D2', filename='6.0_heatmap_ALL_features.pdf', show_rownames=F, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cellwidth = cw)
    pheatmap(m, scale='none', fontsize_row=0, annotation_row = rdesc.feat, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, 
             cluster_cols = F, filename='6.0_heatmap_ALL_features.pdf', show_rownames=F, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cellwidth = cw)
    
    
    
    ## #########################################################
    ##  complex heatmap
    # reorder according to 'omics type
    #cdesc.ord <- cdesc[ , rev(c(class.variable, variable.other))]
    #ord.idx <- order(cdesc.ord$NMF.consensus)
    
    ## column annotations
    #cdesc.ha <- HeatmapAnnotation(df=cdesc[ ord.idx, rev(c(class.variable, variable.other))], col=cdesc.color)
    cdesc.ha <- HeatmapAnnotation(df=cdesc[ , rev(c(class.variable, variable.other))], col=cdesc.color)
    
    #m <- m[, colnames(m)[ord.idx]]
    col.hm <- colorRamp2(seq(-m.max, m.max, length.out=11), rev(brewer.pal (11, "RdBu")))
    
    hm <- Heatmap(m, col=col.hm, 
                  cluster_columns = F,
                  top_annotation = cdesc.ha, 
                  split = sub('-.*','',rownames(m)),
                  row_dend_side = 'right',
                  name='NMF features',
                  show_row_names = F)
    ## plot
    pdf('6.1_ComplexHeatmap_ALL_features.pdf')
    draw(hm)
    for(an in rev(c(class.variable, variable.other)) ) {
      decorate_annotation(an, {
        # annotation names on the right
        grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
        # annotation names on the left
        #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
      })
    }
    dev.off()
    
    ## ####################################################
    ## z-scored matrix
    hm.scaled <- Heatmap(scale(m), col=rev(brewer.pal (11, "RdBu")), 
                         cluster_columns = F,
                         top_annotation = cdesc.ha, 
                         split = sub('-.*','',rownames(m)),
                         row_dend_side = 'right',
                         name='NMF features',
                         show_row_names = F)
    ## plot
    pdf('6.1_ComplexHeatmap_ALL_features_z-score.pdf')
    draw(hm.scaled)
    for(an in rev(c(class.variable, variable.other)) ) {
      decorate_annotation(an, {
        # annotation names on the right
        grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
        # annotation names on the left
        #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
      })
    }
    dev.off()
    
    
    
    ## ############################################################
    ## plot ALL markersbut without CNV
    #if(length( grep('CNV',names(data.str)) ) > 0 ){
    if(length( grep('^CNV-', rownames(m) )) > 0 ){
      
      #if('CNV' %in% names(data.str)){
      m.org <- m
      #m <- expr.org[ unique( unlist(s.acc)[-grep('CNV',  unlist(s.acc))] ), ]
      cnv.idx <- grep('CNV-', rownames(m.org))
      
      if(length(cnv.idx) > 0){
        m <- m[-cnv.idx, ]
        
        if(is.null(max.val)){
          m.max <- ceiling(max(abs(m), na.rm=T))
        } else {
          m.max <- max.val
          m[m > m.max] <- m.max
          m[m < -m.max] <- -m.max
        }
        #pheatmap(m, scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, clustering_distance_col='euclidean', clustering_method='ward.D2', filename='6.2_heatmap_ALL_features_no_CNV.pdf', show_rownames=F,breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
        pheatmap(m, scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], 
                 annotation_colors=cdesc.color, cluster_cols = F, filename='6.2_heatmap_ALL_features_no_CNV.pdf', show_rownames=F,breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
        
        
        ## #####################################
        ## complexheatmap
        cdesc.ha <- HeatmapAnnotation(df=cdesc[ , rev(c(class.variable, variable.other))], col=cdesc.color)
        col.hm <- colorRamp2(seq(-m.max, m.max, length.out=11), rev(brewer.pal (11, "RdBu")))
        
        hm <- Heatmap(m, col=col.hm, 
                      cluster_columns = F,
                      top_annotation = cdesc.ha, 
                      split = sub('-.*','',rownames(m)),
                      row_dend_side = 'right',
                      name='NMF features',
                      show_row_names = F)
        ## plot
        pdf('6.2_ComplexHeatmap_ALL_features_no_CNV.pdf')
        draw(hm)
        for(an in rev(c(class.variable, variable.other)) ) {
          decorate_annotation(an, {
            # annotation names on the right
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
            # annotation names on the left
            #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
          })
        }
        dev.off()
        
        
        
        
        
        ## ##############################################################################
        ## CNV only
        #m <- expr.org[ unlist(s.acc)[grep('CNV',  unlist(s.acc))] , ]
        m <-m.org[cnv.idx, ]
        m[m > 1] <- 1
        m[m < -1] <- -1
        m.max <- ceiling(max(abs(m), na.rm=T))
        #pheatmap(m, scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color, clustering_distance_col='euclidean', clustering_method='ward.D2', filename='6.2_heatmap_ALL_features_CNV_only.pdf', show_rownames=F,breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
        pheatmap(m, scale='none', fontsize_row=3, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color,
                 cluster_cols = F, filename='6.2_heatmap_ALL_features_CNV_only.pdf', show_rownames=T,breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
        
        
        
        #m <- m[, colnames(m)[ord.idx]]
        
        ## #####################################
        ## complexheatmap
        cdesc.ha <- HeatmapAnnotation(df=cdesc[ , rev(c(class.variable, variable.other))], col=cdesc.color)
        col.hm <- colorRamp2(seq(-m.max, m.max, length.out=11), rev(brewer.pal (11, "RdBu")))
        
        hm <- Heatmap(m, col=col.hm, 
                      cluster_columns = F,
                      top_annotation = cdesc.ha, 
                      split = sub('-.*','',rownames(m)),
                      row_dend_side = 'right',
                      name='NMF features',
                      show_row_names = F)
        ## plot
        pdf('6.2_ComplexHeatmap_ALL_features_CNV_only.pdf')
        draw(hm)
        for(an in rev(c(class.variable, variable.other)) ) {
          decorate_annotation(an, {
            # annotation names on the right
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
            # annotation names on the left
            #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
          })
        }
        dev.off()
        
      }
    }
    
    
    
    ## #############################################################################################
    ##
    ##                        positive conrols / genes of interest
    ##
    ## #############################################################################################
    if(!is.null(genes.of.interest))
    {
      
      # pdf('barchart_positive_controls.pdf', 3, 3)
      for(i in genes.of.interest){
        
        ## extract features
        feat <- lapply(s.acc.gn, function(x){ 
          idx=grep(i, x$SYMBOL, value=F, ignore.case=F)
          paste(x$Type[idx], x$Accession[idx], sep='-' )
        })
        
        if(sum(sapply(feat, length)) > 0){
          
          feat.expr <- expr.org[unlist(feat), ]
          
          ## order: subtype + NMF
          ord.idx <- with(cdesc, order( NMF.consensus))
          
          
          ## if there is only a single row
          if(is.null(dim(feat.expr))){
            feat.expr <- feat.expr[ ord.idx ]
            tmp <- feat.expr
            
            dim(tmp) <- c(1, length(feat.expr))
            colnames(tmp) <-  names(feat.expr)
            rownames(tmp) <- unlist(feat)
            m <- tmp
            
            ## heatmap
            m.max <- ceiling(max(abs(m), na.rm=T))
            pheatmap( m,  annotation_col=cdesc[, rev(c(class.variable, variable.other) )], annotation_colors=cdesc.color, cluster_col=F, gaps_col=cumsum(table(cdesc[, 2])), main=paste('GOI:', i), filename=paste('7_heatmap_gene_of_interest_', sub('\\[.*|\\$','',i), '.pdf', sep=''), cellwidth=cw, cellheight=ch, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cluster_row=F)
            
          } else {
            m <- feat.expr[, ord.idx]
            
            ## heatmap
            m.max <-  ceiling(max(abs(m), na.rm=T))
            pheatmap( m,  annotation_col=cdesc[ , rev(c(class.variable, variable.other) )],annotation_colors=cdesc.color, cluster_col=F, gaps_col=cumsum(table(cdesc[, 'NMF.consensus'])),main=paste('GOI:', i), filename=paste('7_heatmap_gene_of_interest_', sub('\\[.*|\\$','',i), '.pdf', sep=''), cellwidth=cw, cellheight=ch, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")))
          }
          
          feat.n <- unlist(lapply(feat, length))
          
          ## barplot
          pdf(paste('8_barplot_',sub('\\[.*|\\$','',i),'_features.pdf', sep=''), 3, 3)
          fancyBarplot(feat.n, srt=45, xlab='NMF basis', names=1:length(s.acc), main=paste(i), ylab='No. features')
          dev.off()
        }
      }
      # 
    }
    
    ## #######################################################################
    ##
    ##                      compare to PCA
    ##
    ## #######################################################################
    pch.vec <- c(15, 19, 17, 18, 8, 4, 3, 2, 1, 5:7, 9:14, 16, 20:25 )
    
    
    ## ############################################
    ## PCA on entire matrix
    pca <- prcomp( t(expr.org) )
    
    pc1 <- pca$x[,1]
    pc2 <- pca$x[,2]
    
    ## correct colors
    col.tmp <- cdesc.color[[class.variable]]
    col <- as.character(cdesc[,class.variable])
    #names(col) <- cdesc[,1]
    col <- sapply(col, function(i) col.tmp[i])
    
    ## plot
    pdf(paste('9.0_PCA.pdf'))
    #plot(pc1, pc2, col=col,  cex=2, main=paste('PCA: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')), xlab='PC 1', ylab='PC2', pch=pch.vec[ kmeans(cbind(pc1, pc2), rank)$cluster] )
    plot(pc1, pc2, col=col,  cex=2, main=paste('PCA: ', paste(data.omes, collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')), xlab='PC 1', ylab='PC2', pch=20 )
    pointLabel(pc1, pc2, labels=cdesc[,1], col=col, offset=50, method='SANN', cex=.5)
    legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
    
    
    #  fancyPlot(cdesc$CSresponseRank, pc1, xlab='CSresponseRank', ylab='Principle component 1', reg='linear', boxplots = 'n', cor = T, main=paste('PCA: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')))
    #  pointLabel(cdesc$CSresponseRank, pc1, labels=cdesc[,1], col=col, offset=50, method='SANN', cex=.5)
    
    dev.off()
    
    
    
    
    
    ## #############################################
    ## PCA on NMF features
    expr.nmf <- expr[feat.all, ]
    pca.nmf <- prcomp( t(expr.nmf) )
    
    
    pc1 <- pca.nmf$x[,1]
    pc2 <- pca.nmf$x[,2]
    
    col.tmp <- cdesc.color[[class.variable]]
    col <- as.character(cdesc[,class.variable])
    #names(col) <- cdesc[,1]
    col <- sapply(col, function(i) col.tmp[i])
    
    ## correct colors
    ## plot
    pdf(paste('9.1_PCA_NMF_features.pdf'))
    #plot(pc1, pc2, col=col,  cex=2, main=paste('PCA on NMF features: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')), xlab='PC 1', ylab='PC2', pch=pch.vec[ kmeans(cbind(pc1, pc2), rank)$cluster])
    plot(pc1, pc2, col=col,  cex=2, main=paste('PCA on NMF features: ', paste(data.omes, collapse=', '), '\n', paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')), xlab='PC 1', ylab='PC2', pch=20)
    pointLabel(pc1, pc2, labels=cdesc[,1], col=col, offset=20, method='SANN', cex=.5)
    legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
    
    # fancyPlot(cdesc$CSresponseRank, pc1, xlab='CSresponseRank', ylab='Principle component 1', reg='linear', boxplots = 'n', cor = T, main=paste('PCA on NMF features: ', paste(names(data.str), collapse=', '), '\n', paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')) )
    #  pointLabel(cdesc$CSresponseRank, pc1, labels=cdesc[,1], col=col, offset=50, method='SANN', cex=.5)
    
    dev.off()
    
    
    ## ################################################################
    ##
    ##   T-SNE plot of NMF markers
    ##
    ## ################################################################
    ## tsne
    tsne <- try(Rtsne(expr.nmf, check_duplicates=F))
    if(!class(tsne) == 'try-error'){
      
      ## tsne components
      tsneY = tsne$Y
      
      ## #################################
      ## color according to kmeans clustering
      tsne.clust.col=kmeans(tsneY, length(data.omes))$cluster
      
      ## ##################################
      ## color according to data type
      tsne.dt.col=factor(sub('-.*','',rownames(expr.nmf)))
      
      ## ################################
      ## color according to NMF consensus
      tsne.nmf.col <- rownames(expr.nmf)
      names(tsne.nmf.col) <- tsne.nmf.col
      for(i in 1:length(s.acc))
        tsne.nmf.col[ s.acc[[i]]$Accession] <- NMF.consensus.col[ i ] %>% as.character #unlist(NMF.basis.col[as.character(i)])
      
      ## just the NMF class
      tsne.nmf.class <- rownames(expr.nmf)
      names(tsne.nmf.class) <- tsne.nmf.class
      for(i in 1:length(s.acc))
        tsne.nmf.class[s.acc[[i]]$Accession ] <- paste(i,  cons.map.max[i], sep=' ')
      
      
      
      d <- data.frame(
        tsne_1=tsneY[,1],
        tsne_2=tsneY[,2],
        Data_type=sub('-.*', '', rownames(expr.nmf)),
        Feature=rownames(expr.nmf),
        NMF_col=tsne.nmf.col,
        NMF_class=tsne.nmf.class
      )
      
      
      ##plot_ly(d, x=~tsne_1, y=~tsne_2, mode='markers', text=~Feature, color=~NMF_class)
      
      symbs <- c('circle', 'x', 'o', 'square', 'diamond', 'triangle-down', 'pentagon')
      
      #col.tmp <- unique(d$NMF_col)
      #names(col.tmp) <- unique(d$NMF_class)
      col.tmp <- NMF.consensus.col
      names(col.tmp) <- paste(1:length(col.tmp), names(col.tmp ))
      
      symb.tmp <- symbs[ 1:length(data.omes) ]
      names(symb.tmp) <- unique(d$Data_type)
      
      ##
      dir.create('tsne')
      save(d, col.tmp, symb.tmp, file='tsne/tsne.RData')
      
      
      ## ####################################################
      ## Rmarkdown
      rmd <- c(paste("## T-SNE plot of NMF markers (n=",nrow(d),")\n
                     ### Colored by NMF classification:\n
                     ```{r echo=F}
                     load('tsne.RData')
                     library(plotly)
                     plot_ly(d, x=~tsne_1, y=~tsne_2, mode='markers', type='scatter', text=~Feature, color=~NMF_class, colors=as.character(col.tmp))\n
                     ```\n
                     ### Colored by data type:\n
                     ```{r echo=F}
                     plot_ly(d, x=~tsne_1, y=~tsne_2, mode='markers', type='scatter', text=~Feature, color=~Data_type)\n
                     ```
                     \n
                     ", sep=''))
      writeLines(rmd, con='tsne/tsne.rmd')
      try(rmarkdown::render('tsne/tsne.rmd'))
      
    } ## ens tsne
    
    setwd('..')  
  } 
} ## end nmf.plot()


