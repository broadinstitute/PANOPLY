

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


