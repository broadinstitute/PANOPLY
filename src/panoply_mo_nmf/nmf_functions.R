#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

###########################################
## returns 'topn.rank' no. of clusters
## with max. 'metric'-scores
GetBestRank <- function(metric,      ## vector of a metric to assess clustering
                        topn.rank=1, ## No. of topN ranks
                        exclude_2=F,  ## should rank=2 be excluded?
                        rel.inc=1e-3
                        ){

  ## https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
  # ## Answer #37
  # localMaxima <- function(x) {
  #   # Use -Inf instead if x is numeric (non-integer)
  #   y <- diff(c(-.Machine$integer.max, x)) > 0L
  #   rle(y)$lengths
  #   y <- cumsum(rle(y)$lengths)
  #   y <- y[seq.int(1L, length(y), 2L)]
  #   if (x[[1]] == x[[2]]) {
  #     y <- y[-1]
  #   }
  #   y
  # }
  #

  if(exclude_2){
     if('2' %in% names(metric))
       metric <- metric[-which(names(metric) == '2')]
  }

  rank.top <- names(metric)[1]

  success <- F
  i <- 1
  while(!success){

    #cat(i,' - rank top:', rank.top,'\n')

    top_k <- metric[rank.top]
    this_k <- metric[i]

    if(i == length(metric)){
      if(this_k > top_k){
        inc_perc <-(this_k - top_k)/abs(top_k)
        if(inc_perc > rel.inc)
          rank.top <- names(metric)[i]
      }
      #rank.top <- names(metric)[i]
      success <- T

    } else {

       next_k <- metric[i+1]

      if(next_k < this_k){
          i <- i+1
        }
      if(next_k >= this_k){

        if(next_k > top_k){

          inc_perc <-(next_k - top_k)/abs(top_k)
          #cat('diff:', inc_perc, '\n')
          if(inc_perc < rel.inc){
            #rank.top <- names(metric)[i]
            success <- T
          } else {
            rank.top <- names(metric)[i + 1]
            i <- i + 1

            }
        } else {
          i <- i + 1
        }
      }
    }
  }

  return(rank.top)
}

##################################################
##  calculate fisher p of levels of 'class.vec'
##  in each level of 'clust.vec'
CalcClustEnrich <- function(clust.vec,  ## vector of cluster labels
                            class.vec   ## vector of class labels, same order and length as 'clust'

){

  clust.vec <- as.character(clust.vec)
  class.vec <- as.character(class.vec)

  res.per.clust <- matrix(0, nrow=length(unique(clust.vec)), ncol=length(unique(class.vec)),
                          dimnames=list(sort(unique(clust.vec)), unique(class.vec)))

  ## loop over clusters
  for(clust in sort(unique(clust.vec))){

    ## in cluster
    class.vec.clust <- class.vec[ which(clust.vec == clust)]
   
    ## not in cluster
    class.vec.not.clust <- class.vec[ which(clust.vec != clust)]
    
    ## loop over class levels
    for(class.level in na.omit(unique(class.vec))){
      
      x11 <- sum(class.vec.clust %in% class.level)
      x12 <- sum(!(class.vec.clust %in% class.level))
      x21 <- sum( class.vec.not.clust %in%  class.level )
      x22 <- sum(!(class.vec.not.clust %in% class.level ))
        
      #n <- length(class.vec.clust)
      #p <- sum( (class.vec %in% class.level)/length(class.vec) )
    
      
      # Fisher's p
      res.per.clust[clust, class.level] <- fisher.test( rbind(c(x11, x12), c(x21, x22)), alternative = 'greater')$p.value
      
      }
  }
  return(res.per.clust)
}

###########################################
CalcClustEnrichBinom <- function(clust.vec,  ## vector of cluster labels
                            class.vec   ## vector of class labels, same order and length as 'clust'
                            
){
  
  clust.vec <- as.character(clust.vec)
  class.vec <- as.character(class.vec)
  
  res.per.clust <- matrix(0, nrow=length(unique(clust.vec)), ncol=length(unique(class.vec)),
                          dimnames=list(sort(unique(clust.vec)), unique(class.vec)))
  
  ## loop over clusters
  for(clust in sort(unique(clust.vec))){
    
    class.vec.tmp <- class.vec[ which(clust.vec == clust)]
    
    ## loop over class levels
    for(class.level in na.omit(unique(class.vec))){
      x <- sum(class.vec.tmp %in% class.level)
      n <- length(class.vec.tmp)
      p <- sum( (class.vec %in% class.level)/length(class.vec) )
      
      # binomial p
      res.per.clust[clust, class.level] <- binom.test(x, n, p, alternative = 'greater')$p.value
    }
  }
  return(res.per.clust)
}

########################################################
## map clusters to annotations
MapCalcClustEnrich <- function(clust.vec,  ## vector of cluster labels
                               class.vec,   ## vector of class labels, same order and length as 'clust'
                               p.max
){
  
  enrich <- CalcClustEnrich(class.vec, clust.vec)
  
  map=apply(enrich , 2, function(xx) paste(names(xx)[xx < p.max], sep='|'))
  map.names <- names(map)
  ## if multiple mappings were significant,
  ## pick the most significant mapping
  map <- lapply(names(map), function(y){ 
    if(length(map[[y]])>1){ 
      return(map[[y]][ which.min( enrich[map[[y]], y]) ])
    } else {
      return(map[[y]])
    }
  })
  names(map) <- map.names
  map=unlist(map)
  map.names=names(map)
  map=make.unique(map, sep = '')
  names(map)=map.names
  
  for(cc in as.character(unique(clust.vec)))
    clust.vec[which(clust.vec == cc)] <- ifelse(!is.na(map[cc]), unlist(map[cc]), paste0('c', cc)  )
  
  return(clust.vec)
}


## #################################################
## parse.colors
##
## function that assigns colors to levels of
## 'class.variable' and 'variable.other'
## opt        - object returned by function 'parse_args'
## blank      - character, levels with "" or NA are set to
##              this string
## blank.col  - color used for levels set to 'blank'
## #################################################
parse.colors <- function(opt, cdesc, blank='N/A', blank.col='white'){

  ## user defined colors
  if(!is.na(opt$cat_colors)){
        cat_colors <- opt$cat_colors
        color.all <- cat_colors %>% strsplit(., '\\|') %>% unlist
      
        ## check dups
        color.all <- color.all[ !duplicated(sub('^(.*?)=.*', '\\1', color.all)) ]
      
        ## put into a format suitable for 'pheatmap'
        cdesc.color <- vector('list', length(color.all))
        names(cdesc.color) <- sub('^(.*?)=.*', '\\1', color.all)
      
        color.all <- sub('.*?=','', color.all)
      
        ## check whether in meta data
        keep.idx <- names(cdesc.color) %in% colnames(cdesc)
        cdesc.color <- cdesc.color[ keep.idx ]
        color.all <- color.all[ keep.idx ]
        
        ## assign levels to colors
        for(i in 1:length(cdesc.color)){
      
          cat('parsing colors for: ', names(cdesc.color)[i], ' ...')
          
          ## cdesc column
          cdesc.column <- names(cdesc.color)[i]
      
          
          ## present in cdesc?
          if(!cdesc.column %in% colnames(cdesc)){
            
            cdesc.color <- lapply(cdesc.color, function(x){names(x)=sub('\\n.*', '',names(x));x})
            #break;
          
            } else {
              cdesc.levels <- cdesc[, cdesc.column] %>% unique
              cdesc.levels[ nchar(as.character(cdesc.levels)) == 0 | is.na(cdesc.levels) ] <- blank
      
              ## extract colors specified for levels of cdesc.column
              color.tmp <- color.all[i] %>% strsplit(., ';') %>% unlist %>% strsplit(. , ':')
              color.names <- sapply(color.tmp, function(x)x[1])
              color.tmp <- sapply(color.tmp, function(x)x[2])
           
                 names(color.tmp) <- color.names
              color.tmp <- gsub('\\\'', '', color.tmp)
              names(color.tmp) <- gsub('\\\'', '', names(color.tmp))
      
              ## check whether all levels were assigned
              ## fill with random colors
              if(sum(!cdesc.levels %in% names(color.tmp)) > 0){
                color.tmp.bck <- color.tmp
                color.tmp <- cdesc.levels
                names(color.tmp) <- cdesc.levels
                color.tmp[names(color.tmp.bck)] <- color.tmp.bck
      
                idx <- which(!cdesc.levels %in% names(color.tmp.bck))
                color.tmp[idx] <- alphabet()[1:length(idx)]
      
                if(blank %in% names(color.tmp))
                  color.tmp[blank] <- blank.col
              }
              cdesc.color[[i]] <- color.tmp
            } ## end else
             #cdesc.color <- lapply(cdesc.color, function(x){names(x)=sub('\\n.*', '',names(x));x})
          
          cat('done\n')
          } # end for loop
        
        cdesc.color <- lapply(cdesc.color, function(x){names(x)=sub('\\n.*', '',names(x));x})
        
        } else { ## end !is.na(opt$cat_colors) (user defined colors)
          
        ###################################
        ## pick colors automatically
          
          ## categories specified
          cdesc.cat <- opt$cat_anno %>% strsplit(., ';') %>% unlist %>% unique
          cdesc.cat <- cdesc.cat[cdesc.cat %in% colnames(cdesc) ]
          
          ## color palettes to pick from
          my.brewer.pal.info <- brewer.pal.info[ c('Set1', 'Set2', 'Set3', 'Paired', 'Dark2', 'Accent'), ]  %>% rownames_to_column('pal') %>% as.tibble
          
          ## list
          cdesc.color <- vector('list', length(cdesc.cat))
          names(cdesc.color) <- cdesc.cat
          
          ## loop over annotation columns
          for(i in names(cdesc.color)){
            
            ## get levels in current category
            cat.lev <- cdesc[, i] %>% unlist
            if(sum(is.na(cat.lev)) > 0){
              cat.lev[which(is.na(cat.lev))] <- blank.anno
            }
            cat.lev <- unique(cat.lev)
            cat.numb <- length(cat.lev)
            
            ## if there are fewer than 12 categories
            if(cat.numb <= 12){
              pal.tmp <- dplyr::filter(my.brewer.pal.info, maxcolors >= cat.numb)
              ## pick randomly
              pick <- base::sample(1:nrow(pal.tmp), 1)
              
              ## assign colors to levels
              color.tmp <- brewer.pal(n=cat.numb, name=unlist(pal.tmp[pick, 'pal']))
              names(color.tmp) <- cat.lev
              
              ## colors for blanks
              if(blank %in% names(color.tmp))
                color.tmp[blank] <- blank.col
              
              cdesc.color[[i]] <- color.tmp
            }
          }
        }
  return(cdesc.color)
}


## #############################################
## extract tarbar and read configuration file
prepare.data.sets <- function( tar.file, tmp.dir){

  ## #################################
  ## extract tar ball
  if(!dir.exists(tmp.dir))
    dir.create(tmp.dir)
  cat('\nextracting archive...')
  untar(tar.file, exdir=tmp.dir)
  cat('done!\n')

  ##  import config file
  conf <- read.delim(paste( tmp.dir, 'nmf.conf', sep='/'), row.names = NULL, stringsAsFactors = F, header=F)

  data.str <- paste(tmp.dir, conf[, 2], sep='/')
  names(data.str) <- conf[,1]

  return(data.str)
}

## #############################################
## import and merge data tables
#import.data.sets <- function(data.str, zscore.cnv=F){
import.data.sets <- function(tar.file, tmp.dir, zscore.cnv=F){

  ##############################################
  ## helper function to harmonize column
  ## data types of columns in df2 having the
  ## same name in df1 and df2
  ##
  ## Found here (second answer):
  ## https://stackoverflow.com/questions/49215193/r-error-cant-join-on-because-of-incompatible-types
  matchColClasses <- function(df1, df2) {

    sharedColNames <- names(df1)[names(df1) %in% names(df2)]
    sharedColTypes <- sapply(df1[,sharedColNames], class)

    for (n in sharedColNames) {
      class(df2[, n]) <- sharedColTypes[n]
    }

    return(df2)
  }


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
      cdesc[] <- lapply(cdesc, function(x){if(sum(is.na(as.numeric(sub('\\|.*', '', x)))) == 0){as.numeric(sub('\\|.*', '', x));} else {x} })

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
        cdesc.tmp[] <- lapply(cdesc.tmp, function(x){if(sum(is.na(as.numeric(sub('\\|.*', '', x)))) == 0){as.numeric(sub('\\|.*', '', x));} else {x} })

        cdesc.tmp <-  matchColClasses(cdesc, cdesc.tmp)

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


## ##############################################
## filter datasets
filter.datasets <- function(gct.comb,
                            sd_filt,
                            mode=c('global', 'separate', 'equal', 'none')
                            ){

  mode <- match.arg(mode)
  
  if(mode == 'none') return(gct.comb)
 
  ## ###################
  ## calculate SD accross samples
  sd.expr <- apply(gct.comb@mat, 1, sd, na.rm=T)
  
  ## data types
  data.type <- sub('^(.*?)-.*','\\1', gct.comb@rdesc$Data.Type.ID)
  names(data.type) <- rownames(gct.comb@rdesc)
    
  #######################
  ## helper function
  sd_filt_func <- function(sd.vec, perc){
    sd.perc <- quantile(sd.vec, c(perc))
    keep <- which( sd.vec > sd.perc)
    return(names(keep))
  }
  
  ## ###################
  ## global filter
  if(mode == 'global'){
    
    if(sd_filt > 0){
      
      cat("\napplying SD-filter:", mode, '\n')
      sd.keep <- sd_filt_func(sd.expr, sd_filt)
      
      #cat("removed", nrow(gct.comb@mat)-length(sd.keep), "features with SD<=", round(sd.perc, 2), "(", names(sd.perc),"-tile)\n\n")
     
    } else{
      sd.keep <- rownames(gct.comb$rdesc)
    }
  }
  
  ## #################
  ## filter each dataset separately
  if(mode == 'separate' ){
    cat("\napplying SD-filter:", mode, '\n')
      
    ## filter each type separetely
    sd.keep <- base::tapply(sd.expr, data.type, sd_filt_func, sd_filt) %>% unlist
  }
  
  ## ######################
  ## equal weights
  if(mode == 'equal' ){
    
    cat("\napplying SD-filter:", mode, '\n')
    
    rdesc <- gct.comb@rdesc
    expr <- gct.comb@mat
    
    ## apply global filter first
    sd.keep <- sd_filt_func(sd.expr, sd_filt)
    ## update
    rdesc <- rdesc[sd.keep, ]
    sd.expr <- sd.expr[sd.keep]
    data.type <- data.type[sd.keep]
    
    ## data types
    #data.type <- sub('^(.*?)-.*','\\1', rdesc$Data.Type.ID)
    
    ## determine smallest dataset
    data.type.dist <- table(data.type)
    
    min.numb <- min(data.type.dist <- table(data.type))
    
    ## filter each type separately
    sd.keep <- base::tapply(sd.expr, data.type, function(x, min.numb){
        keep.idx <- order(x, decreasing = T)    
        keep.idx <- keep.idx[1:min.numb]
        names(x)[keep.idx]
    }, min.numb) %>% unlist
    
   
  }
  ##########################
  ## update 
  
  ## GCT
  gct.filt <- gct.comb
  gct.filt@mat <- data.matrix(gct.filt@mat [sd.keep,])
  gct.filt@rdesc <- gct.filt@rdesc[sd.keep, ]
  gct.filt@rid <- rownames(gct.filt@mat)
  
  ## sd vector
  sd.expr.filt <- sd.expr[sd.keep]
  ## data type vector
  data.type.filt <- data.type[sd.keep]
  
  #####################################################
  ## plot number of feature before/after filtering
  data.type.dist <- table(sub('^(.*?)-.*','\\1', gct.comb@rdesc$Data.Type.ID))
  data.type.dist.filt <- table(sub('^(.*?)-.*','\\1', gct.filt@rdesc$Data.Type.ID))
  
  data.type.dist.comb <- rbind(data.type.dist, data.type.dist.filt)
  rownames( data.type.dist.comb) <- c('all', 'filtered')
  col <- c('grey30', 'grey80')
  
  pdf(paste0('0_barplot_number_features_filt-',mode,'.pdf'), width=max(4, length(data.type.dist)*2), height = 6)
  fancyBarplot( data.type.dist.comb, col=col, srt = 45, ylab='No. features')
  legend('top', legend=rownames(data.type.dist.comb), fill=col, title = paste0('mode: ', mode) )
  dev.off()
  
  #####################################################
  ## plot SD of features before/after filtering
  sd.list <- base::tapply(sd.expr, data.type, function(x)x)
  sd.list.filt <- base::tapply(sd.expr.filt, data.type.filt, function(x)x)
  names(sd.list.filt) <- paste0(names(sd.list), '.filt')
  
  sd.list  <- c(sd.list, sd.list.filt)
  sd.list <- sd.list[order(names(sd.list))]
  
  pdf(paste0('0_boxplot_std-dev_filt-',mode,'.pdf'), width=max(4, length(data.type.dist)*2), height = 6)
  p<- try(fancyBoxplot(sd.list, ylab='Std Dev', show.numb = 'median'))
  dev.off()
  
  return(gct.filt)
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


###################################################################
## heatmap using ComplexHeatmap package
MyComplexHeatmap <- function(m, cdesc, cdesc.color, class.variable, variable.other, max.val=NULL, 
                             symm.col=T, ## symmetric color scale centered at zero
                             ##row_title='', 
                             name='NMF features'){
  library(pacman)  
  p_load(circlize)
  p_load(ComplexHeatmap)
  p_load(RColorBrewer)
    
    ## cap values
    if(is.null(max.val)){
      if(symm.col)
        m.max <- ceiling( max(abs(m), na.rm=T) )
      else
        m.max <- ceiling( max(m, na.rm=T) )
    } else {
      m.max <- max.val
      m[m > m.max] <- m.max
      m[m < -m.max] <- -m.max
    }
    ## #####################################
    ## complexheatmap
  
    if(symm.col){
      col.breaks <- seq(-m.max, m.max, length.out=11)
    } else {
      m.min <- floor( min(m, na.rm=T) )
      col.breaks <- seq(m.min, m.max, length.out=11)
    }
    col.hm <- colorRamp2(col.breaks, rev(brewer.pal (11, "RdBu")))
    
    ## column annotation
    cdesc.ha <- HeatmapAnnotation(df=cdesc[ , rev(c(class.variable, variable.other))], col=cdesc.color,
                                  show_legend = T, show_annotation_name = T, annotation_name_side = 'right',
                                  annotation_legend_param=list(
                                    direction='horizontal',
                                    vt_gap = unit(1, 'cm')
                                    
                                    #title_position = "leftcenter"
                                    )
                                  )
    ## heatmap
    hm <- Heatmap(m, col=col.hm,
                  cluster_columns = F,
                  top_annotation = cdesc.ha,
                  
                  row_split = sub('-.*','',rownames(m)),
                  
                  row_title_rot=0,
                  
                  column_split = paste0('C', cdesc$NMF.consensus),
                  column_title_rot=0,
                  
                  row_dend_side = 'right',
                  #name='NMF features',
                  name=name,
                  show_row_names = F,
                  show_column_names = F,
                  use_raster = FALSE)
    ## plot
    draw(hm, annotation_legend_side='bottom')
    
    ## heatmap without annotation
    hm.noanno <- Heatmap(m, col=col.hm,
                  cluster_columns = F,
                  row_split = sub('-.*','',rownames(m)),
                  
                  row_title_rot=0,
                  
                  column_split = paste0('C', cdesc$NMF.consensus),
                  column_title_rot=0,
                  
                  row_dend_side = 'right',
                  name=name,
                  show_row_names = F,
                  show_column_names = F,
                  use_raster = FALSE)##,
                 ## row_title=row_title)
    
    return(list(hm.anno=hm, hm.noanno=hm.noanno))
}

###################################################################
## heatmap using teh pheatmap package
MyPheatMap <- function(m, cdesc, cdesc.color, rdesc=NULL, class.variable, variable.other, filename, color, cw, ch,
                       max.val=NULL, gaps_row=NULL, gaps_col=NULL, cluster_rows=F, cluster_cols=F,  show_rownames=T, ... ){

  if(is.null(max.val)){
    m.max <- ceiling(max(abs(m), na.rm=T))
  } else {
    m.max <- max.val
    m[m > m.max] <- m.max
    m[m < -m.max] <- -m.max
  }

  if(min(m, na.rm=T) < 0)
    breaks <- seq(-m.max, m.max, length.out=12)
  else
    breaks <- NULL

    #breaks <- seq(0, m.max, length.out=12)
  #pheatmap(m, scale='none', annotation_row = rdesc.feat, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color,
  #         cluster_cols = F, cluster_rows=F, filename='6.0_heatmap_ALL_features.pdf', show_rownames=F, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cellwidth = cw)
  ## sort
  pheatmap(m, scale='none', annotation_row = rdesc, annotation_col=cdesc[ , rev(c(class.variable, variable.other))], annotation_colors=cdesc.color,
           cluster_cols = cluster_cols, cluster_rows=cluster_rows, filename=filename, show_rownames=show_rownames,
           breaks=breaks,
           color=color,
           cellwidth = cw,
           cellheight = ch,
           gaps_col = cumsum(table(cdesc$NMF.consensus)), gaps_row = gaps_row,
           ...)
}

################################################
## boxplots comapring continous variables
##
boxplotPerCluster <- function(cdesc,   ## clin.anno
                              nmf2col,     ## vector af length k mapping colors to clusters
                              cont_anno,     ## continous variables to include
                              core_membership=0.5,
                              blank.anno = 'N/A'
                              ){
  
  ######################
  # nclust 
  nclust <- max(as.numeric(cdesc$NMF.consensus))
  names(nmf2col) <- glue("{1:nclust}")
  
  ## continuous variable to plot
  keep.idx <- which(cont_anno %in% colnames(cdesc))
  if(length(cont_anno) == 0){
    warning('The specified variables could not be found. Skipping boxplots...\n')
    return()
  }
  cont_anno <- cont_anno[keep.idx]
  
  ## pairwise comparisons
  comps <- list()
  cc <- 1
  for(i in 1:(nclust-1))
    for(j in (i+1):nclust){
      comps[[cc]] <- c(i, j)
      cc <- cc + 1
    }
  
  
  ###############################
  ## loop over variables
  pdf(paste0('7.0_boxplots-min-membership-',core_membership,'.pdf'), 5, 5)
  for(var in cont_anno){
    
    
    cdesc.filt <- cdesc[!grepl(blank.anno, cdesc[, var]), ]
    cdesc.filt[, var] <- as.numeric(cdesc.filt[, var])
    
    cdesc.filt <- cdesc.filt %>% filter(NMF.cluster.membership > core_membership) 
    #cdesc.filt <- 
      
    #p <-# cdesc %>% filter(NMF.cluster.membership > core_membership) %>%
        #filter(!grepl('N/A', `var`)) %>%
      
    p <-  ggboxplot(cdesc.filt, x="NMF.consensus", y=var, add='jitter', palette=nmf2col, color = "NMF.consensus", size=1.5) + 
        ggtitle(var) +
        stat_compare_means(comparisons=comps)
    
    
    plot(p)
    
    # p <- try(
    #   cdesc %>% filter(NMF.cluster.membership > core_membership) %>%
    #            filter(!grepl('N/A', var)) %>%
    #       ggboxplot(., x="NMF.consensus", y=var, add='jitter', palette=nmf2col, color = "NMF.consensus", size=1.5) + 
    #       ggtitle(var) +
    #       stat_compare_means(comparisons=comps)
    #   )
    # 
    # if(class(p) != 'try-error') plot(p)
    # 
  }
  dev.off()

}

## #############################################
## main function ...
##
nmf.post.processing <- function(ws,                       ## filename of R-workspace
                                blank.anno = 'N/A',       ## used to replace blanks/NAs in meta data
                                blank.anno.col = 'white', ## color for blanks/NAs usind in heatmap annotation tracks
                                #blank.anno.col = 'grey90', ## color for blanks/NAs usind in heatmap annotation tracks
                                
                                core_membership=0.5,      ## NMF.cluster.membership score to define core memebrship
                                                          ## cluster enrichment of clinical variables will be done on the core set 
                                feature.fdr=0.01,         ## FDR for NMF features after 2-sample mod T (cluster vs. rest)
                                pval.ora= 0.01,           ## p-value for overrepresenation analysis of core cluster with categorial metadata 
                                max.categories=10,        ## max. number of levels in categorial metadata to be included on the 
                                                          ## overrepresentation analysis
                                organism=c('human', 'mouse', 'rat') ## required to annotate features
                      ){

    ## ############################
    ## import workspace
    cat('loading ws...')
    load(ws)
    cat('done.\n')

    ###############################
    ## organism
    organism <- match.arg(organism)
    if(organism == 'human') org.id <- 'Hs'
    if(organism == 'mouse') org.id <- 'Mm'
    if(organism == 'rat') org.id <- 'Rn'
    
    ## heatmap aesthetics
    cw <- opt$hm_cw
    ch <- opt$hm_ch
    if(opt$z_score){
      max.val <- opt$hm_max_val_z
    } else {
      max.val <- opt$hm_max_val
    }
    
    ###################################  
    ## get data type of each column
    cdesc <- as.tibble(cdesc)
    cc <- cdesc %>%
      dplyr::summarise_all(class) %>%
      tidyr::gather(variable, class)

    ## fix blanks in cdesc
   cdesc <- apply(cdesc, 2, function(x){
      x[nchar(x) == 0 | is.na(x)] <- blank.anno
      x
    })
   ## #######################################
   ## harmonize column types in cdesc
   cdesc <- as.tibble(cdesc)
   int.idx <- which(cc$class == 'integer')
   if(length(int.idx) > 0){
    for(i in int.idx)
      cdesc %<>% mutate_at(i, as.integer)
   }
   num.idx <- which(cc$class == 'numeric')
   if(length(num.idx) > 0){
   for(i in num.idx)
     cdesc %<>% mutate_at(i, as.numeric)
   }
   char.idx <- which(cc$class == 'character')
   if(length(char.idx) > 0){
   for(i in num.idx)
     cdesc %<>% mutate_at(i, as.numeric)
  }

  ## ###################################################
  ## class variable/annotation tracks of interest
  cat_anno <- strsplit( opt$cat_anno, ';' ) %>% unlist
  class.variable <- cat_anno[1]
  class.variable <- gsub('\\\'', '', class.variable)
  
  ## other categorial variables of interest
  if(length(cat_anno) > 1){
    variable.other <- cat_anno[2:length(cat_anno)]
    variable.other <- gsub('\\\'', '', variable.other)
    keep.idx <- variable.other %in% colnames(cdesc)
    variable.other <- variable.other[ keep.idx ]
  } else {
    variable.other <- c()
  }
  ## continous variables
  if(!is.na(opt$cont_anno)){
    cont_anno <- strsplit( opt$cont_anno, ';' ) %>% unlist
    cont_anno <- gsub('\\\'', '', cont_anno)
  } else{
    cont_anno <- NA
  }
  

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

  ## ###############################
  ## parse colors
  cdesc.color <- parse.colors(opt, data.frame(cdesc.org), blank = blank.anno, blank.col = blank.anno.col)
  ## keep a copy
  cdesc.color.org <- cdesc.color

  ## ################################
  ## parse genes of interest
  gene_col <- opt$gene_col

  ## ####################################################
  ##              loop over ranks
  for(rank in rank.top) {

    dir.create(paste('K', rank, sep='_'))
    setwd(paste('K', rank, sep='_'))

    ##################################
    ## extract NMF results
    res <- res.rank[[as.character( rank )]]

    #################################
    ## colors
    cdesc.color <- cdesc.color.org
    cdesc <- cdesc.org
    variable.other <- variable.other.org
    class.variable <- class.variable.org

    
    # ## ##########################################
    # ## silhoutte plots
    # if(!opt$bnmf){
    #   pdf(paste('1.0_silhouette_K_', rank, '.pdf', sep=''), 10, 6)
    #   #par(mfrow=c(1,2))
    #   plot(rank.sil[[rank]], main=paste('K=', rank, sep=''), col=palette()[1:as.numeric(rank)+1] )
    #   #plot(rank.sil.random[[rank]], main=paste('Randomized data', sep=''), col=palette()[1:as.numeric(rank)+1])
    #   dev.off()
    # }
    # 
    # 
    ## data matrix
    expr <- expr.org

    #############################################
    ## coefficients matrix H
    if(opt$bnmf){
      H <- res@coeff[[1]]
    } else {
      H <- res@fit@H
    }

    ## ########################################
    ## determine cluster assignments and
    ## add NMF classification to cdesc object
    if(opt$bnmf){
      NMF.basis <- apply(H, 2, which.max)
      NMF.consensus <- NMF.basis
    } else {
      NMF.basis <- predict(res)
      NMF.consensus <- predict(res, 'consensus')
    }

    #########################################
    ## cluster membership score:
    ## - fractional
    NMF.cluster.membership.alldigits <- apply(H, 2, function(x) max(x/sum(x)))
    NMF.cluster.membership <- round( NMF.cluster.membership.alldigits, 3)
    
    
    ## map consensus to basis
    NMF.consensus.map <- NMF.consensus
    for(ii in 1:as.numeric(rank))
      NMF.consensus[ NMF.consensus == ii ] <-  NMF.basis[NMF.consensus == ii]
    
    ##########################################
    ## define core membership
    NMF.consensus.core <- NMF.consensus
    NMF.consensus.core[ which(NMF.cluster.membership.alldigits < core_membership) ] <- NA
    
    ## add to cdesc
    cdesc <- data.frame(cdesc, 
                        NMF.consensus=as.factor(NMF.consensus), 
                        NMF.consensus.core=as.factor(NMF.consensus.core), 
                        NMF.cluster.membership=NMF.cluster.membership)

    ## export
    write.table(cdesc, file=paste('clin_anno_nmf.txt'), sep='\t', quote=F, na='', col.names=NA)

    ## ################################################
    ## characterize clusters by calculating enrichments
    ## of levels in 'cdesc'
    ## - enrichment is calculated on the core set of samples
    variables.all <- unique(c(class.variable, variable.other))
    variables.all <- variables.all[variables.all %in% colnames(cdesc)]
    
    ## fewer than 10 levels
    enrich.idx <- which(sapply( variables.all, function(x) ifelse(length(unique(cdesc[, x])) <= max.categories, 1, 0)) == 1)
    
    ## exclude 'N/A' category
    
    
    clust.class.enrichment <- lapply(variables.all[enrich.idx], function(o){
      core.idx  <- which(!is.na(cdesc$NMF.consensus.core))
      
      tmp=CalcClustEnrich(cdesc$NMF.consensus[core.idx], cdesc[core.idx, o])
      colnames(tmp) <- paste(o,colnames(tmp), sep=':')
      tmp
    })
    names(clust.class.enrichment) <- variables.all[enrich.idx]
    
    ## make matrix
    cons.map.mat <- Reduce(cbind, clust.class.enrichment)
    cons.map.mat.signif <- apply(cons.map.mat, 1, function(x) paste( names(x)[x < pval.ora], collapse='|') )
    cons.map.mat.signif <- data.frame(cluster=names(cons.map.mat.signif), enriched=cons.map.mat.signif)

    write.table(t(cons.map.mat), sep='\t', file=paste('cluster-enrichment.txt', sep=''), quote=F, row.names=T, col.names = NA)
    write.table(cons.map.mat.signif, sep='\t', file=paste('cluster-enrichment-signif-nom-p-', pval.ora, '.txt', sep=''), quote=F, row.names=F)

    ###########################################################
    ## significantly enriched in 'class variable'
    cons.map.max <- clust.class.enrichment[[class.variable]]
    colnames(cons.map.max) <- sub('.*\\:', '', colnames(cons.map.max))
 
    ## export
    write.table(cons.map.max, col.names = NA, sep='\t', file=paste('nmf_vs_', class.variable, '.txt', sep=''), quote=F)
    cons.map.max <- apply(cons.map.max, 1, function(x) paste(sub('.*\\:','', names(x))[x < pval.ora], collapse='|'))
    not.mapped.idx <- which(nchar(cons.map.max) == 0 )
    
    if(length(not.mapped.idx) > 0)
      cons.map.max[not.mapped.idx] <- as.character(not.mapped.idx)

    cons.map.max <- strsplit(cons.map.max, '\\|') %>% lapply(., function(x){
      if(length(x) > 1 & blank.anno %in% x) {
          keep.idx <- which(x != blank.anno)
          return(x[keep.idx])
      } else {
        return(x)
      }
      
      }) %>% unlist
    
    ####################################################################
    ## map colors of 'class.variable' to NMF clusters
    mapped.idx <- which(cons.map.max %in% names(cdesc.color[[class.variable]]))
    NMF.consensus.col <- rep(NA, length(cons.map.max))
    NMF.consensus.col[mapped.idx] <- cdesc.color[[class.variable]][ cons.map.max[mapped.idx] ]
    names(NMF.consensus.col) <- cons.map.max

    ## add NMF to heatmap annotation tracks
    variable.other <- c('NMF.consensus',  'NMF.consensus.mapped', 'NMF.consensus.core', 'NMF.cluster.membership', variable.other)
  
    ######################################################
    ## check whether all NMF cluster have color
    if( sum(is.na(NMF.consensus.col)) > 0 ){
      idx.tmp <- which( is.na(NMF.consensus.col) )
      col.tmp <-  rev( alphabet() )[ (1:length(idx.tmp) ) +  1]
      NMF.consensus.col[idx.tmp] <- col.tmp
    }
    cdesc$NMF.consensus.mapped <-  cdesc$NMF.consensus
    cdesc.color$NMF.consensus.mapped <- NMF.consensus.col
    names(cdesc.color$NMF.consensus.mapped) <- 1:length(NMF.consensus.col)
    
    ##
    cdesc.color$NMF.consensus.core <- c(NMF.consensus.col,  blank.anno.col)
    names(cdesc.color$NMF.consensus.core) <- c(1:length(NMF.consensus.col), blank.anno)
    
    
    ## ##########################################
    ## silhoutte plots
    if(!opt$bnmf){
       pdf(paste('1.0_silhouette_K_', rank, '.pdf', sep=''), 10, 6)
       plot(rank.sil[[rank]], main=paste('K=', rank, sep=''), col=NMF.consensus.col)#col=palette()[1:as.numeric(rank)+1] )
       dev.off()
    }
    
    #######################################################
    ##
    ##     boxplots comparing continuous variables
    ##
    if(!is.na(cont_anno)){
         try( boxplotPerCluster(cdesc=cdesc,   ## clin.anno
                            nmf2col=cdesc.color$NMF.consensus.mapped,     ## vector af length k mapping colors to clusters
                            cont_anno=cont_anno,     ## continuous variables to include
                            core_membership=core_membership)
              )
    }
    

    ## ########################################
    ## heatmap coefficients

    ## order according to nmf clustering and
    ## cluster membership
    nmf.ord <-  rownames(cdesc)[order( -1*as.numeric(cdesc$NMF.consensus), cdesc$NMF.cluster.membership, decreasing = T)]

    expr <- expr[, nmf.ord]
    cdesc <- cdesc[nmf.ord, ]
    H <- H[, nmf.ord]
    H.norm <- apply(H, 2, function(x)x/max(x))

    #############################
    ## export H and H.norm
    gct.H <- new('GCT')
    gct.H@mat <- H
    gct.H@cdesc <- cdesc
    gct.H@rid <- paste0('C',1:nrow(H))
    gct.H@cid <- rownames(cdesc)
    write.gct(gct.H, ofile=as.character(glue('pattern-matrix-H')))
    gct.H@mat <- H.norm
    write.gct(gct.H, ofile=as.character(glue('pattern-matrix-H-norm-by-max')))
    
        
    ## cluster
    ## pdf
    MyPheatMap(H.norm, cdesc=cdesc, cdesc.color=cdesc.color, rdesc=NULL, class.variable=class.variable,
               variable.other=variable.other, filename='2.2_coefmap_pheatmap_clustered.pdf', cluster_cols=T,
               color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), cw=cw, ch=ch, max.val=NULL)
    ## png
    MyPheatMap(H.norm, cdesc=cdesc, cdesc.color=cdesc.color, rdesc=NULL, class.variable=class.variable,
               variable.other=variable.other, filename='2.2_coefmap_pheatmap_clustered.png', cluster_cols=T,
               color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), cw=cw, ch=ch, max.val=NULL)
    
    ## disable clustering
    MyPheatMap(H.norm, cdesc=cdesc, cdesc.color=cdesc.color, rdesc=NULL, class.variable=class.variable,
               variable.other=variable.other, filename='2.2_coefmap_pheatmap_sorted.pdf',
               color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), cw=cw, ch=ch, max.val=NULL)
    ## png
    MyPheatMap(H.norm, cdesc=cdesc, cdesc.color=cdesc.color, rdesc=NULL, class.variable=class.variable,
               variable.other=variable.other, filename='2.2_coefmap_pheatmap_sorted.png',
               color=colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100), cw=cw, ch=ch, max.val=NULL)
    
    ## #######################################
    ## consensus map
    if(!opt$bnmf){
      cons <- res@consensus[nmf.ord, nmf.ord]

      res2 <- res
      res2@consensus <- cons
      consensusmap(res2, filename=paste('3.0_consensusmap_nrun_', opt$nrun, '.pdf', sep=''))

      ## pheatmap version
      MyPheatMap(cons, cdesc=cdesc, cdesc.color=cdesc.color, rdesc=NULL, class.variable=class.variable,
                 variable.other=variable.other, filename=paste('3.1_consensusmap_nrun_', opt$nrun,'_pheatmap.pdf', sep=''),
                 color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), cw=cw, ch=ch, max.val=NULL)
      ## png
      MyPheatMap(cons, cdesc=cdesc, cdesc.color=cdesc.color, rdesc=NULL, class.variable=class.variable,
                 variable.other=variable.other, filename=paste('3.1_consensusmap_nrun_', opt$nrun,'_pheatmap.png', sep=''),
                 color=colorRampPalette(c('blue', 'blue4','darkred', 'red'))(100), cw=cw, ch=ch, max.val=NULL)
      
      }

    ## ##############################################################
    ##
    ##           extract meta feature matrix
    ##
    ## used as ranking for cluster marker selection
    ## ##############################################################
    if(opt$bnmf){
      W <- res@basis[[1]]
      colnames(W) <- 1:ncol(W)
    } else {
      W <- res@fit@W
      colnames(W) <- 1:ncol(W)
    }
    
    ## row-normalize
    W.norm <- t(apply(W, 1,  function(x)x/sum(x)))

    ## combine up/down: take highest normalized score
    feat.comb <- unique( sub('^(.*-.*)_(up|down)$' , '\\1', rownames(W.norm)) )

    W.norm.dir <- W.norm
    direction.tmp <- sub('.*_(.*)$', '\\1', rownames(W.norm.dir))

    ## sign the coefficients
    direction.tmp[direction.tmp == 'down' ] <- -1
    direction.tmp[direction.tmp == 'up' ] <- 1
    direction.tmp <- as.numeric(direction.tmp)
    W.norm.dir <- direction.tmp*W.norm.dir

    ## ########################################
    ## combine up/down
    cores <- detectCores() - 1
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    W.tmp2 <- foreach(i = 1:length(feat.comb)) %dopar% {
        library(glue)
        ii=feat.comb[i]
        idx.tmp <- grep(glue("^{ii}_(up|down)$"), rownames(W.norm))
        W.tmp <- W.norm.dir[ idx.tmp, ]
        if(length(idx.tmp) > 1){
          W.tmp <- apply(W.tmp, 2, function(x) x[ which.max(abs(x)) ])
        }
        ## if grep didn't find the current id
        ## try to remove special characters etc.
        if(length(idx.tmp) == 0){
          idx.tmp <- grep(glue("^{ii}_(up|down)$"), rownames(W.norm))
          W.tmp <- W.norm.dir[ idx.tmp, ]
        }
        W.tmp
    }
    on.exit(stopCluster(cl))

    ## matrix
    W.norm.comb <- Reduce(rbind, W.tmp2)
    dimnames(W.norm.comb) <- list(feat.comb, colnames(W))

    ## create GCT file
    if(opt$gene_col %in% colnames(rdesc)){
      w.comb.rdesc <- data.frame(
        Type=sub('^(.*?)-.*','\\1', rownames(W.norm.comb)),
        geneSymbol=gct.comb@rdesc[rownames(W.norm.comb), opt$gene_col],
        stringsAsFactors = F)

      ## if 'opt$gene_col' was not present in all tables
      ## use rownames to extract gene symbols...
      if(sum(is.na(w.comb.rdesc$geneSymbol)) > 0){
        na.idx <- which(is.na(w.comb.rdesc$geneSymbol))
        w.comb.rdesc$geneSymbol[na.idx] <- sub('^(.*?)-(.*)$','\\2', rownames(W.norm.comb)[ na.idx ])
      }

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
    write.gct(w.gct, ofile="matrix_W_combined_signed.gct", appenddim = F)

    ################################################################
    ##
    ##                      extract features
    ##
    ## ##############################################################
    feature.methods <- list('kim', 'max', 1, 10L)
    feature.methods.scores <- list('kim', 'max', 'max', 'max')

    success <- F
    pdf('4.0_histogram_feature_scores.pdf')
    j <- 1
    while(!success){

      feature.method=feature.methods[[j]]
      feature.method.score=feature.methods.scores[[j]]

      ## feature scores
      s.i <- featureScore(W, method=feature.method.score)
      hist( unlist(s.i), main=paste("Feature scores (method:", feature.method,")"), xlab='Feature score' )

      ## extract features
      s <- extractFeatures(W, method=feature.method)

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
    names(s) <- paste(1:length(s))

    #################################
    ## cluster with no features
    rm.idx <- which( sapply(s, function(x) sum(is.na(x))) > 0)
    if(length(rm.idx) > 0){
      s <- s[-c(rm.idx)]
      cons.map.mat <- cons.map.mat[-c(rm.idx)]
      cons.map.max <- cons.map.max[-c(rm.idx)]
    }

    
    ## ###########################################################
    ##                annotate the features
    ## ###########################################################
    
    ## get accession numbers
    s.acc <- lapply(s, function(x) data.frame(Accession=sub('_up|_down','',rownames(expr.comb)[x]) , id=rownames(expr.comb)[x], stringsAsFactors = F) )
    
    ## features plus scores
    s.acc.scores <- lapply(s.acc, function(x) data.frame(x, Score=s.i[x$id], stringsAsFactors = F))

    ## ####################################
    ## gene names
    if(opt$gene_col %in% colnames(rdesc)){
      s.gn <- lapply(s.acc, function(x) rdesc[x$Accession, opt$gene_col] )
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
        NMF.Score=s.acc.scores[[i]]$Score,
        stringsAsFactors=F)
      
      ## replace NAs in SYMBOL by Accession 
      na.idx <- which(is.na(tmp$SYMBOL))
      if(length(na.idx) >0 )
        tmp$SYMBOL[na.idx] <- tmp$Accession[na.idx] 
      
      
      s.acc.gn[[i]] <- tmp[order(tmp$SYMBOL),]
      
    }
    s.acc.gn.anno <- s.acc.gn

    ## ################################################################
    if(opt$gene_col %in% colnames(rdesc)){

      ## add description and enzyme codes
      s.acc.gn.anno <- lapply( s.acc.gn, function(x)
        AnnotationDbi::select(eval(parse(text=paste0("org.",org.id,".eg.db"))), keys=x$SYMBOL , column=c( 'GENENAME',  'ENZYME'), keytype='SYMBOL', multiVals='first')
      )

      ## add cytoband
      map <- org.Hs.egMAP
      map.genes <- mappedkeys(map)
      map <- as.list(map[map.genes])
      
      s.acc.gn.anno.cyto <- lapply( s.acc.gn.anno, function(x){
        
        entrez.id <- mapIds(eval(parse(text=paste0("org.",org.id,".eg.db"))), keys = x$SYMBOL, keytype = 'SYMBOL', column = 'ENTREZID')
        
        ## check for NULL
        null.idx <- sapply(entrez.id, is.null)
        if(sum(null.idx) > 0){
          names.tmp <- names(entrez.id)
          entrez.id[null.idx] <- 'NA'
          names.tmp[null.idx] <- 'NA'
          entrez.id <- make.unique(unlist(entrez.id))
          names(entrez.id) <- make.unique(names.tmp)
        }
        
        #if(class(entrez.id) == 'try-error') return(data.frame(SYMBOL=x$SYMBOL, ENTREZID='', CYTOBAND=''))
        
        cytoband <- map[entrez.id] 
        cytoband[which(sapply(cytoband, is.null))] <- NA
        cytoband <- lapply(cytoband, paste, collapse='|')
        cytoband <- unlist(cytoband)
        #cytoband
        data.frame(SYMBOL=x$SYMBOL, ENTREZID=entrez.id, CYTOBAND=cytoband)
        
      })
      
      ## join in a single data frame
      for(i in 1:length(s.acc)){
        #tmp <- merge(s.acc.gn[[i]], s.acc.gn.anno.cyto[[i]], 'SYMBOL' )
        tmp <- left_join(s.acc.gn[[i]], s.acc.gn.anno.cyto[[i]], 'SYMBOL' )
        #tmp <- merge(tmp,  s.acc.gn.anno[[i]], 'SYMBOL')
        tmp <- left_join(tmp,  s.acc.gn.anno[[i]], 'SYMBOL')
        
        #############################
        ## mapping to annotation dataabses might result
        ## in 1:many mappings and will introduce redunancy
        ## remove identical rows
        tmp <- tmp[which(!duplicated(apply(tmp, 1, paste, collapse=' '))), ]
        
        dup.idx <- which(duplicated(tmp$Accession))
        if(length(dup.idx) > 0){
          tmp <-  aggregate(tmp, list(tmp$Accession), function(x) paste(unique(x), collapse='|'))
          tmp <- tmp[, -c(1)]
        }
        
        s.acc.gn.anno[[i]] <- tmp
      }

    } ## end if(opt$gene_col %in% colnames(rdesc))

    ## sort
    s.acc.gn.anno <- lapply(s.acc.gn.anno, function(x)x[order(x$NMF.Score, decreasing=T), ])

    #############################################################
    ##
    ##   for each cluster run a one-vs-rest 2-sample mod T
    ##
    #############################################################
    s.acc.gn.anno.modT <- lapply(names(s.acc.gn.anno), function(cl){
      
      feat.cl <- s.acc.gn.anno[[cl]]
      feat.ids <- paste(feat.cl$Type, feat.cl$Accession, sep='-')
      
      feat.ids <- feat.ids[feat.ids %in% rownames(expr)]
      
      m <- data.frame(Accession=feat.ids, expr[ feat.ids, ])
      m$Accession <- sub('^.*?-', '', m$Accession)
      
      ## class vector
      class.vec <- paste0("C",as.character(cdesc$NMF.consensus))
      names(class.vec) <- rownames(cdesc)
      
      ## direction of test is determined by sorting the class labels
      ## alphanumerically -> 'CX over allothers' 
      class.vec[ class.vec != paste0('C', cl) ] <- 'allother'
      
      ## test
      res.tmp <-  try(modT.test.2class( m, groups=class.vec, id.col='Accession', label=glue("C{cl}.over.rest") )$output)
      if(class(res.tmp) == 'try-error')
        return(feat.cl)
    
      colnames(res.tmp)[which(colnames(res.tmp) == 'id')] <- 'Accession'
      
      feat.cl <- full_join(res.tmp[, c(1, grep(glue("C{cl}.over.rest$"), colnames(res.tmp)))], feat.cl, 'Accession')
     
      feat.cl[ order( feat.cl[, glue("P.Value.C{cl}.over.rest")] ), ]
    })
    names(s.acc.gn.anno.modT) <- names(s.acc.gn.anno)
    
    #############################################################
    ## re-calculate FDR
    #############################################################
    all.pval <- lapply( s.acc.gn.anno.modT , function(x) x[, grep('^P\\.Value\\..*', colnames(x)) ])
    all.pval.fdr <- p.adjust(unlist(all.pval), method = 'BH')
    all.pval.fdr2 <- all.pval
    
    ## split the cluster
    start <- 1
    for(i in 1:length(all.pval)){
      all.pval.fdr2[[i]] <- all.pval.fdr[start:(start+length(all.pval[[i]])-1)]
      start <- start + length(all.pval[[i]])
    }
    ## replace in table
    for(i in 1:length( all.pval.fdr2))
      s.acc.gn.anno.modT[[i]][, grep('^adj\\.P\\.Val\\..*', colnames(s.acc.gn.anno.modT[[i]]))] <- all.pval.fdr2[[i]]
      
    ############################################
    ## filter
    s.acc.gn.anno.modT <- lapply( s.acc.gn.anno.modT, function(x)x[which(x[, grep('^adj\\.P\\.Val\\..*', colnames(x))] < feature.fdr  ), ])
    keep.idx <- which(sapply(s.acc.gn.anno.modT, nrow) > 0)
    s.acc.gn.anno.modT <- s.acc.gn.anno.modT[ keep.idx ]
    
    ## ##########################################################
    ## export as Excel file
    ## all features
    list2env(s.acc.gn.anno.modT , envir=.GlobalEnv)
    all.features <- lapply(s.acc.gn.anno.modT, function(x)x$Accession) %>% unlist %>% unique
  
    nmf.marker.str <- paste('NMF_features_N_', length(all.features),'.xlsx', sep='')
    ##nmf.marker.str <- paste('NMF_features_N_', length(unique(unlist(s.acc2))),'.xlsx', sep='')
    WriteXLS(names(s.acc.gn.anno.modT), ExcelFileName=nmf.marker.str, FreezeRow=1, FreezeCol=1, SheetNames=make.unique(substr(names(s.acc.gn.anno.modT), 1, 20)), row.names=F, BoldHeaderRow=T, AutoFilter=T)
   
    
    ## ##########################################################
    ## extract protein kinases and phosphatases
    s.acc.gn.kinase <- lapply(s.acc.gn.anno.modT, function(x) x[ grep('2\\.7\\.1[0-2]|3\\.1\\.3\\.(16|48)', x$ENZYME), ])
    names(s.acc.gn.kinase) <- paste(names(s.acc.gn.kinase), 'KINASE', sep='_')

    s.acc.gn.kinase <- lapply(s.acc.gn.kinase, function(x) x[which(!duplicated(apply(x, 1, paste, collapse=' '))),] )

    ## kinases
    list2env(s.acc.gn.kinase , envir=.GlobalEnv)
    all.kinase.features <- lapply(s.acc.gn.kinase, function(x)x$Accession) %>% unlist %>% unique
    
    WriteXLS(names(s.acc.gn.kinase), ExcelFileName=paste('NMF_KINASE_PHOSPHATASE_features_', length(all.kinase.features),'.xlsx', sep=''), FreezeRow=1, FreezeCol=1, SheetNames=make.unique(substr(names(s.acc.gn.kinase), 1, 20)), row.names=F, BoldHeaderRow=T, AutoFilter=T)

    
    
    ## #################################################
    ## UpSet plot of features per cluster
    ##
    if(length(s.acc.gn.anno.modT) > 1) {
              upset.mat <- matrix(0, ncol=length(s.acc.gn.anno.modT), nrow=length(unique(unlist(sapply(s.acc.gn.anno.modT, function(x)x$Accession )))), 
                                  dimnames=list(unique(unlist(sapply(s.acc.gn.anno.modT, function(x)x$Accession ))), paste0('C',names(s.acc.gn.anno.modT))))
              for(ii in names(s.acc.gn.anno.modT))
                upset.mat[s.acc.gn.anno.modT[[ii]]$Accession, paste0('C',ii)] <- 1
              
              pdf('4.1_upset_NMF_markers.pdf')
              upset(data.frame(upset.mat), point.size = 4, text.scale = 1.5, nintersects = NA, nsets = length(s.acc.gn.anno.modT))
              dev.off()
              png('4.1_upset_NMF_markers.png')
              upset(data.frame(upset.mat), point.size = 4, text.scale = 1.5, nintersects = NA, nsets = length(s.acc.gn.anno.modT))
              dev.off()
    }


    ## #####################################################################################
    ##
    ##               contribution of each data type to NMF coefficients
    ##
    ## ######################################################################################
    ##data.types <- lapply(s.acc.gn.anno.modT, function(x) table(sub('-.*','', x$Accession)))
    data.types <- lapply(s.acc.gn.anno.modT, function(x) table(x$Type))
    
    tmp <- lapply(data.types, names)
    tmp <- unique(unlist(tmp))
    data.types.mat <- matrix(0, nrow=length(data.types), ncol=length(tmp), dimnames=list(names(s.acc.gn.anno.modT), tmp) )
    for(i in rownames(data.types.mat))
      data.types.mat[i, names(data.types[[i]])] <- data.types[[i]]

    ## add mapping to class variable
    idx.tmp <- match( names(cons.map.max[names(s.acc.gn.anno.modT)]), rownames(data.types.mat))
    ## commented out on 2020/09/04 
    #rownames(data.types.mat)[idx.tmp] <- paste(rownames(data.types.mat)[idx.tmp], ' (', cons.map.max, ')', sep='' )

    col <- RColorBrewer::brewer.pal(length(tmp), "Set2")[1:length(tmp)]


    ## #############
    ## barplot
    names(col) <- tmp
    pdf('4.2_barplot_features_per_basis2.pdf')
    fancyBarplot(t(data.types.mat), srt=30, xlab='NMF basis', col=col, main='Proteogenomic features', ylab='No. features')
    legend('top', legend=names(col), fill=col, bty='n')
    dev.off()
    png('4.2_barplot_features_per_basis2.png')
    fancyBarplot(t(data.types.mat), srt=30, xlab='NMF basis', col=col, main='Proteogenomic features', ylab='No. features')
    legend('top', legend=names(col), fill=col, bty='n')
    dev.off()
    
    
    ## ##############
    ## table
    m <- data.types.mat
    ## pdf
    try(pheatmap(m, cluster_rows=F, cluster_cols=F, display_numbers=T, legend=F, number_format='%.0f', cellwidth=30, cellheight=30, fontsize=15, number_color='grey20', main='NMF features', scale='none', col='white', border_color='blue', filename='4.3_matrix_features_per_NMF.consensus.pdf'))
    ## png
    try(pheatmap(m, cluster_rows=F, cluster_cols=F, display_numbers=T, legend=F, number_format='%.0f', cellwidth=30, cellheight=30, fontsize=15, number_color='grey20', main='NMF features', scale='none', col='white', border_color='blue', filename='4.3_matrix_features_per_NMF.consensus.png'))
    

    ## ###########################################################
    ## plot the actual expression values for the extracted features
    cc <- 1
    hm.list <- list()
    for(i in names(s.gn[ names(s.acc.gn.anno.modT) ] )){

      acc <- paste(s.acc.gn.anno.modT[[i]]$Type, s.acc.gn.anno.modT[[i]]$Accession, sep='-')
      m <- expr[acc , ]

      if( !is.null(nrow(m))){

        ## plot pdf
        pdf(paste('5.',cc-1,'_ComplexHeatmap_expression_C', i, '.pdf', sep=''), 13, 12)
        hm.list[[cc]] <- MyComplexHeatmap(m, cdesc, cdesc.color, class.variable, variable.other, max.val) #, row_title=paste0('C',i))
        dev.off()

        ## plot png
        png(paste('5.',cc-1,'_ComplexHeatmap_expression_C', i, '.png', sep=''),  width=13, height=12, units = 'in', res=100)
        draw(hm.list[[cc]]$hm.anno)
        dev.off()
        
        ###############################################
        ## assemble GCT
        gct <- new('GCT')
        gct@mat <- data.matrix(m)
        gct@rdesc <- data.frame(s.acc.gn.anno.modT[[i]])
        gct@cdesc <- data.frame(cdesc)
        gct@cid <- colnames(m)
        gct@rid <- rownames(m)
        write.gct(gct, ofile = as.character(glue("5.{cc-1}_data_matrix_C{i}")))
        
        names(hm.list)[cc] <- i
        cc <- cc + 1
        
      } # end !is.null(nrow(m))

    }
    #############################
    ## concatenate heatmaps
    if(length(hm.list) > 1){
      hm.concat <- hm.list[[1]]$hm.anno
      
      for(i in 2:length(hm.list))
        hm.concat <- hm.concat %v% hm.list[[i]]$hm.noanno
      
      ## draw heatmap
      pdf('6.0_ComplexHeatmap_ALL_features-concat.pdf', 16, 16)
      draw(hm.concat, annotation_legend_side='bottom')
      dev.off()
      
      png('6.0_ComplexHeatmap_ALL_features-concat.png', width=16, height=16, units = 'in', res=100)
      draw(hm.concat, annotation_legend_side='bottom')
      dev.off()
      
    }
    ## ############################################################
    ##
    ##                plot ALL markers
    ##
    feat.all <- unique(unlist(lapply(s.acc.gn.anno.modT, function(x)paste(x$Type, x$Accession, sep='-')))) 
    m <- expr[ feat.all, ]
   
    ## #########################################################
    ##  complex heatmap
    pdf('6.1_ComplexHeatmap_ALL_features.pdf', 12, 16)
    MyComplexHeatmap(m, cdesc, cdesc.color, class.variable, variable.other, max.val)$hm.anno
    dev.off()

    png('6.1_ComplexHeatmap_ALL_features.png', width=12, height=16, units='in', res=100)
    MyComplexHeatmap(m, cdesc, cdesc.color, class.variable, variable.other, max.val)$hm.anno
    dev.off()
    
    ## ############################################################
    ## plot ALL markersbut without CNV
    if(length( grep('^CNV-', rownames(m) )) > 0 ){

      m.org <- m
      cnv.idx <- grep('CNV-', rownames(m.org))

      if(length(cnv.idx) > 0){

        m <- m[-cnv.idx, ]
        ## #####################################
        ## complexheatmap
        pdf('6.2_ComplexHeatmap_ALL_features_no_CNV.pdf', 11, 16)
        MyComplexHeatmap(m, cdesc, cdesc.color, class.variable, variable.other, max.val)$hm.anno
        dev.off()

        png('6.2_ComplexHeatmap_ALL_features_no_CNV.png', width=11, height=16, units='in', res=100)
        MyComplexHeatmap(m, cdesc, cdesc.color, class.variable, variable.other, max.val)$hm.anno
        dev.off()
        
        ## ##############################################################################
        ## CNV only
        if(length(cnv.idx) == 1) {
          m <- data.frame(m.org[cnv.idx, ]) %>% t
        } else {
          m <- m.org[cnv.idx, ]
        }
        
        ## #####################################
        ## complexheatmap
        pdf('6.2_ComplexHeatmap_ALL_features_CNV_only.pdf', 11, 16)
        MyComplexHeatmap(m, cdesc, cdesc.color, class.variable, variable.other, max.val=1)
        dev.off()

      }
    }

    ## #############################################################################################
    ##
    ##                        positive conrols / genes of interest
    ##
    ## #############################################################################################
    s.acc.gn.anno.valid <- lapply(s.acc.gn.anno.modT, function(x) x[paste(x$Type, x$Accession, sep='-') %in% feat.all,])
    genes.of.interest <- try(lapply(s.acc.gn.anno.valid, function(x) x$SYMBOL) %>% unlist %>% unique)
    
    #feat.all 
    #genes.of.interest <- genes.of.interest[1:10]

    if(class(genes.of.interest) != 'try-error'){
    ##if(!is.null(genes.of.interest)){

      dir.create('nmf-features')
      setwd('nmf-features')

      # pdf('barchart_positive_controls.pdf', 3, 3)
      for(i in genes.of.interest){

        ## extract features
        feat <- lapply(s.acc.gn, function(x){
          idx=grep(glue("{i}"), x$SYMBOL, value=F, ignore.case=F)
          paste(x$Type[idx], x$Accession[idx], sep='-' )
        })

        if(sum(sapply(feat, length)) > 0){

          feat <- lapply(feat, function(x) x[x %in% feat.all])
          
          feat.expr <- expr.org[unlist(feat), ]

          ## order: subtype + NMF
          ord.idx <- with(cdesc, order( NMF.consensus))

          #####################################
          ## heatmap

          ## if there is only a single row
          if(is.null(dim(feat.expr))){

            names(feat.expr) <- colnames(expr.org)
            feat.expr <- feat.expr[ rownames(cdesc)[ord.idx] ]
            tmp <- feat.expr

            dim(tmp) <- c(1, length(feat.expr))
            colnames(tmp) <-  names(feat.expr)
            rownames(tmp) <- unlist(feat)
            m <- tmp

            ## heatmap
            m.max <- ceiling(max(abs(m), na.rm=T))
            #m.max <- max.val
            pheatmap( m,  annotation_col=cdesc[ord.idx, rev(c(class.variable, variable.other) )],
                      annotation_colors=cdesc.color,
                      cluster_col=F,
                      gaps_col=cumsum(table(cdesc[ord.idx, 'NMF.consensus'])),
                      main=paste(i),
                      #filename=paste('7_heatmap_gene_of_interest_', sub('\\[.*|\\$','',i), '.pdf', sep=''),
                      filename=paste(sub('\\[.*|\\$','',i), '_heatmap_K', rank, '.pdf', sep=''),
                      cellwidth=cw, cellheight=ch, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")),
                      cluster_row=F)

          } else {
            m <- feat.expr[, rownames(cdesc)[ord.idx]]

            ## heatmap
            m.max <-  ceiling(max(abs(m), na.rm=T))
            #m.max <- max.val
            
            pheatmap( m,  annotation_col=cdesc[ord.idx , rev(c(class.variable, variable.other) )],
                      annotation_colors=cdesc.color, cluster_col=F,
                      gaps_col=cumsum(table(cdesc[ord.idx, 'NMF.consensus'])),main=paste( i),
                      #filename=paste('7_heatmap_gene_of_interest_', sub('\\[.*|\\$','',i), '.pdf', sep=''),
                      filename=paste(sub('\\[.*|\\$','',i), '_heatmap_K', rank,' .pdf', sep=''),
                      cellwidth=cw, cellheight=ch, breaks=seq(-m.max, m.max, length.out=12), color=rev(brewer.pal (11, "RdBu")), cluster_cols = F)
          }

          feat.n <- unlist(lapply(feat, length))

          ##########################
          ## barplot
          #pdf(paste('8_barplot_',sub('\\[.*|\\$','',i),'_features.pdf', sep=''), 3, 3)
          #pdf( paste(sub('\\[.*|\\$','',i), '_barplot_K', rank,'.pdf', sep=''), 3, 3)
          #fancyBarplot(feat.n, srt=45, xlab='NMF basis', names=1:length(s.acc), main=paste(i), ylab='No. features')
          #dev.off()
        }


      }

      ##########################
      ## boxplots
      try(BoxplotNMFmarkers(nmf.marker.str=dir('../', pattern='^NMF_features_N_[0-9]*\\.xlsx', full.names = T),
                        mo.mat.str=dir('../../', pattern='^mo-data-matrix_.*.gct', full.names = T),
                        nmf.map.str=dir('../', pattern='^nmf_vs_.*\\.txt', full.names = T),
                        clust.str=dir('../', pattern='^clin_anno_nmf.txt', full.names = T),
                        nmf.res.dir='../',
                        K=as.numeric(rank),
                        gene_col=opt$gene_col))
      
      setwd('..')
    }

    ## #######################################################################
    ##
    ##                      compare to PCA
    ##
    ## #######################################################################
    pch.vec <- c(15, 19, 17, 18, 8, 4, 3, 2, 1, 5:7, 9:14, 16, 20:25 )


    ## ############################################
    ## PCA on entire matrix
    pca <- prcomp( t(expr.org[, names(NMF.consensus) ]) )
    pca.var <- summary(pca)$importance
    
    pc1 <- pca$x[,1]
    pc2 <- pca$x[,2]

    ## correct colors
    col.tmp <- cdesc.color[[class.variable]]
    col <- as.character(cdesc[names(NMF.consensus), class.variable])
    
    #names(col) <- cdesc[,1]
    col <- sapply(col, function(i) col.tmp[i])

    #################################
    ## plot PDF
    pdf(paste('9.0_PCA.pdf'))
    plot(pc1, pc2, col=col,  cex=2, main=paste('PCA: ', paste(data.omes, collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')), 
         xlab=glue("PC 1 ({round(pca.var[2,1]*100, 1)} %)"), 
         ylab=glue("PC 2 ({round(pca.var[2,2]*100, 1)} %)"), 
         pch=pch.vec[ NMF.consensus ] )
    pointLabel(pc1, pc2, labels=cdesc[ names(NMF.consensus), 1], col=col, offset=50, method='SANN', cex=.5)
    
    legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
    legend('bottomright', legend=paste0("C",sort(unique(NMF.consensus))),
           pch=pch.vec[sort(unique(NMF.consensus))], title='cluster')
    
    dev.off()
    #################################
    ## plot PNG
    png(paste('9.0_PCA.png'))
    plot(pc1, pc2, col=col,  cex=2, main=paste('PCA: ', paste(data.omes, collapse=', '), '\n', paste(paste(dim(expr), collapse=' x '), 'matrix')), 
         xlab=glue("PC 1 ({round(pca.var[2,1]*100, 1)} %)"), 
         ylab=glue("PC 2 ({round(pca.var[2,2]*100, 1)} %)"), 
         pch=pch.vec[ NMF.consensus ] )
    pointLabel(pc1, pc2, labels=cdesc[ names(NMF.consensus), 1], col=col, offset=50, method='SANN', cex=.5)
    
    legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
    legend('bottomright', legend=paste0("C",sort(unique(NMF.consensus))),
           pch=pch.vec[sort(unique(NMF.consensus))], title='cluster')
    
    dev.off()
    
    

    ## #############################################
    ## PCA on NMF features
    expr.nmf <- expr[feat.all, names(NMF.consensus) ]
    pca.nmf <- prcomp( t(expr.nmf) )
    pca.nmf.var <- summary(pca.nmf)$importance

    pc1 <- pca.nmf$x[,1]
    pc2 <- pca.nmf$x[,2]

    col.tmp <- cdesc.color[[class.variable]]
    col <- as.character(cdesc[names(NMF.consensus), class.variable])
    #names(col) <- cdesc[,1]
    col <- sapply(col, function(i) col.tmp[i])

    ########################################
    ## plot: PDF
    pdf(paste('9.1_PCA_NMF_features.pdf'))
    
    plot(pc1, pc2, col=col,  cex=2, 
         main=paste('PCA on NMF features: ', 
                    paste(data.omes, collapse=', '), '\n', 
                    paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')), 
         xlab=glue("PC 1 ({round(pca.nmf.var[2,1]*100, 1)} %)"), 
         ylab=glue("PC 2 ({round(pca.nmf.var[2,2]*100, 1)} %)"), 
         pch=pch.vec[ NMF.consensus ])
    pointLabel(pc1, pc2, labels=cdesc[ names(NMF.consensus), 1], col=col, offset=20, method='SANN', cex=.5)
    legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
    legend('bottomright', legend=paste0('C',sort(unique(NMF.consensus))),
           pch=pch.vec[sort(unique(NMF.consensus))], title='cluster')
    dev.off()

    ## plot: PNG
    png(paste('9.1_PCA_NMF_features.png'))
    plot(pc1, pc2, col=col,  cex=2, 
         main=paste('PCA on NMF features: ', 
                    paste(data.omes, collapse=', '), '\n', 
                    paste(paste(dim(expr.nmf), collapse=' x '), 'matrix')), 
         xlab=glue("PC 1 ({round(pca.nmf.var[2,1]*100, 1)} %)"), 
         ylab=glue("PC 2 ({round(pca.nmf.var[2,2]*100, 1)} %)"), 
         pch=pch.vec[ NMF.consensus ])
    pointLabel(pc1, pc2, labels=cdesc[ names(NMF.consensus), 1], col=col, offset=20, method='SANN', cex=.5)
    legend('topright', legend=names(col.tmp), pch=16, cex=1, col=col.tmp, bty='n')
    legend('bottomright', legend=paste0('C',sort(unique(NMF.consensus))),
           pch=pch.vec[sort(unique(NMF.consensus))], title='cluster')
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
        tsne.nmf.col[ s.acc[[i]]$Accession[ s.acc[[i]]$Accession %in% names(tsne.nmf.col) ] ] <- NMF.consensus.col[ i ] %>% as.character #unlist(NMF.basis.col[as.character(i)])

      ## just the NMF class
      tsne.nmf.class <- rownames(expr.nmf)
      names(tsne.nmf.class) <- tsne.nmf.class
      for(i in 1:length(s.acc))
        tsne.nmf.class[ s.acc[[i]]$Accession[ s.acc[[i]]$Accession %in% names(tsne.nmf.class)] ] <- paste(i,  cons.map.max[i], sep=' ')



      d <- data.frame(
        tsne_1=tsneY[,1],
        tsne_2=tsneY[,2],
        Data_type=sub('-.*', '', rownames(expr.nmf)),
        Feature=rownames(expr.nmf),
        NMF_col=tsne.nmf.col,
        NMF_class=tsne.nmf.class
      )
      
      ## symbols for data type
      symbs <- c('circle', 'x', 'o', 'square', 'diamond', 'triangle-down', 'pentagon')

      ## colors
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



###################################################
# multi-omics boxplots based on nmf clustering
BoxplotNMFmarkers <- function(nmf.marker.str, ## path to Excel sheet containing markers
                              mo.mat.str,     ## path to multi-omic data matrix
                              nmf.map.str,    ## path to 'nmf_vs_XX.txt' file
                              clust.str,      ## path to 'clin_anno_nmf.txt'
                              K,              ## number of clusters
                              nmf.res.dir,    ##
                              gene_col        ## column name containing gene symbols
                              ){

  p.clust.enrich <- 0.01

  ## import clustering
  clust <- read.delim(clust.str, row.names=1)

  ## import data matrix
  mo <- parse.gctx(mo.mat.str)
  mo.mat <- mo@mat
  mo.rdesc <- mo@rdesc
  mo.cdesc <- mo@cdesc

  # nmf2class vectors
  nmf.map.str <- dir(nmf.res.dir, pattern='^nmf_vs_.*\\.txt', full.names = T)

  nmf.map <- read.delim(nmf.map.str, stringsAsFactors = F, row.names = 1)
  nmf.map <- apply(nmf.map, 1, function(x) paste(names(x)[x <  p.clust.enrich], collapse='|' ) )
  idx.tmp <- which(nchar(nmf.map) == 0)
  nmf.map[idx.tmp] <- names(nmf.map)[idx.tmp]

  ## make mappings unique
  nmf.map.names <- names(nmf.map)
  nmf.map <- make.unique(nmf.map)
  names(nmf.map) <- nmf.map.names
  
  ## ########################################
  ## import nmf markers
  sheets <- excel_sheets(nmf.marker.str)
  nmf.map <- nmf.map[ sheets ]
  
  #nmf.marker.list <- lapply(1:K, function(i) read_xlsx(nmf.marker.str, sheet = i) %>% as.data.frame )
  nmf.marker.list <- lapply(sheets, function(i) read_xlsx(nmf.marker.str, sheet = i) %>% as.data.frame )
  names(nmf.marker.list) <- nmf.map

  #valid.idx <- which(sapply(nmf.marker.list, nrow) > 0)
  
  goi <- lapply(nmf.marker.list, function(x)x$SYMBOL) %>% unlist %>% paste(., collapse=';')
  goi <- goi %>% strsplit(., ';') %>% unlist %>% unique

  ## signed scores
  # nmf.marker.list <- lapply(nmf.marker.list[valid.idx], function(x) x %>% mutate(Score = if_else(Direction == 'down', -1*Score, Score)))
  nmf.marker.list <- lapply(nmf.marker.list, function(x) 
    if(nrow(x) > 0)
      x %>% mutate(NMF.Score = if_else(Direction == 'down', -1*NMF.Score, NMF.Score))
    )
  
  ## all 'omics data types
  omes.all <- sub('^(.*?)-.*','\\1',mo.rdesc$Data.Type.ID) %>% unique

  ## all nmf cluster
  clust.all <- names(nmf.marker.list)
  clust2samples <- lapply(nmf.map, function(x) rownames(clust)[clust$NMF.consensus == as.numeric(names(nmf.map)[nmf.map == x])])
  names(clust2samples) <- nmf.map

  for(g in goi){

    ## goi in data set?
    g.mo.idx <- grep(glue("^{g}$"), mo.rdesc[, gene_col])
    
    if( length( g.mo.idx ) > 0 ){

      ## all occurences
      mo.mat.g <- matrix(mo.mat[g.mo.idx, ], nrow=length(g.mo.idx))
      rownames(mo.mat.g) <- mo@rid[g.mo.idx]
      colnames(mo.mat.g) <- colnames(mo.mat)

      ## split by ome
      g.omes <- lapply(omes.all, function(x){
        idx=grep(glue("^{x}-.*"), rownames(mo.mat.g))
        if(length(idx) > 0)
          matrix( mo.mat.g[idx, ], nrow=length(idx), dimnames=list(rownames(mo.mat.g)[idx], colnames(mo.mat.g)) )
      })
      names(g.omes) <- omes.all

      ## remove empty elements
      keep.tmp <- which( sapply(g.omes, length) > 0)
      g.omes <- g.omes[keep.tmp]

      ## split by NMF cluster
      g.omes.nmf <- lapply(g.omes, function(x){
        lapply(clust2samples, function(xx) matrix(x[, xx], nrow=nrow(x), dimnames=list(rownames(x), xx)) )
      })

      pdf(glue("{g}_boxplot_K{K}.pdf"))

      #########################################
      ## loop over omes
      for(o in names(g.omes)){

        og <- g.omes[[o]]

        ## loop over each feature
        feat <- rownames(og)
        for(f in feat){

          ## split into nmf cluster
          og.nmf <-  lapply(clust2samples, function(x) og[f, x] )

          ## feature among NMF markers?
          nmf.idx <- sapply( nmf.marker.list, function(x) {
            xo <- x[x$Type == o, ]
            sub('^.*?-','',f) %in% xo$Accession
          })

          col.tmp=rep('black', length(og.nmf))

          #title <- glue("{o}-{g} ({sub('.*?-','',f)})")
          title <- glue("{g}\n{f}")
          nmf.score <- NA
          if(sum(nmf.idx) > 0){
            nmf.idx <- which(nmf.idx == 1)[1]
            col.tmp[nmf.idx] <- 'red'

            ## feature score
            nmf.tmp <- nmf.marker.list[[nmf.idx]]
            nmf.tmp <- nmf.tmp[nmf.tmp$Type == o, ]
            nmf.tmp <- nmf.tmp[which(nmf.tmp$Accession == sub('^.*?-','',f)), ]
            #title <- glue("{title}\n{nmf.tmp$GENENAME}")
            nmf.score <- glue("score={round(nmf.tmp$NMF.Score,3)}")
            nmf.modt.adjp <- glue("adj.p={format(nmf.tmp[, grep('^adj.P', colnames(nmf.tmp))], digits=3)}")
          }

          fancyBoxplot(og.nmf, main=title, vio.alpha = 0, show.numb = 'numb.bottom',numb.cex = 1,
                       col = 'white', lwd=3, box.border=col.tmp, ylab=glue("{o}"), las=2)
          abline(h=0, lty='dashed', col='grey', lwd=2)
          if(!is.na(nmf.score))
            legend('topright', legend=c(nmf.score, nmf.modt.adjp), text.col='red')
        }

      } ## end loop over omes
      dev.off()
    }
  }

  return(0)
}


