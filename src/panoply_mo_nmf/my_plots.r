#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
## 20151013 collection of plotting functions
library(pacman)

# 
# OS <- Sys.info()['sysname']
# if(OS == 'Windows')
#   url.stem <- 'C:/Users/karsten/Dropbox/Devel/R-code/'
#   ##  url.stem <- '//flynn-cifs/prot_proteomics/'
# if(OS == 'Linux')
#   url.stem <- '~/karsten/'
# if(OS == 'Darwin')
#   url.stem <- '/Volumes/prot_proteomics/'
#   
# 
# 
# source(paste(url.stem, "misc.r", sep='/') )

##source(paste("c:/Users/Karsten/Dropbox/Devel/Code/misc.r", sep='/') )

## ###########################################################################
## ssGSEA volcano
##
##
gsea_volc <- function(x, y, s=NULL, pch.cex=NULL, ids=NULL, l=NULL, lcol=NULL, fdr.max=0.05, alpha=150, ...){
  # x - logFC
  # y - log p-val
  # s - index of significant features
  # pch.cex - point size
  # ids - unique ids for features
  # l - named list of characters of labels
  # l.col - names list (same as l) defining colors
  # fdr.max -only used as label in legend; ought to be consistent to 's'
  # alpha - passed to 'my.col2rgb'
  
  col=my.col2rgb('grey', alpha)
  col.b <- my.col2rgb('grey 20', alpha) # border
  
  col.dn <- my.col2rgb('lightblue', alpha)
  col.up <- my.col2rgb('coral1', alpha)
  
  pch <- 21
  pch.org <- pch
  #pch.b <- 21 # border
  
  pch.sig <- 23
  pch.sig.org <- pch.sig
  #pch.sig.b <- 23 # border
  
  if(is.null(pch.cex)) pch.cex <- 1
  
  ## pch and col
  if(!is.null(s)){
    # color points
    col=rep(col, length(x))
    col[ s[x[s] > 0] ] <- col.up
    col[ s[x[s] < 0] ] <- col.dn
    
    # border border
    col.b=rep(col.b, length(x))
    #col.b[ s[x[s] > 0] ] <- col.up
    #col.b[ s[x[s] < 0] ] <- col.dn
    
    ## pch
    pch=rep(pch, length(x))
    pch[s] <- rep(pch.sig, length(s))
    ## pch border
    #pch.b=rep(pch.b, length(x))
    #pch.b[s] <- rep(pch.sig.b, length(s))
  }
  ## map labels
  if(is.null(ids))
    ids <- 1:length(x) %>% as.character()
  
  ## plot  
  plot(x, y, col=col.b, bg=col, pch=pch, cex=pch.cex, xlab=expression(log(FC)), ylab=expression(-10*log[10](p-value)), ...)
  abline(v=0, lty='dashed', lwd=2, col='grey')
  #points(x, y, col=col.b, bg=col.b ,pch=pch.b, cex=pch.cex)
  
  ## labels
  if(!is.null(l)){
    ## loop over label groups
    for(j in 1:length(l)){
      ## match to 'ids'
      l.idx <- grep(l[[j]], ids)
      l.lab <- sub('.*?_','',grep(l[[j]], ids, value=T))
      ## label if found
      if(length(l.idx) > 0){
        
        ## right side
        l.idx.up <- l.idx[ which(x[l.idx] > 0) ] 
        l.lab.up <- l.lab[ which(x[l.idx] > 0) ]
        if(length(l.idx.up)>0)
          #text(x=xlim[2], y=y[l.idx.up], labels = l.lab.up, cex=0.6, pos=2, col=lcol[[j]])
          text(x=x[l.idx.up]+0.5, y=y[l.idx.up], labels = l.lab.up, cex=0.6, pos=4, col=lcol[[j]])
        ## left side
        l.idx.dn <- l.idx[which(x[l.idx] < 0)]
        l.lab.dn <- l.lab[ which(x[l.idx] < 0) ]
        if(length(l.idx.dn))
          text(x=x[l.idx.dn]-0.5, y=y[l.idx.dn], labels = l.lab.dn, cex=0.6, pos=2, col=lcol[[j]])
        #text(x=xlim[1], y=y[l.idx.dn], labels = l.lab.dn, cex=0.6, pos=4, col=lcol[[j]])
        
        points(x[l.idx], y[l.idx], col=lcol[[j]], bg=my.col2rgb(lcol[[j]], alpha), pch=pch[l.idx], cex=pch.cex[l.idx])
      }
    }
  }
  
  legend('topright', legend=names(lcol), col=unlist(lcol), pch=20, cex=0.8, pt.cex=2, bty='n')
  legend('topleft', legend=c(paste('FDR >=', fdr.max), paste('FDR <', fdr.max)), col='black', pch=c(pch.org, pch.sig.org), bty='n', cex=0.8, pt.cex=1.6)
  legend('top', legend=c('100%', '50%', '10%') , pch=21, pt.cex=c(2, 1, 0.2), bty='o', title='Signature coverage', ncol = 3, cex=.8, bg='white', box.col='white')
}
## ###############################################################
## gsea_volc_ly
gsea_volc_ly <- function(fc, logP, fdr, id, perc, perc.scale=3, fdr.max, main=''){
  require(plotly)
  dat <- data.frame(fc, logP, fdr, id, perc, 
                    cex=perc/perc.scale,
                    info=paste(id,'\nFDR: ', round(fdr,3), '\n% detect: ', perc, sep=''))
  sig <- which(fdr < fdr.max)
  dat.bg <- dat[-sig, ]
  dat.sig <- dat[sig, ]
  
  p <- plot_ly(x=dat.bg$fc, y=dat.bg$logP, text=dat.bg$info, 
               type='scatter', mode='markers', marker=list(size=dat.bg$cex, color='grey'), name=paste('FDR>=', fdr.max, sep='') ) 
  p <- p %>% add_trace(x=dat.sig$fc, y=dat.sig$logP, text=dat.sig$info, type='scatter', mode='markers', marker=list(size=dat.sig$cex, color='red'), name=paste('FDR<', fdr.max, sep=''))
  p <- p %>% layout(title=main)
  return(p)
}

#################################################################################################
##            multiscatterplot using hexagonal binning
## - mat    numerical matrix of expression values, rows are features, columns are samples
##
## changelog: 2015116 implementation
#################################################################################################
my.multiscatter <- function(mat, hexbin=30, hexcut=5, cor=c('pearson', 'spearman', 'kendall')){

    size.title.hist=10
    size.cor=35
    xlim=c(-10, 10)
    ylim=c(-10, 10)
    
    p_load(hexbin)
    p_load(ggplot2)
    p_load(Hmisc)
    p_load(grid)
    
    ###########################################################################
    ## original code from:  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
    multiplot <- function(plots, cols=1) {
      ## Make a list from the ... arguments and plotlist
      ##plots <- c(list(...))
      ## number of plots
      numPlots = length(plots)
      ## layout matrix
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
      ## Set up the page
      grid.newpage()
      ## grid layout
      la <-  grid.layout(nrow(layout), ncol(layout))
      pushViewport(viewport(layout = la))
      ## Make each plot, in the correct location
      for (i in numPlots:1) {
        ## Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col)
        
        ## textplot: correlation coefficient
        if(matchidx$row < matchidx$col){
          numb = plots[[i]]
          ##col=
          ##size = max(abs(90*as.numeric(numb)), 30)
          size=size.cor
          grid.rect(width=unit(.85, 'npc'), height=unit(.85, 'npc'), vp=vp, gp=gpar(fill='grey95', col='transparent'))
          grid.text(numb, vp=vp, gp=gpar(fontsize=size))
        } else {
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    } ## end function 'multiplot'
    #########################################################################
    
    ## cor method
    corm = match.arg(cor)
    ## correlation
    cm = cor(mat, use='pairwise.complete', method=corm)
    
    ## number of samples to compare
    N = ncol(mat)
    
    ## list to store the plots
    plotList=vector('list', N*N)
    count=1
    for(i in 1:N)
      for(j in 1:N){
        
        dat <- data.frame(x=mat[,i], y=mat[,j])
        
        ###########################
        ## lower triangle
        if(i < j){
          
          ## hexbin
          hex <- hexbin(dat$x, dat$y, hexbin)
          gghex <- data.frame(hcell2xy(hex), c = cut2(hex@count, g = hexcut))
          
          p <- ggplot(gghex) + geom_hex(aes(x = x, y = y, fill = c), stat = "identity") + guides(fill=FALSE) + scale_fill_manual( values=paste('grey', ceiling(seq(70, 20, length.out=hexcut)), sep='')) + xlab('') + ylab('') + xlim(xlim[1], xlim[2]) + ylim(ylim[1], ylim[2])
          
          marg=c(0, -.25, -.4, -.25)
          
          if(i == 1 & j < N){
            ##p = p + theme( plot.margin=unit(rep(-.2, 4), 'cm'), axis.ticks.x=element_blank(), axis.text.x=element_blank() )
            p = p + theme( plot.margin=unit(marg, 'cm'), axis.ticks.x=element_blank(), axis.text.x=element_blank() )
          } else if(i == 1 & j == N){
            ##p = p + theme( plot.margin=unit(rep(-.2, 4), 'cm') )
            p = p + theme( plot.margin=unit(marg, 'cm') )
          } else if(i > 1 & j == N) {
            ##                    p = p + theme( plot.margin=unit(rep(-.2, 4), 'cm'), axis.ticks.y=element_blank(), axis.text.y=element_blank()  )
            p = p + theme( plot.margin=unit(marg, 'cm'), axis.ticks.y=element_blank(), axis.text.y=element_blank()  )
          } else if(i > 1 & j < N) {
            ##p = p + theme( plot.margin=unit(rep(-.2, 4), 'cm'), axis.ticks=element_blank(), axis.text=element_blank()  )
            p = p + theme( plot.margin=unit(marg, 'cm'), axis.ticks=element_blank(), axis.text=element_blank()  )
          }
          
          
        }
        ###########################
        ## diagonal
        if(i == j){
          p = ggplot(dat, aes(x=x)) + geom_histogram(fill='grey70', colour='black', binwidth=sum(abs(range(dat$x, na.rm=T)))/50) + ggtitle(colnames(mat)[i]) + theme(plot.title=element_text(size=size.title.hist)) + theme( panel.background = element_blank(), plot.margin=unit(rep(0, 4), 'cm')) + ylab('') + xlim(xlim[1], xlim[2]) + xlab('') ## + xlab(paste('N',sum(!is.na(dat$x)), sep='=')) + annotate('text', label=sum(!is.na(dat$x)), x=unit(0, 'npc'), y=unit(0, 'npc'))
          
        }
        ###########################
        ## upper triangle
        if(i > j){
          cortmp = cm[i,j]
          
          p=paste(round(cortmp,2))
          
        }
        
        plotList[[count]] <- p
        count=count+1
      }
    
    multiplot( plotList, cols=N)
}


##########################################################################################################
## translate a color name into rgb space
##
## changelog:  20100929 implementation
##########################################################################################################
my.col2rgb <- function(color, alpha=80, maxColorValue=255){
  
  out <- vector( "character", length(color) )
  
  for(col in 1:length(color)){
    
    col.rgb <- col2rgb(color[col])
    
    out[col] <- rgb(col.rgb[1], col.rgb[2], col.rgb[3], alpha=alpha, maxColorValue=maxColorValue)
    
  }
  return(out)
}

#####################################################
##
##                   fancyBarplot
##
## a more adjustable version of 'barplot'
##
## 20110217 added support for matrices
## 20140201 'cex.numb'
#####################################################
fancyBarplot <- function(x, space=0.2, 
                         ylim=c(0, max(x, na.rm=T)+ 0.15*max(x, na.rm=T)), 
                         ndigits=3, 
                         add.numb=T,
                         numb=c('counts', 'dev.median.perc'),
                         cex.numb=.9, srt=90, border=NA, ...)
{
  numb <- match.arg(numb)
  
  #########################################
  #         vector
  #########################################
  if(length(dim(x)) <= 1){
    # make the plot
    barplot(x, space=space, ylim=ylim, border=border, ...)
    
    # add numbers
    if(add.numb){
      if(numb=='counts'){
        numb.to.add <- x
        ## add to plot
        text(  seq( 0.5+space, (0.5+space)+((length(numb.to.add)-1)*(1+space)), 1+space ),
               x+(0.05*max(x)), round(numb.to.add, ndigits), srt=srt, pos=4, offset=-0.05, cex=cex.numb )
      }
      
      if(numb == 'dev.median.perc'){
        numb.median <- median(x, na.rm=T)
        numb.to.add <- sapply(x, function(y) 100*(y-numb.median)/numb.median )
        abline(h=numb.median, col='darkblue', lty='dashed')
        #legend('topleft', legend=c(paste('Median:', round(numb.median, ndigits))), bty='n', text.col = 'darkblue')
        mtext(c(paste('Median:', round(numb.median, ndigits))),side = 2, at = numb.median, col = 'darkblue', cex = 0.8)
        ## add to plot
        text(  seq( 0.5+space, (0.5+space)+((length(numb.to.add)-1)*(1+space)), 1+space ),
               x+(0.05*max(x)), paste(round(numb.to.add, ndigits),'%', sep=''), srt=srt, pos=4, offset=-0.05, cex=cex.numb )
      }
        }
      
    } else{
    
    #########################################
    #        matrix
    #########################################
    offset = -0.1
    barplot(x, ylim=ylim, beside=T, border=border, ...)
    
    if(add.numb){
      p.tmp=1
      for(p in 1:dim(x)[2]){
        if( sum(is.na(x[, p])) < nrow(x) ){
          text( (p.tmp : (p.tmp + dim(x)[1]-1) ) + offset, max(x[, p], na.rm=T)+0.05*( max(x[, p], na.rm=T)), round(x[,p], ndigits), adj=c(1,1), pos=4, cex=cex.numb, srt=srt  )
        }
        p.tmp <- p.tmp+dim(x)[1]+1
      }
    }
    
  }
}


############################################################################################
## code written by mani dr
##
##
## Scatter plot with marginal histograms, prediction confidence interval calculation
## (plus other features ... see below)
scatterhist.ci <- function (x, y, id=NULL, xlab="", ylab="", title="",
                            scatterplot.color='black', scatterplot.transparency=1, cex=1,
                            marginal.histograms=TRUE, draw.axes=FALSE,
                            special.subsets=NULL, subset.colors=NULL,
                            subset.transparencies=NULL, subset.cex=NULL,
                            limits=NULL, n.bars=100, pch=19, regr.type='linear',
                            plot.foldchange=NULL, regr=FALSE, ci.level=0.99,
                            regr.line.color='blue', ci.line.color='magenta',
                            ...) {
  #
  # 1. combines scatter plot + histogram plotting function
  #      a scatter plot of (x,y) is created, with marginal histograms for
  #      each of the x and y data points
  # 2. if special.subset is specified (a boolean vector of the same length
  #      as x and y, or a list of indices), those points are replotted in subset.color
  # 3. can calculate / plot linear or loess regression line, along with predcition CIs
  # 4. plot.foldchange: if specified as a single number, that limit is plotted in x and y
  #                     if a pair is specified, the x and y limits are plotted accordingly
  #                     (units identical to x and y)
  # Returns a table with (id, x, y) in addition to predict fit, ci and indication of
  #     whether the point lies with in the precition confidence interval
  #
  # Some extra features:
  # i.  In addition to the main plot, multiple subsets can be specified (as an index list
  #     using special.subsets. Color of the main plot and subsets can be directly
  #     specified (using
  #     scatterplot.color and subset.color
  # ii. scatterplot.transparency and subset.transparencies can be set to control
  #     transparency of main and subset plot points, respectively
  #     (semi-transparent points provide a better indication of how many points overlap)
  # iii. Similar control of relative size of point is provided by cex and subset.cex
  # iv. marginal.histograms=FALSE will disable plotting marginal histograms
  #
  
  
  
  if ( length(x) != length(y) ) stop ('x and y must contain the same number of points')
  
  
  if (is.null (id)) {
    # use row numbers as id if not specified
    id <- 1:length(x)
  }
  
  if (is.null (limits)) {
    # use full range of x, y if limits not specified
    r <- c (range (x, na.rm=TRUE), range (y, na.rm=TRUE))
    limits <- c ( min(r), max(r) )
  }
  
  
  keep <- x >= limits[1] & x <= limits[2] & y >= limits[1] & y <= limits[2] & is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  id <- id[keep]
  
  
  return.value <- data.frame (id,x,y)
  
  if (marginal.histograms) {
    zones <- matrix (c(2,0,1,3), ncol=2, byrow=TRUE)
    layout (zones, widths=c(7/10,3/10), heights=c(3/10,7/10))
  }
  
  xhist <- hist (x, breaks=seq(limits[1], limits[2], length.out=n.bars), plot=FALSE)
  yhist <- hist (y, breaks=seq(limits[1], limits[2], length.out=n.bars), plot=FALSE)
  top <- max (c(xhist$counts, yhist$counts))
  
  if (marginal.histograms) par (mar=c(5,6,1,1))
  else par (mar=c(5,6,2,2))
  plot (x, y, xlab=xlab, ylab=ylab, xlim=limits, ylim=limits, pch=pch,
        col=alpha (scatterplot.color, scatterplot.transparency), cex=cex, ...)
  if (draw.axes) {
    abline (h=0)
    abline (v=0)
  }
  if (!is.null (special.subsets)) {
    n <- length (special.subsets)
    # if a list of subset vectors is provided, re-plot those points using subset.colors
    special.subsets <- lapply (special.subsets, function (s) { s[keep] })
    # create colors if necessary; if colors provided, recycle if not sufficient
    if (is.null (subset.colors)) subset.colors <- brewer.pal (n, 'Dark2')
    else subset.colors <- rep (subset.colors, length.out=n)
    # create transparency if necessary; if provided, recycle if not sufficient
    if (is.null (subset.transparencies)) subset.transparencies <- rep (1, length.out=n)
    else subset.transparencies <- rep (subset.transparencies, length.out=n)
    # create cex vector if necessary; if provided, recycle if not sufficient
    if (is.null (subset.cex)) subset.cex <- rep (1, length.out=n) * cex
    else subset.cex <- rep (subset.cex, length.out=n) * cex
    
    for (i in 1:n)
      points ( x[special.subsets[[i]]], y[special.subsets[[i]]], pch=pch,
               col=alpha (subset.colors[i], subset.transparencies[i]), cex=subset.cex[i], ... )
  }
  if ( !is.null (plot.foldchange) ) {
    # plot foldchange lines at specified (x and/or y) value
    if ( length (plot.foldchange)!=2 ) f.x <- f.y <- plot.foldchange [1]     # use only first value
    else {
      # plot.foldchange = (x,y)
      f.x <- plot.foldchange [1]
      f.y <- plot.foldchange [2]
    }
    points (c(-f.x,-f.x), limits, type='l', lty=2, col=4)
    points (c(f.x,f.x), limits, type='l', lty=2, col=4)
    points (limits, c(f.y,f.y), type='l', lty=2, col=4)
    points (limits, c(-f.y,-f.y), type='l', lty=2, col=4)
    rect (-f.x,-f.y,f.x,f.y, border=4)
  }
  if (regr) {
    x.values <- seq (limits[1]/0.9, 1.1*limits[2], (limits[2]-limits[1])/100)
    
    if (regr.type=='linear') {
      # build linear regression model
      m <- rlm (y ~ x, method='MM', weight=abs(x)+abs(y), wt.method='case')
      print (summary (m))
      # plot model and prediction confidence interval
      abline (m, col=regr.line.color)
      m.ci <- predict (m, data.frame (x=x.values), interval='prediction', level=ci.level)
      points (x.values, m.ci[,'lwr'], type='l', col=ci.line.color)
      points (x.values, m.ci[,'upr'], type='l', col=ci.line.color)
      
      # derive prediction and CI for input data points
      x.pred <- predict (m, data.frame (x=x), interval='prediction', level=ci.level)
      return.value <- cbind (return.value, x.pred)
    }
    
    if (regr.type=='loess') {
      m <- loess (y ~ x, degree=1, span=0.5, family='symmetric', surface='direct')
      m.pred <- predict (m, data.frame(x=x.values), se=TRUE)
      m.se <- sqrt ( m.pred$se.fit^2 + as.numeric (m$s)^2 )
      m.lwr <- m.pred$fit + qt ((1-ci.level)/2, df=as.numeric (m$enp)) / 2 * m.se
      m.upr <- m.pred$fit + qt (1 - (1-ci.level)/2, df=as.numeric (m$enp)) / 2 * m.se
      lines (x.values, m.pred$fit, col=regr.line.color)
      lines (x.values, m.lwr, col=ci.line.color)
      lines (x.values, m.upr, col=ci.line.color)
      m.ci <- data.frame (fit=m.pred$fit, lwr=m.lwr, upr=m.upr)
      
      # derive prediction and CI for input data points
      x.pred <- predict (m, data.frame (x=x), se=TRUE)
      x.se <- sqrt ( x.pred$se.fit^2 + as.numeric (m$s)^2 )
      x.lwr <- x.pred$fit + qt ((1-ci.level)/2, df=as.numeric (m$enp)) / 2 * x.se
      x.upr <- x.pred$fit + qt (1 - (1-ci.level)/2, df=as.numeric (m$enp)) / 2 * x.se
      return.value <- cbind (return.value, fit=x.pred, lwr=x.lwr, upr=x.upr)
    }
    # determine points which lie in the confidence interval limit
    within.interval <- with (return.value, y >= lwr & y <= upr)
    return.value <- cbind (return.value, within.interval)
  }
  
  if (marginal.histograms) {
    par (mar=c(0,6,3,1))
    barplot (xhist$counts, axes=T, ylim=c(0, top), space=0, main=title)
    text (limits[1]+1, top/3, paste ('mean=', round (mean(x),digits=2), '\nsd=', round (sd(x),digits=2), sep=''), pos=4)
    par (mar=c(5,0,1,3) )
    barplot (yhist$counts, axes=T, xlim=c(0, top), space=0, horiz=TRUE)
    text (top/3, limits[2]+1, paste ('mean=', round (mean(y),digits=2), '\nsd=', round (sd(y),digits=2), sep=''), pos=4, srt=90)
  }
  
  invisible (return.value)
}




##########################################################################
##
##                                 fancyPlot
##
##  produce fancy scatterplot
##
##
## x
## y
## rug
## grid
## loess
## cor
##
##
## changelog:
## 20090703     - implementation
## 20100216     - checking missing values also for NAN and infinity
## 20100311     - added transparancy
## 20100325     - number of datapoints is now reported
## 20100510     - added parameter 'groups'
## 20101011     - changed parameter 'loess' to 'reg' and
##                added linear regression
##              - parameter 'cex.text': controls 'cex' for lagend texts like
##                correlation or regression
## 20110615     - changed color vector: in previous R version one had to
##        ?       supply a color vector starting with an arbitrary color followed
##                by the actual colors (function: scatterpot)
##              - now the color vector starts with the actual color
## 20121004     - included 'reset.par' ain the argument list
##              - before it was set to FALSE as default
## 20130523     - included missing values using function 'rug'
## 20140122     - parameter 'reg.col' to specify color of regression line
##              - parameter 'reg.legend'
##              - loess regression plottted as lines
##
############################################################################
fancyPlot <- function(x, y, rug=F, grid=T, reg=c("none", "linear", "loess"), 
                      reg.col="darkred", reg.legend=T, cor=T, 
                      cor.method = c( "pearson", "spearman", "kendall"), pch=19, col="black", 
                      alpha=60, boxplots=c("xy","x", "y", "n"), groups=1, cex.text=1, reset.par=F, cex.pch=1, ...)
{
  p_load(car)
  
  # regression
  reg = match.arg(reg)
  
  ############################
  # color vector?
  ############################
  col.asis = F
  if(length(col) == length(x))
    col.asis = T
  
  ############################
  # color
  ############################
  if(col.asis){
    
    color <- list()
    
    for(i in 1:length(col)){
      color[[i]] <- col2rgb(col[i])
      color[i] <- rgb(color[[i]][1], color[[i]][2], color[[i]][3], alpha, maxColorValue=255)
    }
    #col=c("black", unlist(color))
    col=c(unlist(color)) # 20110615
    #return(col)
  }
  
  # plotting characters: the same number as groups
  if(length(pch) == 1 & length(unique(groups)) > 1) pch <- rep(pch,length(unique(groups)))
  
  #############################
  # correlation coefficient
  #############################
  cor.method=match.arg(cor.method)
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # store orginal data
  x.org <- x
  y.org <- y
  
  ##############################
  # remove missing values
  ##############################
  #na.idx <- union(union( which(is.na(x)), which(is.na(y)) ), union( which(is.nan(x)), which(is.nan(y)) ))
  na.idx <- union( which(is.na(x)), which(is.na(y)) )
  if(length(na.idx) > 0){
    x <- x[-na.idx]
    y <- y[-na.idx]
    
    if(length(groups)>1) groups <- groups[-na.idx]
    if(col.asis) col <- col[-na.idx]
  }
  inf.idx <- union( which(is.infinite(x)), which(is.infinite(y))  )
  if(length(inf.idx)>0){
    x <- x[-inf.idx]
    y <- y[-inf.idx]
    if(length(groups) > 1) groups <- groups[-inf.idx]
    if(col.asis) col <- col[-inf.idx]
  }
  
  #if(col.asis) col <-  c("black", col)
  #if(col.asis) col <-  c("black", col) # 20110615
  if(col.asis) col <-  c("#FFFFFF", col) # 20180411
  
  ###############################
  # boxplot on the axes
  ###############################
  boxplots=match.arg(boxplots)
  
  
  palette.org <- palette()
  ########################################
  #            plot
  #########################################
 # rver <- Rversion() 
  rver <-  paste(version$major, version$minor, sep='.')
  if(rver != "3.5.0")
    scatterplot(x,y, col=col, pch=pch, boxplots=boxplots, reg.line=F, smooth=F, sub=paste("N=", length(x),sep=""), groups=groups, reset.par=reset.par ,legend.plot=ifelse(length(groups) == 1, F, T), cex=cex.pch, ...)
  else{
    
    scatterplot(x,y, col=col, pch=pch, boxplots=boxplots, regLine=F, smooth=F, sub=paste("N=", length(x),sep=""), reset.par=reset.par , cex=cex.pch, ...)
  }
  
  # grid
  if(grid) grid()
  
  # all points
  if(rug) rug(x, col="indianred")
  
  # missing values
  rug(x.org[is.na(y.org)], side=1, col=my.col2rgb("darkblue"), ticksize=.02)
  rug(y.org[is.na(x.org)], side=2, col=my.col2rgb("darkblue"), ticksize=.02)
  
  # regression
  if(reg == "loess"){
    LOESS <- loess( y~x  )
    lines(x[order(x)], LOESS$fitted[order(x)], col=reg.col, type="l",  lwd=2)
    
    if(reg.legend)
      legend("topright", legend="loess", col=reg.col, lwd=2, bg="white", cex=cex.text)
    
  } else if(reg == "linear"){
    
    fit = lm( formula=y~x)
    abline(fit, col=reg.col, lwd=1.5)
    
    if(reg.legend)
      legend("bottomright", legend=paste("y = ", round(fit$coefficient[2], 3), "x ", ifelse(fit$coefficient[1] < 0, "- ", "+ "), round(abs(fit$coefficient[1]), 3), sep="" ) , bty="n", text.col=reg.col, cex=cex.text)
  }
  
  
  # correlation coefficient
  if(cor){
    minY <- min(y, na.rm=T)
    legend("top", legend=paste( toupper(substr(cor.method, 1,1)), substring(cor.method, 2) ,"-correlation = ",  round(cor(x,y, method=cor.method), 2), sep=""), bty="n", cex=cex.text)
    
  }
  #return(col)
}



##############################################################################
#
#  - perform principle component analysis
#  - calculate variances explained by components
#  - plot the results
#
# changelog: 20131001 implementation
#            20131007 3D plot now plot pc1 vs. pc3 vs. pc2
#                     instead of pc1 vs. pc2 vs. pc3
##############################################################################
my.prcomp <- function(x, col=NULL, cor=T, plot=T, rgl=F, scale=T, pch=20, cex.points=3, rgl.point.size=30, main="PCA", ...){
  
  # color
  if( is.null(col) ) col="black"
  
  # perform pca
  pca <- prcomp(x, scale=scale)
  
  # calculate variance
  comp.var <- eigen(cov(pca$x))$values
  
  
  ##############
  # rgl plot
  ##############
  if(rgl){
    p_load(rgl)
    plot3d(pca$x[,1], pca$x[,2], pca$x[,3], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), type="s", col=col, expand=1.2, size=rgl.point.size)
  }
  
  ###############
  # scatterplot 2D/3D
  ###############
  if( plot){
    p_load(scatterplot3d)
    
    par(mfrow=c(1,2))
    
    ## PC 1-2
    plot(pca$x[,1], pca$x[,2], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), pch=pch, main=main, col=col, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:2]/sum(comp.var)),1),"%", sep=""), cex=cex.points )
    
    
    ## PC 1-3
    scatterplot3d( pca$x[,1], pca$x[,3], pca$x[,2], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), color=col,  cex.symbols=cex.points, pch=pch, main=main, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:3]/sum(comp.var)),1),"%", sep=""), type="h" )
    
    par(mfrow=c(1,1))
  }
  
  return(pca)
}


#######################################################
##
##
## changelog: 20140722 implementation
##            20141104 fixed names for x-axis
##            20150209 parameter 'grid'
##            20151027 parameter 'grid.at'
##            20160308 replaced axes by xaxis and yaxis
##            20170303 'vio.wex'
#######################################################
fancyBoxplot <- function(x, 
                         ylim=NULL, 
                         xlim=NULL, 
                         at=NULL, 
                         col="grey80", 
                         vio.alpha=100, 
                         box.border="black", 
                         box.pch=20, 
                         drawRect=F, 
                         xlab="", 
                         ylab="", 
                         xaxis=T, 
                         yaxis=T, 
                         main="boxplot", 
                         names=NULL, 
                         las=1,
                         cex.names=1, 
                         grid=T, 
                         grid.at=NULL,
                         grid.lwd=1,
                         grid.lty='dotted',
                         vio.wex=1.2, 
                         show.numb=c('none', 'median', 'mean', 'median.top', 'median.bottom', 'numb.top', 'numb.bottom'), 
                         numb.cex=.6, 
                         numb.col='black', 
                         numb.pos=3, 
                         ...){
  
  p_load(vioplot)
  
  show.numb <- match.arg(show.numb)
  
  ####################################
  ## names for x-axis
  if(is.null(names) & is.null(names(x)))
    names.x=1:length(x)
  else if(!is.null(names(x)))
    names.x=names(x)
  else
    names.x=names
  
  
  ##################################
  # remove NAN/Inf
  x <- lapply(x, function(x) {x=na.omit(x);x=x[is.finite(x)]})
  
  # ylim
  if(is.null(ylim))
    ylim=range(unlist(x))
  
  ##################################
  if(length(x) == 1){
    x = unlist(x)
    plot(NA, axes=F, xlab=xlab, ylab=ylab, ylim=ylim, type="n", main=main, ...)
    
    try(vioplot( x, add=T, at=1, col=col, drawRect=drawRect, wex=vio.wex, ...))
    boxplot( x, add=T, at=1, border=box.border, pch=box.pch, ...)
  }
  
  #################################
  if(length(x) > 1){
    
    ## xlim
    if(is.null(xlim))
      xlim=c(0.5, length(x)+0.5)
    
    ## at
    if(is.null(at))
      at=1:length(x)
    if( !is.null(at) ){
      xlim=c(0, max(at)+1)
    }
    
    ## colour
    if(length(col) == 1)
      col=rep(col, length(x) )
    if(length(box.border) == 1)  
      box.border <- rep(box.border, length(x))
    
    ## initialize plot
    plot(NA, axes=F, ylim=ylim, xlim=xlim, type="n", xlab=xlab, ylab=ylab, main=main, ...)
    
    ## axes
    if(xaxis){
      axis(1, at=at, labels=names.x, las=las, cex.axis=cex.names)
    }
    if(yaxis){
      axis(2, las=2)
    }
    ## grid
    if(grid){
      if(is.null(grid.at))
        abline( h=floor(ylim[1]):ceiling(ylim[2]), col='grey', lty=grid.lty, lwd=grid.lwd )
      else
        abline( h=grid.at, col='grey', lty=grid.lty, lwd=grid.lwd )
      
    }
    
    ## plot each box/violin
    for(i in 1:length(x)){
      
      if(length(x[[i]]) > 0){
        
        try(vioplot( x[[i]], add=T, at=at[i], col=my.col2rgb(col[i], vio.alpha), drawRect=drawRect, border=F, wex=vio.wex,...))
        boxplot( x[[i]], add=T, at=at[i], border=box.border[i], pch=box.pch, col=col[i], axes=F,...)
        
        if(show.numb=='median')
          text(at[i], median(x[[i]]), round(median(x[[i]]),2), pos=numb.pos, cex=numb.cex, offset=0.1, col=numb.col )
        if(show.numb=='mean')
          text(at[i], mean(x[[i]]), round(median(x[[i]]),2), pos=numb.pos, cex=numb.cex, offset=0.1, col=numb.col )
        if(show.numb=='median.top')
          text(at[i], ylim[2], round(median(x[[i]]),2), pos=1, cex=numb.cex, offset=0.1, col=numb.col )
        if(show.numb=='median.bottom')
          text(at[i], ylim[1], round(median(x[[i]]),2), pos=3, cex=numb.cex, offset=0.1, col=numb.col )
        if(show.numb=='numb.top')
          text(at[i], ylim[2], length(x[[i]]), pos=1, cex=numb.cex, offset=0.1, col=numb.col )
        if(show.numb=='numb.bottom')
          text(at[i], ylim[1], length(x[[i]]), pos=3, cex=numb.cex, offset=0.1, col=numb.col )
        
        ##text(at[i], mean(x[[i]]), "*", adj=c(0,0) )
      }
      
    }
    
  }
  
}

###################################################################
#                  Bland-Altman-Plot
#
# 20140122 implementation
# 20140214 added +/- sd in MA plot
###################################################################
ba.plot <- function( x, y, type=c("BA", "MA"), ylim=NULL, xlim=NULL, main="", ... ){
  
  type=match.arg(type)
  
  ######################
  ## Bland-Altman
  if(type == "BA"){
    Y <- x - y
    X <- (x+y)/2
    
        #fancyPlot(X, Y, ...)
        plot(X, Y, ylim=ylim, xlim=xlim,...)
    }

    ######################
    ## MA plot
    if(type == "MA"){

        M <- log(x/y, 2)
        A <- .5*log(x*y, 2)

        rm.idx <- union( which(is.infinite(M)), which(is.na(M)) )
        rm.idx <- union( rm.idx, which( is.infinite(A)) )
        rm.idx <- union( rm.idx, which( is.na(A)) )

        if(length(rm.idx)>0){
            M <- M[-rm.idx]
            A <- A[-rm.idx]
        }

        if(is.null(ylim))
            ylim=c(-max(abs(M[!is.infinite(M)]), na.rm=T), max( abs(M[!is.infinite(M)]), na.rm=T))

        fancyPlot(A, M, xlim=xlim, ylim=ylim, xlab=expression(A~~(~0.5~log[2](x~y))), ylab=expression(M~~(log[2](x/y))), cor=F, boxplots="y", main=paste("MA plot", main, sep="\n"), ...)
        abline( h=0, lwd=2, col="darkblue", lty="dashed")
        legend("bottomright", legend=c(paste("Median:", round(median(M, na.rm=T), 3)), paste("IQR:", round(IQR(M, na.rm=T), 3)) ) )

        abline(h=mean(M, na.rm=T), lwd=2, lty="dashed", col="darkred")
        
        upper.sd <- mean(M, na.rm=T)+1.96*sd(M, na.rm=T)
        lower.sd <- mean(M, na.rm=T)-1.96*sd(M, na.rm=T)
        
        abline(h=upper.sd, lwd=2, lty="dashed", col="darkred")
        abline(h=lower.sd, lwd=2, lty="dashed", col="darkred")
        legend("bottomleft", legend=c(paste("+/- 1.96 StdDev")), col="darkred", lty="dashed", lwd=2)
        
        out <- list()
        out[[1]] <- names(M)[which(M > upper.sd)]
        out[[2]] <- names(M)[which(M < lower.sd)]
        names(out) <- c('upper.bound', 'lower.bound')
    
        return(out)    

    }


}


################################################################################################
#                                ratio vs intensity plot
#
# rat            - numeric, ratios on raw scale
# int            - numeric, intensities on raw scale
# p              - numeric, significance values
#
# table          - NULL or a MaxQuant table. currently only the protein groups table is supported
#                - if specified, the rat/int/p vectors are not used but the respective
#                  columns are extracted from that table
# experiment     - character, name of an experiment                                  ## only considered if 'table' != NULL
# silac.state    - character, specifies which SILAC ratio to plot                    ## only considered if 'table' != NULL
# norm           - logical, if TRUE the normalized ratios are used                   ## only considered if 'table' != NULL
# signif         - character, specifies which signifiance values are used            ## only considered if 'table' != NULL
#
# label          - character vector of same length as rat/int/p
# label.string   - character, the string is used to search 'label' via grep
# pcut           - numeric, cutoff of p
# adjust         - character, method of p-value adjustment, "none" means no correction
# col            - charater, color used for unsignificant ratios
# sig.col        - character, color used for significant ratios
# label.col      - character, color used for the label
# pch            - numeric vector, specifying the plotting symbols for unregulated and regulated features
# label.ph       - numeric, plotting symbol for the label
# alpha          - numeric, transparency of points
# boxplots       - character, see function scatterplot in the car package
# legend.left    - logical, if TRUE a legend in the top left corner is shown
# legend.right   - logical, if TRUE a legend in the top right corner is shown
# xlim           - NULL or numeric vector of length 2
# reset.par      - logical, see function scatterplot in the car package
# label.legend   - logical, if TRUE a legend at the top will be plotted depicting the
#                  marked labels
#
# value:
#      list containing the indices of up and downregulated ratios
#
# changelog: 20100219 implementation
#            20100223 some documentation
#            20100311 added transparancy to the points -> 'alpha'
#                     switched to function 'scatterplot' from the 'car' package
#            20100610 added parameter 'table', 'experiment', 'silac.state', 'norm', 'signif'
#                     the idea was to replace the old 'sigRatioPlot' function that was used by
#                     the 'summaryPDF' function to plot quantified protein groups
#            20100909 several different labels can be plotted in different colors
#            20100915 if 'norm=F'meaning the unnormalized ratios are used, no regulated
#                     protein groups are reported since the significance values are based
#                     based on normalized ratios
#            20101021 compatibility for MQ v1.1.x:
#                      - if there are no significance column in the table, all p-values are set
#                        to one
#            20110505 if no valid values are found a dummy plot will be created and the functions stops.
#            20110701 if no significance column can be found the significance B values are
#                     calculated using my implementation
#            20120705 - tickmarks at x-axis based on xlim values
#                                        #            20121121 -  'ylim' as parameter
#
#################################################################################################
sigRatioPlot <- function(rat, int, p=NULL, signif=c("B", "A"), label=NULL, pcut = 0.01, adjust=c("none", "BH"), col="black", sig.col="red", label.col="darkred", label.pch=17, label.all=F, alpha=100, pch=c(16, 16), boxplots=c("xy", "x", "y", "n"), legend.left=T, legend.right=T, xlim=NULL, ylim=NULL, reset.par=TRUE, label.legend=F, main="", ...){

    # original parameter
    opar <- par(no.readonly=TRUE)

    signif=match.arg(signif)


    # make it all numeric
    if(is.factor(rat))
        rat <- as.numeric(as.character(rat))
    if(is.factor(int))
        int <- as.numeric(as.character(int))

    #########################
    ## estimate significance
    if(is.null(p)){
        if(signif =='A')
            p <- significanceA(rat, log=ifelse(sum(rat<0 , na.rm=T)>0, F, T))
        if(signif =='B')
            p <- significanceB(rat, int, log=ifelse(sum(rat<0 , na.rm=T)>0, F, T))[[1]]
    }
    if(is.factor(p))
        p <- as.numeric(as.character(p))

    p_load(car)

    ##################################
    # store the orignial index of
    # the data values before removing
    # missing values
    ##################################
    data.index <- 1:length(rat)
    names(data.index) <- 1:length(rat)

    ##################################
    # remove missing values
    ##################################
    rat.na <- which( is.na(rat) )
    int.na <- which( is.na(int))
    p.na <- which( is.na(p) )
    na.idx <- union(  union(rat.na, int.na), p.na )

   if(length(na.idx)>0){
       rat <- rat[-na.idx]
       int <- int[-na.idx]
       p <- p[-na.idx]
       data.index <- data.index[-na.idx]

       if(!is.null(label))
           label <- label[ -na.idx ]
   }


    ###################################
    # boxplots on axes
    ###################################
    boxplots <- match.arg(boxplots)

    ###################################
    # log
    ###################################
    rat.log <- log(rat, 2)
    int.log <- log(int, 10)
    int.log[is.infinite(int.log)] <- NA




    #########################################################################################
    #
    #          - if there are no valid values, stop right here.....
    #          - produce an empty plot for compatibility with PDF script
    #
    #########################################################################################
    if(length(intersect(which(!is.na(rat)),which(!is.na(int)))) < 1  ){

        plot(1,1, col='white', axes=F, xlab='', ylab='')

        return(1)
    }


    ###################################
    # number of quantified items
    ###################################
    quant.numb <- sum(!is.na(rat.log), na.rm=T)

    ###################################
    # plotting symbols
    ###################################
    pch.vec <- pch

    ###################################
    #   adjust p-values
    ###################################
    adjust <- match.arg(adjust)
    p <- p.adjust(p, adjust)

    ###################################
    # determine regulated
    ###################################
    rat.log.reg <- rat.log[ p <= pcut]
    int.log.reg <- int.log[ p <= pcut]

    # index of regulated
    regulated <- ifelse(p <= pcut, 1, 0)

    # in numbers
    down.numb <- sum(rat.log.reg < 0, na.rm=T)
    up.numb <- sum(rat.log.reg > 0, na.rm=T)

    # index of regulated
    down.idx <- which((rat.log < 0) & (p <= pcut))
    up.idx <- which((rat.log > 0) & (p <= pcut))

    ###################################
    # colors
    ###################################
    col.org <- col
    sig.col.org <- sig.col

    col <- my.col2rgb(col, alpha)
    sig.col <- my.col2rgb(sig.col, alpha)


   ## if(length(label.col) < length(label.string) )
   ##     label.col <- rep(label.col[1], length(label.string))
    for(cc in 1:length(label.col)){

        label.col[cc] <- my.col2rgb(label.col[cc], alpha)
    }

    # set the color palette: the first entry is just a dummy -> used in function 'scatterplot'
    #palette( c("black", col, sig.col, label.col[as.numeric(levels(regulated)[3:length(levels(regulated))]) - 1]))     # R version 2.10.1
    palette( c(col, sig.col, label.col[as.numeric(levels(regulated)[3:length(levels(regulated))]) - 1]))               # new R version

    ###################################
    # xlim  &  ylim
    ###################################

    if(is.null(xlim)){
        xmax <- max(abs(rat.log), na.rm=T)

        xlim <- c(-1*xmax, xmax)
    }

    # ylim
    if(is.null(ylim))
    {
        ylim <- c(min(int.log, na.rm=T),max(int.log, na.rm=T)+1)
    }

    ##########################################
    # original plotting parameters
    ##########################################
    #par.org <- par(no.readonly = TRUE)

    ##########################################
    # plot
    ##########################################
    scatterplot( rat.log ,int.log, xlim=xlim, ylim=ylim, pch=pch.vec, xlab=expression(log[2](Ratio)), ylab=expression(log[10](Intensity)), boxplots=boxplots, reg.line=F, smooth=F, groups=factor(regulated), legend.plot=F, main=main, axes=F, reset.par=F, ...)
    axis(1, at=seq( ceiling(xlim[1]), floor(xlim[2]), 1 )  )
    axis(2)
    # legends
    if(legend.left){
        legend("topleft", legend=c(paste("down:", down.numb), paste("up:", up.numb), paste("quantified:", quant.numb) ), bty="n", inset=c(0.02,0))
        ##else legend("topleft", legend=c(paste("quantified:", quant.numb) ), bty="n", inset=c(0.02,0))
    }
    if(legend.right)
        legend("topright", legend=c( paste("p adjustment:", adjust), paste("p cutoff:", pcut)), bty="n"  )

    if(label.legend)
        legend("top", legend=label.string, col=label.col, pch=label.pch, bty="n"  )


    out <- vector("list", 2)
    names(out) <- c("down", "up")


    out[[1]] <- data.index[down.idx]
    out[[2]] <- data.index[up.idx]


    # default color palette
    palette("default")

    # default plotting parameters
    if(reset.par)
        par(opar)


    return(out)
}



###############################################
## Plot a 2-Way, 3-Way or 4-Way Venn Diagram ##
###############################################
## Author: Thomas Girke
## Last update: Nov 6, 2008
## Utility: Plots a non-proportional 2-, 3- or 4-way venn diagram based on overlaps among data sets (vectors)
## Detailed instructions for running this script are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn

## Define venndiagram function
venndiagram <- function(x=x, y=y, z=z, w=w, unique=T, title="Venn Diagram", labels=c("x", "y", "z", "w"), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE, ...) {
	## Remove duplicates and NA fields in x, y, z and w
	if(unique==T) {
		x <- unique(x); x <- as.vector(na.omit(x))
		y <- unique(y); y <- as.vector(na.omit(y))
		if(!missing("z")) {
			z <- unique(z); z <- as.vector(na.omit(z))
		}
		if(!missing("w")) {
			w <- unique(w); w <- as.vector(na.omit(w))
		}
	}

	## Check valid type selection
	if(!type %in% c("2", "2map", "3", "3map", "4", "4map", "4el", "4elmap")) {
		return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map, 4el, 4elmap.")
	}

	## Plot a 2-way venn diagram
	if(type=="2") {
		# Define ovelap queries
		q1 <- x[x %in% y]
		q2 <- x[!x %in% y]
		q3 <- y[!y %in% x]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=countDF$count)
		if(printsub==TRUE) {mysub <- paste(paste("N unique: xy =", length(unique(c(x,y)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), sep="")} else {mysub <- ""}
		if(plot==T) {
			## Plot the 2-way venn diagram
			symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), labels[1:2], col=lcol, ...)
		}

		## Return query list
		return(qlist)
	}

	## Plot 2-way mapping venn diagram
	if(type=="2map") {
		olDFdebug <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=paste("q", 1:3, sep=""), ...)
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), paste(labels[1:2], "=", c("x","y")), col=lcol, ...)
	}

	## Plot a 3-way venn diagram
	if(type=="3") {
		## Define ovelap queries
		q1 <- x[x %in% y & x %in% z]
		q2 <- x[x %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- y[y %in% z]; q3 <- q3[!q3 %in% x]
		q4 <- x[x %in% y]; q4 <- q4[!q4 %in% z]
		q5 <- x[!x %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- y[!y %in% z]; q6 <- q6[!q6 %in% x]
		q7 <- z[!z %in% x]; q7 <- q7[!q7 %in% y]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=countDF$count)
		if(printsub==TRUE) {mysub <- paste(paste("N unique: xyz =", length(unique(c(x,y,z)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), sep="")} else { mysub <- "" }
		if(plot==T) {
			## Plot the 3-way venn diagram
			symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels[1:3], col=lcol, ...)
		}

		## Return query list
		return(qlist)
	}

	## Plot 3-way mapping venn diagram
	if(type=="3map") {
		olDFdebug <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=paste("q", 1:7, sep=""), ...)
		symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), paste(labels[1:3], "=", c("x","y","z")), col=lcol, ...)
	}

	## Overlap queries for 4-way venn diagram
	if(type=="4" | type=="4el" | type=="4elmap") {
		## Define ovelap queries
		xy <- x[x %in% y]; xz <-x[x %in% z]; xw <- x[x %in% w]; yz <- y[y %in% z]; yw <- y[y %in% w]; zw <- z[z %in% w]
		q1 <- xy[xy %in% zw]
		q2 <- xw[xw %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- yz[yz %in% w]; q3 <- q3[!q3 %in% x]
		q4 <- yz[yz %in% x]; q4 <- q4[!q4 %in% w]
		q5 <- xw[xw %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- xy[!xy %in% z]; q6 <- q6[!q6 %in% w]
		q7 <- zw[!zw %in% x]; q7 <- q7[!q7 %in% y]
		q8 <- xz[!xz %in% y]; q8 <- q8[!q8 %in% w]
		q9 <- yw[!yw %in% x]; q9 <- q9[!q9 %in% z]
		q10 <- x[!x %in% c(y,z,w)]
		q11 <- y[!y %in% c(x,z,w)]
		q12 <- z[!z %in% c(x,y,w)]
		q13 <- w[!w %in% c(x,y,z)]
		q14 <- xw[!xw %in% y]; q14 <- q14[!q14 %in% z]
		q15 <- yz[!yz %in% x]; q15 <- q15[!q15 %in% w]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, q11=q11, q12=q12, q13=q13, q14=q14, q15=q15)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=countDF$count)

		if(printsub==TRUE) {mysub <- paste(paste("N unique: xyzw =", length(unique(c(x,y,z,w)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), paste("; w =", length(unique(w))), sep="") } else { mysub <- "" }

	## Plot 4-way venn diagram as circles
		if(plot==T & type=="4") {
			symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], col=tcol, ...) # rows 14-15 of olDF are printed in last step
			text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels, col=lcol, ...)
			text(c(3.8, 3.8), c(1.0, 0.4), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDF$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDF$count[15], sep="")), col=diacol, ...)
		}

	## Plot 4-way venn diagram as ellipses
	if(plot==T & (type=="4el" | type=="4elmap")) {
		olDF <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=countDF$count)
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 4-way venn diagram
		ellipseVenn <- function(lines=lines, olDF, title=title, labels=labels, sub=mysub, main, lcol=lcol, tcex=1.3, ...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=title, sub=mysub, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, ...)
			text(olDF[1:15,1], olDF[1:15,2], olDF[1:15,3], col=tcol, ...)
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels, col=lcol, ...)
			close.screen(all=TRUE)
		}
		## Plot 4-way ellipse venn diagram
		if(type=="4el") {
			ellipseVenn(olDF=olDF, lcol=lcol, lines=lines, labels=labels, title=title, ...)
		}

		## Plot 4-way ellipse mapping venn diagram
		if(type=="4elmap") {
			olDFdebug <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=paste("q", 1:15, sep=""), ...)
			ellipseVenn(olDF=olDFdebug, lcol=lcol, lines=lines, labels=paste(labels, "=", c("x","y","z","w")), title="Mapping Venn Diagram", ...)
		}
	}

	## Return query list
	return(qlist)
	}

	## Plot 4-way circle mapping venn diagram
	if(type=="4map") {
		olDFdebug <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=paste("q", 1:15, sep=""), ...)
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x[1:13], olDFdebug$y[1:13], olDFdebug$count[1:13], col=tcol, ...); text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), paste(labels, "=", c("x","y","z","w")), col=lcol, ...)
		text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDFdebug$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDFdebug$count[15], sep="")), col=tcol, ...)
	}

}


######################################################################
##
## http://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
## downaloaded: 20140303
######################################################################
my.plotcorr <- function (corr, outline = FALSE, col = "grey", upper.panel = c("ellipse", "number", "none"), lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"), digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...)
{
# this is a modified version of the plotcorr function from the ellipse package
# this prints numbers and ellipses on the same plot but upper.panel and lower.panel changes what is displayed
# diag now specifies what to put in the diagonal (numbers, ellipses, nothing)
# digits specifies the number of digits after the . to round to
# unlike the original, this function will always print x_i by x_i correlation rather than being able to drop it
# modified by Esteban Buz
  if (!p_load('ellipse', quietly = TRUE, character = TRUE)) {
    stop("Need the ellipse library")
  }
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ #diag behavior
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ #lower half of plot
      if (lower.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { #upper half of plot
      if (upper.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}



#######################################################################
#                        colorGradient
#
# arguments
#   x          numeric, value between 'mi' and 'ma'
#   scheme     color scheme to use
#
# value
#   RGB color corresponding to x
#
#
#######################################################################
colorGradient <- function(x, mi=0, ma=1, scheme=c("heat", "topo", "terrain", "grey", "cm"), eps=0.02 )
{
    # intervall
    int <- seq(mi, ma, length.out=1000  )

    # color scheme
    if(length(scheme) == 1){

      scheme=match.arg(scheme)

      if(scheme=="heat")
         cols <- heat.colors(length(int))
      if(scheme=="topo")
        cols <- topo.colors(length(int))
      if(scheme=="terrain")
        cols <- terrain.colors(length(int))
      if(scheme=="grey")
        cols <- grey.colors(length(int))
      if(scheme=="cm")
        cols <- cm.colors(length(int))
    }
    if(length(scheme) == 3){
        p_load(gplots)
        cols <- colorpanel(length(int), scheme[1], scheme[2], scheme[3])
    }
    if(length(scheme) > 3){

        int <- seq(mi, ma, length.out = length(scheme))
        cols <- scheme

    }

    return( cols[ max( which( x >= (int-eps) ) )  ] )
}



####################################################################################################################
##
##
## parameter
##    data - feature matrix without NAs, raw scale
##    norm - either 'NULL' or a character specifying the normalization method to use
##           transformation:
##               Zscore      - transform to Z-score row-wise
##           distance measures:
##               euclidean   - usual square distance bewteen two rows ( 2 norm )
##               maximum     - maximum distance between two components of the rows ( supremum norm )
##               canberra
##           correlation:
##               spearman    - Spearman's rank correlation coefficient
##               pearson     - pearson correlation coefficient
##    log       - logical indicating whether the data should be transformed to log-scale
##    logBase   - numerical, the base of logarithm
##    plotOnRawScale   - logical, if TRUE the y-axis of the time profile plot is the non-normalized scale of 'data'
##    ylim
##    algorithm
##    nClust
##    file
##    pdf
##    nrows     - numeric, number of rows in the plot
##    ...
## value
##   the object returned by
##   the cluster algorithm
##
##
## changelog:
##             20100721 added parameter 'hc'; can be an 'hclust' object
##                      Zscore transformation is done by function 'scale'
##             20120320 parameter 'alpha': for non-colorgradients, e.g. hclust, kmeans, ...
##             20120221 parameter 'legend': if TRUE a legend will be plotted, at least in case
##                      of 'cmeans' and 'conensusclust'
##             20130528 changed color of averaged profile to darkred
##             20140123
##
#############################################################################################################
timeProfileClusterPlot <- function(data, hc=NULL, norm=c("none", "Zscore"), log=T, logBase=2, plotOnRawScale=F, addAverageProfile=T, y.lim=NULL, algorithm=c("cmeans", "mycmeans", "kmeans", "consensus", "hclust"), nClust=3, file="clusterTimeProfiles.pdf", pdf=F, ylab="Change in expression", col=c("green", "magenta", "darkred"), nrows=2, minMembShip=0.5, cmeansm=2, alpha=40, legend=T,...)
{
    ##alpha=80

    if(!is.null(hc) && class(hc) == "hclust")
        algorithm="hclust"
    if(!is.null(hc) && class(hc) == "fclust")
        algorithm="mycmeans"

    ####################################
    #  numeric matrix
    ####################################
    if(mode(data) != "numeric")
        data <- data.matrix(data)

    ####################################
    #      log transformation
    ####################################
    if(log)
        data <- log(data, logBase)
    # data before normalization
    dataRaw <- data

    #####################################
    #        normalize the data
    #####################################
    norm=match.arg(norm)

    if(norm != "none")
    {
       ################################
       #    Z-score
       ################################
       if(norm == "Zscore")
           data <- t( scale(t(data), center=T, scale=T))
       else if(norm == "pearson" | norm == "spearman")
           data <- cor(t(data), method=norm)
       else{
           data <- dist(data, method=norm)
       }
    }
    ######################################################
    #
    #  different colors for each cluster
    #
    ######################################################
    if(nClust == length(col)){
           COLORS=col
    } else {
           COLORS=rep("black", nClust)
    }
    for(ii in 1:length(COLORS))
           COLORS[ii] <- my.col2rgb( as.character(COLORS[ii], alpha=alpha)  )

    #######################################################################################
    #
    #                            perform the clustering
    #
    #######################################################################################
    if(is.null(hc)){

       # determine cluster algorithm
       algorithm=match.arg(algorithm)

       ################################
       # kmeans consensus clustering
       ################################
       if(algorithm == "consensus")
       {
           p_load(compdiagTools)
           clustering <- consensusCluster( t(data),  nclass=nClust, ...)
       }
       ################################
       # kmeans consensus clustering
       ################################
       if(algorithm == "consensus+") {
           p_load(ConsensusClusterPlus)
       }
       #################################
       #
       #################################
       if(algorithm == "mycmeans"){

           clustering <- my.cmeans(data, ...)
       }
       #################################
       #         cmeans
       #################################
       if(algorithm == "cmeans"){

           p_load(e1071)
           #clustering <- eval(parse( text=paste(algorithm, "(data, centers=", nClust,", ...)", sep=""  )  ))

           N <- dim(data)[1]
           D <- dim(data)[2]

           cmeansm <-  1 + ( (1418/N) + 22.05) * D^-2 + ( (12.33/N) + 0.243 )*D^(-0.0406*log(N, exp(1))-0.1134)

           cat("optimal fuzzyfication parameter 'm'=", cmeansm, "\n")


           clustering = cmeans(data, centers=nClust, m=cmeansm, ...)

           ##

       }
       #################################
       #        k-means
       #################################
       if(algorithm == "kmeans"){
           clustering = kmeans(data, centers=nClust, ...)
       }

       ########################################
       # cluster assigment
       ########################################
       cl <- clustering$cluster

    } else{ # end if is.null(hc)

       ##################################################
       #         hierarchical clustering
       ##################################################
       if(class(hc) == "hclust"){
           cl <- cutree(hc, nClust)
       }
       if(class(hc) == "fclust"){
           cl <- hc$cluster
           clustering=hc
       }

    }

    ##########################################################
    #
    #                   plot the stuff
    #
    ##########################################################

    # determine the number of clusters
    nC <- length(levels(as.factor(cl)))

    # determine whether the normalized data shall be plotted
    if( plotOnRawScale )
        dataPlot <- dataRaw
    else
        dataPlot <- data

    #################################
    # ylim, equal for all clusters
    #################################
    if(is.null(y.lim)){
            y.lim = rep(max(abs(dataPlot), na.rm=T), 2)
            y.lim[1] <- -1*y.lim[1]
        }



    #################################
    # open pdf file
    #################################
    if(pdf) pdf(file=file, width=ceiling(nC/nrows)*6, height= 5*nrows)

    par(mfrow=c(nrows, ceiling(nC/nrows)))

    #################################
    #        loop over the cluster
    #################################
    ccount=1
    for( cluster in levels(as.factor(cl)) ){

        ## determine the members of the current cluster
        clusterMembers = rownames(dataPlot)[ which( cl == cluster ) ]


        ################################################
        #
        ################################################
        #COLOR <- my.col2rgb( "black", alpha=alpha )
        COLOR = COLORS[ccount]

        ################################################
        # color gradient respective to membership
        ################################################
        if(algorithm=="cmeans" | algorithm=="mycmeans"){

            ## sort according to cluster mebership
            cl.ord <- order(clustering$membership[ clustering$cluster == cluster, cluster])
            clusterMembers <- clusterMembers[cl.ord]

            item.score.range <- c(0,1)
            ## item.score.range = range()
            COLOR <- colorGradient( max(clustering$membership[ clusterMembers[1], ]) , mi=item.score.range[1], ma=item.score.range[2], scheme=col )
            #COLOR <- colorGradient( max(clustering$membership[ clusterMembers[cl.ord[1]], ]) , mi=item.score.range[1], ma=item.score.range[2], scheme=col )

            COLOR <- my.col2rgb( COLOR, alpha=alpha )

            ## calculate the mean cluster membership score
            maxScore <- apply( clustering$membership[clusterMembers, ],1, max  )

        }

        ################################################
        # color gradient respective to membership
        ################################################
        if(algorithm=="consensus")
        {

            # min and max item consensus
            #item.score.range <- range(clustering$itemConsensus[ clusterMembers ])
            item.score.range <- range(clustering$itemConsensus)
            #item.score.range <- c(0,1)
            COLOR <- colorGradient( clustering$itemConsensus[ clusterMembers[1] ] , mi=item.score.range[1], ma=item.score.range[2], scheme=col )

            COLOR <- my.col2rgb( COLOR, alpha=clustering$itemConsensus[ clusterMembers[1] ] * alpha   )


        }


        # plot the first profile....
        plot(dataPlot[clusterMembers[1], ], type="l", ylim=y.lim, main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , col=COLOR, xaxt="n", ylab=ylab, lwd=2, xlab="" )

        if((algorithm=="cmeans" | algorithm=="mycmeans") & legend)
            legend("topright", legend=c( paste( "mean: ", round(mean(maxScore),2) ) , paste("median: ", round(median(maxScore),2)), paste("N <", minMembShip, ":", sum(maxScore < minMembShip)) ), title="Membership" )

        if(algorithm=="consensus" & legend)
            legend("topright", legend=c( paste( "cluster consensus: ", round(clustering$clusterConsensus[as.numeric( cluster) ], 3) ) ) )

        axis(1, at=seq(1,dim(data)[2]), labels=colnames(data), las=2, cex.axis=0.7)


        # if there are other cluster memebrs....
        if( length(clusterMembers) > 1 ) {
        # ... and loop over the remaining members
           for( i in 2:length(clusterMembers)){

                if(algorithm=="cmeans" | algorithm=="mycmeans"){
                   COLOR <- colorGradient( max(clustering$membership[clusterMembers[i],]), mi=item.score.range[1], ma=item.score.range[2], scheme=col )
                      COLOR <- my.col2rgb( COLOR, alpha=alpha )
                }
                if(algorithm=="consensus"){
                   COLOR <- colorGradient( clustering$itemConsensus[clusterMembers[i]], mi=item.score.range[1], ma=item.score.range[2], scheme=col )
                   COLOR <- my.col2rgb( COLOR, alpha= clustering$itemConsensus[clusterMembers[i]] * alpha )
                 }
                lines(dataPlot[clusterMembers[i], ], type="l", main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , col=COLOR, lwd=2)

              }

           # add the average profile
           if(addAverageProfile)
                  lines( apply(dataPlot[clusterMembers, ], 2, mean, na.rm=T), lwd=2, col="darkred", lty="solid" )
        }
        ccount = ccount + 1
    }
    par(mfrow=c(1,1))

    # close pdf file
    if(pdf)dev.off()

    if(is.null(hc)){
        return(clustering)
    } else {
        return(cl)
    }
}
####################################################################################################################
#
#
# changelog: 20111215 implementation
####################################################################################################################
my.cmeans <- function(x, centers=NULL, iter.max=1000, maxK=NULL, verbose=FALSE, dist="euclidean", method="cmeans", m=NULL, rate.par=NULL, weights=1, control=list(reltol=(.Machine$double.eps))){


    # x <- rbind( matrix(rnorm(100,sd=0.3),ncol=2), matrix(rnorm(100,sd=0.3, mean=3),ncol=2), matrix(rnorm(100,sd=0.3, mean=10), ncol=2), matrix(rnorm(100,sd=0.3, mean=5),ncol=2) )


    p_load(e1071)

    #########################################################
    # estimate optimal value for fuzzyfication parameter 'm'
    # Schw?mmle et.al 2010
    #########################################################
    if(is.null(m)){

        D <- dim(x)[2]
        N <- dim(x)[1]

        mt <- 1 + ( (1418/N) + 22.05) * D^-2 + ( (12.33/N) + 0.243 )*D^(-0.0406*log(N, exp(1))-0.1134)

        cat("optimal fuzzyfication parameter 'm' estimated as:", mt, "\n")
    }
    #########################################################
    # determine the number of clusters
    #########################################################
    if(is.null(centers)){

        # maximal cluster number to test
        if(is.null(maxK))
            K = floor( sqrt(dim(x)[1]) )
        else
            K = maxK

        # loop over different number of clusters
        min.centroid.dist <- vector("numeric", K-1)
        names(min.centroid.dist) <- 2:K

        clustering <- vector("list", K-1)
        names(clustering) <- 2:K

        for(cc in 2:K){

            # cluster the data
            cl.tmp <- cmeans(x, centers=cc, iter.max=iter.max, verbose=verbose, dist=dist, method=method, m=mt, rate.par=rate.par, weights=weights, control=control)
            clustering[[as.character(cc)]] <- cl.tmp

            # get centroids
            centroids.tmp <- cl.tmp$centers
            rownames(centroids.tmp) <- paste( "center_", rownames(centroids.tmp), sep=""  )

            # estimate distances to the centroids
            #dist.centroids <- as.matrix( dist(rbind(x, centroids.tmp), method=dist ))
            dist.centroids <- as.matrix( dist(centroids.tmp, method=dist ))
            diag(dist.centroids) <- NA

            # minimal centroid distance
            min.centroid.dist[as.character(cc)] <- min( dist.centroids^2, na.rm=T )
        }
        ####################################
        # calculate K_opt
        ####################################
        dist.decay <- vector("numeric", length(min.centroid.dist)-1)

        for(i in 1:(length(min.centroid.dist)-1))
            dist.decay[i] <- abs(min.centroid.dist[i] - min.centroid.dist[i+1])


        X11()
        # par(mfrow=c(1,2))
        plot(2:K, min.centroid.dist, main="Minimum centroid distance", xlab="cluster number", type="b", pch=20)
        #barplot(dist.decay)

        #
        dist.decay.norm <- dist.decay/max(dist.decay)
        #copt <- max(which(dist.decay.norm > 0.05)) + 1

        copt <- (dist.decay.norm > 0.05)
        copt[min(which(copt == FALSE)):length(copt) ] <- FALSE
        copt <- max(which(copt == TRUE)) + 1


        cat("optimal number of clusters:", copt, "\n")


    }
    return(clustering[[as.character(copt)]])
    #return(copt)
    #return(dist.centroids)
    #return(dist.decay.norm)
    #return(min.centroid.dist)

}


##For the code to be effective in R v1.6.x this function has to replace plot.histogram() in the base environment, i.e.
##source("plot.histogram.R")
##assign("plot.histogram", plot.histogram, pos=which(search()=="package:base"))
##Code

#########################################################################/**
# \name{plot.histogram}
# \alias{plot.histogram}
#
# \title{Plots a histogram}
#
# \usage{
#   plot.histogram(x, freq=equidist, col=NULL, border=par("fg"), lty=NULL,
#      width=1.0, offset=(1.0-width)/2, main=paste("Histogram of", x$xname),
#      xlim=range(x$breaks), ylim=NULL, xlab=x$xname, ylab, axes=TRUE,
#      labels=FALSE, add=FALSE, ...)
# }
#
# \description{
#   This function redefines the \code{plot.histogram} function
#   in the \R base package by adding the two arguments \code{width} and
#   \code{offset}. The function is modified in such a way that it is
#   backward compatible, i.e. if you do not use the arguments \code{width}
#   and \code{offset} the plot will look the same as the plot generated by
#   the original function. Note that \code{plot.histogram} is called by
#   \code{hist}.
# }
#
# \arguments{
#   \item{x}{a `histogram' object, or a list with components
#       \code{intensities}, \code{mid}, etc, see \code{\link[base]{hist}}.}
#   \item{freq}{logical; if \code{TRUE}, the histogram graphic is to present
#       a representation of frequencies, i.e, \code{x$counts}; if
#       \code{FALSE}, relative frequencies ("probabilities"), i.e.,
#       \code{x$intensities}, are plotted. The default is true for
#       equidistant \code{breaks} and false otherwise.}
#   \item{col}{a colour to be used to fill the bars.  The default of
#       \code{NULL} yields unfilled bars.}
#   \item{border}{the color of the border around the bars.}
#   \item{width}{The relative width of each bar compared to the full width.
#       \code{1.0} is full width. Default value is \code{1.0}.}
#   \item{offset}{The relative horisontal offset of each bar compared to the
#       full width. A value of \code{0.0} places each bar to the very left.
#       A value of \code{1.0-width} places each bar to the very right.
#       Default value is \code{(1.0-offset)/2}, i.e. the bars are centered.}
#   \item{lty}{the line type used for the bars, see also \code{lines}.}
#   \item{xlim, ylim}{the range of x and y values with sensible defaults.}
#   \item{main, xlab, ylab}{these arguments to \code{title} have useful
#       defaults here.}
#   \item{axes}{logical, indicating if axes should be drawn.}
#   \item{labels}{logical or character.  Additionally draw labels on top of
#       bars, if not \code{FALSE}; if \code{TRUE}, draw the counts or
#       rounded intensities; if \code{labels} is a \code{character}, draw
#       itself.}
#   \item{add}{logical. If \code{TRUE}, only the bars are added to the
#       current plot. This is what \code{lines.histogram(*)} does.}
#   \item{...}{further graphical parameters to \code{title} and \code{axis}.}
# }
#
# \author{
#   Modified by Henrik Bengtsson,
#   \url{http://www.braju.com/R/}, from the original \R plot.histogram.
# }
#
# \examples{
#   x1 <- rnorm(1000,  0.4, 0.8)
#   x2 <- rnorm(1000,  0.0, 1.0)
#   x3 <- rnorm(1000, -1.0, 1.0)
#   hist(x1, width=0.33, offset=0.00, col="blue", xlim=c(-4,4),
#    main="Histogram of x1, x2 & x3", xlab="x1 - blue, x2 - red, x3 - green")
#   hist(x2, width=0.33, offset=0.33, col="red", add=TRUE)
#   hist(x3, width=0.33, offset=0.66, col="green", add=TRUE)
# }
#
# \seealso{
#   See also the original \code{\link[base]{hist}/plot.histogram} function
#   in the \R base package.
# }
#
# @visibility public
#*/#########################################################################

plot.histogram2 <-
function (x, freq = equidist, density = NULL, angle = 45, col = NULL,
    border = par("fg"), lty = NULL, main = paste("Histogram of",
        x$xname), xlim = range(x$breaks), ylim = NULL, xlab = x$xname,
    ylab, axes = TRUE, labels = FALSE, add = FALSE, width=1.0, offset=(1.0-width)/2, ...)
{
    equidist <- if (is.logical(x$equidist))
        x$equidist
    else {
        h <- diff(x$breaks)
        diff(range(h)) < 1e-07 * mean(h)
    }
    if (freq && !equidist)
        warning("the AREAS in the plot are wrong -- rather use `freq=FALSE'!")
    y <- if (freq)
        x$counts
    else x$intensities
    nB <- length(x$breaks)
    if (is.null(y) || 0 == nB)
        stop("`x' is wrongly structured")
    if (!add) {
        if (is.null(ylim))
            ylim <- range(y, 0)
        if (missing(ylab))
            ylab <- if (!freq)
                "Density"
            else "Frequency"
        plot.new()
        plot.window(xlim, ylim, "")
        title(main = main, xlab = xlab, ylab = ylab, ...)
        if (axes) {
            axis(1, ...)
            axis(2, ...)
        }
    }

    if (width != 1.0 || offset != 0) {
      # Calculates the width of each bar in the histogram
      delta.breaks <- x$breaks[-1] - x$breaks[-nB];
      x.offset <- offset * delta.breaks;
      x.width <- width * delta.breaks;
      x <- x$breaks[-nB]+x.offset;
      rect(x, 0, x+x.width, y, col=col, border=border, angle = angle, density = density, lty=lty);
    } else {
      rect(x$breaks[-nB], 0, x$breaks[-1], y, col = col, border = border,
          angle = angle, density = density, lty = lty)
    }

    if ((logl <- is.logical(labels) && labels) || is.character(labels))
        text(x$mids, y, labels = if (logl) {
            if (freq)
                x$counts
            else round(x$density, 3)
        }
        else labels, adj = c(0.5, -0.5))
    invisible()
} # plot.histogram



## ####################################################################
## 
ssGSEA.plot <-  function(data.expr, gs.name, samp.name, signat, score, pval){


            plot( data.expr, pch=20, col='darkgrey', lwd=2, type='l', xlab='Rank', ylab='Expression', main=paste(gs.name, samp.name, sep='\n'), ylim=range(data.expr), yaxs='i')


            ## #########################################################
            ##  ptm signatures?
            if(length(grep(';u$|;d$', signat[[gs.name]], value=T)) > 0){

                ## locations
                gsea.tmp.u <- sub(';u$','',grep(';u$', signat[[gs.name]], value=T))
                loc.u <- na.omit(match(gsea.tmp.u, gn))

                gsea.tmp.d <- sub(';d$','',grep(';d$',  signat[[gs.name]], value=T))
                loc.d <- na.omit(match(gsea.tmp.d, gn))

                if(!is.null(loc.u)){
                    ##rug(loc.u, col='red', side=3, lwd=2)
                    rug(loc.u, col=my.col2rgb('red',  50), side=3, lwd=2, ticksize=0.02)

                }
                if(!is.null(loc.d)){
                    rug(loc.d, col='green', lwd=2)
                }
                ## some info
                legend('bottom', legend=c(paste('No. down-regulated in signature:', length(grep(';d$', signat[[gs.name]]))),
                                          paste('No. found in data set:', length(loc.d))
                                          ), inset=.05, bty='n', text.col='darkgreen')

                legend('top', legend=c(paste('No. up-regulated in signature:', length(grep(';u$', signat[[gs.name]]))),
                                       paste('No. found in data set:', length(loc.u))
                                       ), inset=.05, bty='n', text.col='darkred')
            } else {## end if signature

                ## ####################################################
                ## regular gene set
                loc <- which(gn %in% signat[[gs.name]])
                rug(loc, col=my.col2rgb('red',  50), side=3, lwd=2, ticksize=0.02)

                ## box plot
                loc.quart <- quantile(loc)
                rug(loc.quart, col='darkblue', side=3, lwd=2, ticksize=0.03)
                rect( loc.quart[2], max(data.expr)-0.04*max(data.expr-min(data.expr)), loc.quart[4], max(data.expr), border='darkblue', lwd=2, col=NA )
                rect( loc.quart[1], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[2], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )
                rect( loc.quart[4], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[5], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )

                ## some info
                legend('bottom', legend=c(paste('No. in signature:', length( signat[[gs.name]])),
                                          paste('No. found in data set (non-redund.):', sum(signat[[gs.name]] %in% gn)),
                                          paste('No. found in data set (redundant):', length(loc))
                                       ), inset=.05, bty='n', text.col='darkred')

            }
            legend('center', legend=paste('NES=', round(score, 3), ' (p=', round(pval, 5), ')', sep=''), bty='n', inset=.2, cex=1.5)
}

## ############################################################
## 20180417
##         Generate heatmaps of significant gene sets
## - gene sets significant ('fdr.max') in at least 'fdr.rep' samples per category in 'class_vector'
##   
##
gene.set.hm <- function(gct.str,             ## path ssGSEA/ssPSEA GCT file (combined version)
                        class_vector='pert', ## must be present in 'cdesc' 
                        class_vector_levels=NULL, ## specify levels of interest 
                        hm_color_scheme=c('cyan2brown','blue2red'),
                        fdr.max=0.01,        ## max FDR; cells will be marked with 'cellnote_pch'  
                        fdr.rep=4,           ## min. number of replicates in 'class_vector'
                        regex=NULL,          ## regular expression applied to @rid
                        cellnote_pch='*',
                        cellnote_size=15,
                        cellnote_color='grey30',
                        cutree_cols=3,
                        anno_colors = list(pert=c('DMSO'='chocolate1', 'EGF'='chartreuse4', 'Nocodazole'='darkmagenta')),
                        main='',
                        max.val=NULL,
                        min.val=NULL,
                        filename=NA,
                        sort_col=F,        ## if TRUE columns will be sorted according to 'class vector' and NOT be clustered
                        cluster_rows=T,
                        cluster_cols=T,
                        exclude.from.cdesc=c('SM.id','Age', 'QC.status', 'QC.status.1', 'Race', 'Ethnicity', 'Ischemia.Time', 'normalization', 'Total.Cellularity', 'Tumor.Cellularity','Sample.Type', 'Channel', 'Sample.ID', 'Experiment', 'Necrosis'),
                        ...  ## further arguments passed to pheatmap
                        ){
  
  require(pacman)
  p_load(cmapR)
  p_load(magrittr)
  p_load(pheatmap)
  p_load(dichromat)
  p_load(RColorBrewer)
  
  hm_color_scheme <- match.arg(hm_color_scheme)
  
  # import
  gct <- parse.gctx(gct.str)
  
  ## extract scores
  score <- gct@mat %>% data.matrix
  rownames(score) <- gct@rid
  colnames(score) <- gct@cid
  
  ## regex
  if(!is.null(regex)){
    idx <- grep(regex, gct@rid)
    
    if(length(idx) ==  0)
      stop(paste('\n\nNo signatures found:', regex, '\n\n'))
    
    score <- score[idx, ]
    gct@rdesc <- gct@rdesc[idx, ]
    gct@mat <- gct@mat[idx, ]
    
    ## filename suffix
    fn.suffix <- gsub('\\^|\\$\\|', '', regex) %>% sub('(-|_)$','', .)
    
    if(!is.na(filename))
      filename <- sub('\\.pdf$', paste('_', fn.suffix, '.pdf', sep=''), filename)
    
  } 
  
  ## extract FDRs
  rdesc <- gct@rdesc
  fdr <-  rdesc[, grep('fdr.pvalue', colnames(rdesc))] %>% data.matrix
  
  score.bin <- score
  score.bin[score < 0] <- -1
  score.bin[score > 0] <- 1
  fdr <- fdr * score.bin
  
  ## extract percent overlaps
  perc.ol <- rdesc[ , grep('Signature.set.overlap.percent', colnames(rdesc))]

  
  
  ## extract colum meta data
  cdesc <- gct@cdesc
  if(nrow(cdesc) > 0){
    cdesc.rownames <- rownames(cdesc)
    cdesc.names <- colnames(cdesc)
    if('id' %in% cdesc.names)
      cdesc <-  data.frame(cdesc[, -which(cdesc.names == 'id')])
    colnames(cdesc) <- cdesc.names[-which(cdesc.names == 'id')]
    rownames(cdesc) <- cdesc.rownames
  }
  ############################################
  #extract columns of interest
  if(!is.null(class_vector_levels)){
    
    cvl.idx <- grep( paste( paste('^',class_vector_levels, '$', sep=''), collapse='|'),  cdesc[, class_vector])
    if(length(cvl.idx) == 0)
      stop('\n\n', paste(class_vector_levels, collapse='; '), 'not found in "', class_vector, '"\n\n')
    score <- score[, cvl.idx] 
    fdr <- fdr[, cvl.idx]
    cdesc <- cdesc[cvl.idx, ]
  
  }
  
  ##############################################
  #order columns?
  gaps_col <- NULL
  if(sort_col){
    cluster_cols <- F
    cutree_cols <- NA
    ord.idx <-  order(cdesc[, class_vector])
    cdesc <- cdesc[ord.idx, ]
    fdr <- fdr[, ord.idx]
    score <- score[, ord.idx ]
    gaps_col <- table(cdesc[, class_vector]) %>% cumsum
  }
  
  #############################################
  #extract significant signatures
  colnames(fdr) <- paste(cdesc[, class_vector], colnames(fdr))
  #gs.signif <- which(apply(fdr, 1, function(x) { sapply( cdesc[, class_vector] %>% unique , function(xx) ifelse( sum(x[grep(paste('^',xx, ' ', sep=''), names(x))] < fdr.max, na.rm=T) >= fdr.rep, 1, 0 )) }) %>% t %>% apply(., 1, sum) > 0) %>% names
 # gs.signif <- which( apply(fdr, 1, function(x) { 
  gs.signif <- apply(fdr, 1, function(x) { 
    
    sapply( 
      
      cdesc[, class_vector] %>% unique ,
      function(class){
        # extract FDRs of current class
        class.fdrs <- x[grep(paste('^',class, ' ', sep=''), names(x))]
        class.fdrs.sig.idx <- which(abs(class.fdrs) < fdr.max)
        
        if(length(class.fdrs.sig.idx) > 0) {
          class.fdrs.sig <- class.fdrs[ class.fdrs.sig.idx ]
        } else {
          return(0)
        }
        # require 'min.rep' significant samples regulated in the same direction
        ifelse( length( class.fdrs.sig)  >= fdr.rep & prod(range( class.fdrs.sig[ class.fdrs.sig != 0] )) > 0, 1, 0 ) 
        #ifelse( length( class.fdrs.sig)  >= fdr.rep, 1, 0 ) 
      })
    }) ## %>% t %>% apply(., 1, sum) > 0) %>% names
  
  if( !is.null(dim(gs.signif)) ){
    gs.signif <- which(t(gs.signif) %>% apply(., 1, sum) > 0) %>% names
  } else {
    gs.signif <- which(gs.signif > 0)
  }
  
  ##gs.signif <- which(apply(fdr, 1, function(x) { sapply( cdesc[, class_vector] %>% unique , function(xx) ifelse( sum( abs(x[grep(paste('^',xx, ' ', sep=''), names(x))]) < fdr.max & prod(range( x[grep(paste('^',xx, ' ', sep=''), names(x))]  )) > 0, na.rm=T) >= fdr.rep, 1, 0 )) }) %>% t %>% apply(., 1, sum) > 0) %>% names
  
  if(length(gs.signif) == 0)
    stop(paste('\n\nNo gene sets meeting the criteria: FDR <', fdr.max, ' in min.', fdr.rep, 'samples in "', class_vector,'"\n\n' ))
  
  score.signif <- score[gs.signif, ]
  fdr.signif <- data.frame(abs(fdr[gs.signif, ]))
  
  ## use cell notes to mark significant signatures
  cellnote <- apply( fdr.signif, 1, function(x){xx=x;xx[x >= fdr.max]='';xx[x < fdr.max]=cellnote_pch;xx[is.na(x)]='';xx }) %>% t %>% data.frame
  
  ###########################################
  ## colors
  if(is.null(max.val))
    max.val <- max(abs(score), na.rm=T)
  if(is.null(min.val))
    min.val <- -max.val
  
  score.signif[score.signif < min.val] <- min.val
  score.signif[score.signif > max.val] <- max.val
  
  if(hm_color_scheme == 'cyan2brown'){
    color.breaks = seq( min.val, max.val, length.out=13 )
    color.hm <- colorschemes$BluetoDarkOrange.12
 
  }
  if(hm_color_scheme == 'blue2red'){
    color.breaks = seq( min.val, max.val, length.out=9 )
    color.hm = rev(brewer.pal (length(color.breaks)-1, "RdBu"))
  }
  
  
  if(!is.null( exclude.from.cdesc)){
    
    rm.idx <- which(colnames(cdesc)  %in% exclude.from.cdesc)
    if(length(rm.idx) > 0)
      cdesc <- cdesc[, -rm.idx]
    
  }
    
  ## plot
  #pheatmap(score.signif, annotation_col=cdesc, annotation_row=perc.ol,
  #         annotation_color=anno_colors, display_numbers = cellnote, cutree_cols = cutree_cols, fontsize_number = cellnote_size, number_color= cellnote_color, main=main, color = color.hm, breaks = color.breaks, filename=filename, 
  #         cluster_cols=cluster_cols, cluster_rows=cluster_rows, gaps_col=gaps_col, ...)
  pheatmap(score.signif, annotation_col=cdesc,
         annotation_color=anno_colors, display_numbers = cellnote, cutree_cols = cutree_cols, fontsize_number = cellnote_size, number_color= cellnote_color, main=main, color = color.hm, breaks = color.breaks, filename=filename, 
         cluster_cols=cluster_cols, cluster_rows=cluster_rows, gaps_col=gaps_col, ...)

  return(0)
  }
