#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
##tar.file <- args[1]
tar.file <- '/media/sf_Dropbox/Devel/PGDAC/test/input/inputSM-output.tar'

## used to extract tar file
tmp.dir <- 'tmp2' 
type='proteome'

label='pipeline-test'

##source('/prot/proteomics/r-util')
source('/media/sf_Karsten/R-code/gct-io.r')



## ###################################################################
##      create a Rmarkdown report for normalization
##
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of
## tmp.dir    - folder used to untar and write output
## type       - 
rmd_normalize <- function(tar.file, label='pipeline-test', type='proteome', tmp.dir){    
    
    ## prepare log file
    logfile=paste(label, '_rmd-normalize.log', sep='')
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_normalize\'--\n\n', file=logfile ) 
    cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
    
    require(pacman)
    p_load(rmarkdown)

    ## #################################
    ## extract tar ball
    if(!dir.exists(tmp.dir))
        dir.create(tmp.dir)
    cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
    untar(tar.file, exdir=tmp.dir)
  
    ## ##################################
    ## gather data for markdown file
    ## normalized data
    gct.norm.str <- paste(tmp.dir, '/', label, '/normalized-data/', type, '-ratio-norm.gct', sep='')
    cat('## Importing file ', gct.norm.str, '\n', file=logfile, append=T)

    gct.norm <- parse.gctx( gct.norm.str )
    save(gct.norm, file=paste(tmp.dir, 'norm.RData', sep='/'))                                  
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- c(paste("# Normalization\n
This document desrcibes the results of sample-wise normalization applied to expression values.
```{r echo=F, warning=F}\n
cat(getwd())
require(plotly)
## ##############################
## prepare data set
load(paste('norm.RData', sep='/')) ## import data 
m=gct.norm@mat     ## expression
rd=gct.norm@rdesc  ## row annotations
cd=gct.norm@cdesc  ## column annotations
d=data.frame(m, stringsAsFactors=F)
colnames(d) <- cd$id
rownames(d) <- rd$id

## applied normalization
norm.type <-  unique(cd$normalization)
```

### Data set
Below find a summary of the data set:\n
No. of features | No. of samples
--------------- | ---------------
`r nrow(m)`     | `r ncol(m)`


### Normalized data\n
The box-whisker plots depict the distribution of expression values after `r norm.type`-normalization. Please note that the graph is interactive.
```{r echo=F, warning=F}
## ##############################
## box plots
##plot(m[,1], m[,2])
##plot_ly(d, x=~x, y=~y, mode='markers', type='scatter', text=~Feature, colors=as.character('blue'))\n
p <- plot_ly(d, y=d[, 1],  type='box', name=colnames(d)[1])\n
for(i in 2:ncol(d))
   p <- p %>% add_trace(y=d[, i], name=colnames(d)[i])
p
```\n

### Normaliziation coefficients
```{r echo=F, warnings=F}
norm.coeff <- data.frame(
   sample=cd$id,
   center=cd$normalization.center,
   scale=cd$normalization.scale, stringsAsFactors=F)
plot_ly(norm.coeff, y=~center, type='scatter', mode='lines+markers', text=~sample, hoverinfo='y+text', name='center') %>% add_trace(y=~scale, name='scale')
```

", sep=''))

    cat('## Rendering Rmarkdown file\n', file=logfile, append=T)
    writeLines(rmd, con=paste(tmp.dir, 'norm.rmd', sep='/'))
    rmarkdown::render(paste(tmp.dir, 'norm.rmd', sep='/'))
    cat('## Output written to:', paste(tmp.dir, 'norm.html', sep='/'), file=logfile, append=T)

    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T ) 
}


## ################################################
## run
rmd_normalize(tar.file=tar.file, label=label, type=type, tmp.dir=tmp.dir)

