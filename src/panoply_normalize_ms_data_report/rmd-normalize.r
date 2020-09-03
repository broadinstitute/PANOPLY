#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
tar.file <- args[1]
label <- args[2]
type <- args[3]
tmp.dir <- args[4]

require(pacman)
p_load(rmarkdown)
p_load(cmapR)

## ###################################################################
##      create a Rmarkdown report for normalization
##
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output
## type       - proteom/phospho
rmd_normalize <- function(tar.file, label='pipeline-test', type='proteome', tmp.dir){    

    wd <- getwd()
    
    ## prepare log file
    logfile=paste(label, '_rmd-normalize.log', sep='')
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_normalize\'--\n\n', file=logfile ) 
    cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
    
 
    ## #################################
    ## extract tar ball
    if(!dir.exists(tmp.dir))
        dir.create(tmp.dir)
    cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
    untar(tar.file, exdir=tmp.dir)
  
    ## ##################################
    ## gather data for markdown file
    ## normalized data
    gct.norm.str <- paste(tmp.dir, '/', label, '/normalized-data/', type, '-ratio-norm-NArm.gct', sep='')
    cat('## Importing file ', gct.norm.str, '\n', file=logfile, append=T)

    gct.norm <- parse.gctx( gct.norm.str )
    save(label, gct.norm, file=paste(tmp.dir, 'norm.RData', sep='/'))                                  
    
    ## unnormalized data
    gct.unnorm.str <- paste(tmp.dir, '/', label, '/parsed-data/', type, '-ratio.gct', sep='')
    cat('## Importing file ', gct.unnorm.str, '\n', file=logfile, append=T)
    
    gct.unnorm <- parse.gctx( gct.unnorm.str )
    save(gct.unnorm, file=paste(tmp.dir, 'unnorm.RData', sep='/'))                                  
    
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- c(paste("\n# Normalization - ", label,"\n
This document describes the results of sample-wise normalization applied to expression data. The goal of this normalization step is to make the data comparable across all samples. Typically the sample distributions are ```centered``` around common value and ```scaled``` to harmonize the variance in the data.
```{r echo=F, warning=F, message=F}\n
require(plotly)
## ##############################
## prepare data set
load(paste('unnorm.RData', sep='/')) ## import data 
mu=gct.unnorm@mat     ## expression
rdu=gct.unnorm@rdesc  ## row annotations
cdu=gct.unnorm@cdesc  ## column annotations
du=data.frame(mu, stringsAsFactors=F)
colnames(du) <- cdu$id
rownames(du) <- rdu$id

load(paste('norm.RData', sep='/')) ## import data 
m=gct.norm@mat     ## expression
rd=gct.norm@rdesc  ## row annotations
cd=gct.norm@cdesc  ## column annotations
d=data.frame(m, stringsAsFactors=F)
colnames(d) <- cd$id
rownames(d) <- rd$id
## applied normalization
norm.type <-  unique(cd$normalization)
if(norm.type == 'median') norm.type <- 'Median-MAD'
if(norm.type == 'mean') norm.type <- 'Z-score'

## filenames
input.fn <-'test'
```
\n\n

***

### Data set

Label  | No. of features | No. of samples  | Normalization method
------------ | --------------- | --------------- | ---------------------
`r label` | `r nrow(m)`     | `r ncol(m)`     | `r norm.type`
\n\n

***

### Data distribution before normalization\n
The box-whisker plots depict the distribution of expression values before normalization. Please note that the graphs in this document are interactive.
```{r echo=F, warning=F, fig.width=12}
## ##############################
## box plots
pu <- plot_ly(du, y=du[, 1],  type='box', name=colnames(du)[1],  boxpoints = FALSE)\n
for(i in 2:ncol(du))
   pu <- pu %>% add_trace(y=du[, i], name=colnames(du)[i],  boxpoints = FALSE)
pu
```\n
\n\n

***

### Data distribution after normalization\n
The box-whisker plots depict the distribution of expression values after ``` `r norm.type` ```-normalization.
```{r echo=F, warning=F, fig.width=12}
## ##############################
## box plots
##plot(m[,1], m[,2])
##plot_ly(d, x=~x, y=~y, mode='markers', type='scatter', text=~Feature, colors=as.character('blue'))\n
p <- plot_ly(d, y=d[, 1],  type='box', name=colnames(d)[1],  boxpoints = FALSE)\n
for(i in 2:ncol(d))
   p <- p %>% add_trace(y=d[, i], name=colnames(d)[i],  boxpoints = FALSE)
p
```\n
\n\n

***

### Normaliziation coefficients
The graph below depicts the normalization coeffiecients applied to each sample, i.e. the center and scale values.
\n\n
```{r echo=F, warnings=F, fig.width=12}
norm.coeff <- data.frame(
   sample=cd$id,
   center=cd$normalization.center,
   scale=cd$normalization.scale, stringsAsFactors=F)
plot_ly(norm.coeff, y=~center, type='scatter', mode='lines+markers', text=~sample, hoverinfo='y+text', name='center') %>% add_trace(y=~scale, name='scale')
```
\n\n

", sep=''))
    
    ## ####################################
    ## Footer
    rmd <- paste(rmd, "\n\n
<br>
\n
\n***
\n
\n**_Created on ", Sys.time(),"_**
\n
\n***
\n<br>
\n", sep='')
    
    

    cat('## Rendering Rmarkdown file\n', file=logfile, append=T)
    writeLines(rmd, con=paste(tmp.dir, paste('norm_',label,'.rmd', sep=''), sep='/'))
    rmarkdown::render(paste(tmp.dir, paste('norm_',label,'.rmd', sep=''), sep='/'))
    file.copy(paste(tmp.dir, paste('norm_',label,'.html', sep=''), sep='/'), paste(wd, paste('norm_',label,'.html', sep=''), sep='/'))
    
    cat('## Output written to:', paste(wd, paste('norm_',label,'.html', sep=''), sep='/'), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )


}


## ################################################
## run
rmd_normalize(tar.file=tar.file, label=label, type=type, tmp.dir=tmp.dir)

