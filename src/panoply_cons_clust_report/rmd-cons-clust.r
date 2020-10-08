#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
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
p_load(morpheus)
p_load(dplyr)

## ###################################################################
##      create a Rmarkdown report for consensus clustering results
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output

rmd_cons_clust <- function(tar.file, label='pipeline-test', label.rmd='cons-clust', type, tmp.dir, res.dir='clustering'){    

    wd <- getwd()
    
    ## prepare log file
    logfile=paste0(label, '_rmd-', label.rmd)
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_', label.rmd,'\'--\n\n', file=logfile ) 
    cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
    
 
    ## #################################
    ## extract tar ball
    if(!dir.exists(tmp.dir))
        dir.create(tmp.dir)
    cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
    untar(tar.file, exdir=tmp.dir)
  
    
    ## ##################################
    ## gather data for markdown file
    ## - figures
    
    ## copy consenus matrix
    fig.str <- file.path(tmp.dir, label, res.dir, 'acetylome_consensus_matrix_K3.png')
    
    
    file.copy( fig.str, tmp.dir)
  
    # save(label, est.cna.gct, est.rna.gct, est.pome.gct, est, signatures, 
    #      data.rna, data.pome, data.cna, 
    #      cm.rna.cna, cm.pome.cna, cm.pome.rna, 
    save(label, label.rmd, fig.str, 
           file=file.path( tmp.dir, 'data.RData'))                                  
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- paste("\n# Consensus Clustering - ", label,"\n
This document describes the results of the consenus clustering module.
```{r echo=F, warning=F, message=F}\n
library(pacman)
p_load(plotly)
p_load(dplyr)

## ##############################
## prepare data set
load(paste('data.RData', sep='/')) ## import data 

```
\n\n

***

\n\n### Consensus matrix

```{r, include=TRUE, fig.align='center', fig.cap=c('your caption'), echo=FALSE}
#knitr::include_graphics('./acetylome_consensus_matrix_K3.png')
knitr::include_graphics(fig.str)
```

XXX
", sep="")
    
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
    writeLines(rmd, con=paste(tmp.dir, paste0( label.rmd, '_', label,'.rmd'), sep='/'))
    rmarkdown::render(paste(tmp.dir, paste(label.rmd, '_', label,'.rmd', sep=''), sep='/'))
    file.copy(paste(tmp.dir, paste(label.rmd, '_', label,'.html', sep=''), sep='/'), paste(wd, paste(label.rmd, '_', label,'.html', sep=''), sep='/'))
    
    cat('## Output written to:', paste(wd, paste(label.rmd, '_', label,'.html', sep=''), sep='/'), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )


}


## ################################################
## run
rmd_sample_qc(tar.file=tar.file, label=label, tmp.dir=tmp.dir, type=type)

