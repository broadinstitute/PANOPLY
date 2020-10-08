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
##      create a Rmarkdown report for normalization
##
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output

rmd_sample_qc <- function(tar.file, label='pipeline-test', type, tmp.dir, harmonize.dir='harmonized-data'){    

  
  
    wd <- getwd()
    
    ## prepare log file
    logfile=paste(label, '_rmd-sample-qc.log', sep='')
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_sample-qc\'--\n\n', file=logfile ) 
    cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
    
 
    ## #################################
    ## extract tar ball
    if(!dir.exists(tmp.dir))
        dir.create(tmp.dir)
    cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
    untar(tar.file, exdir=tmp.dir)
  
    
    ## ##################################
    ## gather data for markdown file
    
    ## estimate scores
    est.cna.gct.str <- paste(tmp.dir, '/', label, '/sample-qc/cna-estimate-scores.gct', sep='')
    est.rna.gct.str <- paste(tmp.dir, '/', label, '/sample-qc/rna-estimate-scores.gct', sep='')
    est.pome.gct.str <- paste(tmp.dir, '/', label, '/sample-qc/pome-estimate-scores.gct', sep='')
    
    ## data tables
    data.cna.str <-paste(tmp.dir, '/', label, '/', harmonize.dir, '/cna-matrix.csv' ,sep='') 
    data.rna.str <-paste(tmp.dir, '/', label, '/', harmonize.dir, '/rna-matrix.csv' ,sep='') 
    data.pome.str <-paste(tmp.dir, '/', label, '/', harmonize.dir, '/', type, '-matrix.csv' ,sep='') 
   
    #save(est.cna.gct.str, est.rna.gct.str, est.pome.gct.str, file='test.RData') 
    
    ## ####################################
    ## import data
    cat('## Importing estimate scores...', file=logfile, append=T)
    est.cna.gct <- parse.gctx( est.cna.gct.str )
    est.rna.gct <- parse.gctx( est.rna.gct.str )
    est.pome.gct <- parse.gctx( est.pome.gct.str )
    cat(' DONE\n', file=logfile, append=T)
    
    
    ## signatures
    signatures <- est.cna.gct@rid
    
    ## merge
    est.cna=t(est.cna.gct@mat)
    colnames(est.cna) <- paste( sub('.*\\.','',rownames(est.cna)[1]),  colnames(est.cna), sep='.')
    rownames(est.cna) <- sub('^(.*)\\..*', '\\1', rownames(est.cna))
    est.cna=data.frame(id=rownames(est.cna), est.cna, stringsAsFactors=F)       ## cna
    
    est.rna=t(est.rna.gct@mat)
    colnames(est.rna) <- paste( sub('.*\\.','',rownames(est.rna)[1]),  colnames(est.rna), sep='.')
    rownames(est.rna) <- sub('^(.*)\\..*', '\\1', rownames(est.rna))
    est.rna=data.frame(id=rownames(est.rna), est.rna, stringsAsFactors=F)       ## rna
    
    est.pome=t(est.pome.gct@mat)
    colnames(est.pome) <- paste( sub('.*\\.','',rownames(est.pome)[1]),  colnames(est.pome), sep='.')
    rownames(est.pome) <- sub('^(.*)\\..*', '\\1', rownames(est.pome))
    est.pome=data.frame(id=rownames(est.pome), est.pome, stringsAsFactors=F)       ## pome

    ## combine
    est <- full_join(est.pome, est.rna, 'id')
    est <- full_join(est, est.cna, 'id')
    est <- est[, -which(colnames(est) == 'id')]
    
    ## import data files
    cat('## Importing data files...', file=logfile, append=T)
    data.cna <- read.csv(data.cna.str, row.names = 1) 
    colnames(data.cna) <- paste(colnames(data.cna), 'CNA', sep='.')
    data.rna <- read.csv(data.rna.str, row.names = 1) 
    colnames(data.rna) <- paste(colnames(data.rna), 'RNA', sep='.')
    data.pome <- read.csv(data.pome.str, row.names = 1) 
    colnames(data.pome) <- paste(colnames(data.pome), 'PROT', sep='.')
    cat(' DONE\n', file=logfile, append=T)
    
    ## calculate correlations
    cm.rna.cna <- cor(data.rna, data.cna, use='pairwise.complete')
    cm.pome.cna <- cor(data.pome, data.cna, use='pairwise.complete')
    cm.pome.rna <- cor(data.pome, data.rna, use='pairwise.complete')
    
    save(label, est.cna.gct, est.rna.gct, est.pome.gct, est, signatures, 
         data.rna, data.pome, data.cna, 
         cm.rna.cna, cm.pome.cna, cm.pome.rna, 
         file=paste(tmp.dir, 'estimate.RData', sep='/'))                                  
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- c(paste("\n# Sample QC - ",label,"\n
This document describes the results of the sample QC module.
```{r echo=F, warning=F, message=F}\n
library(pacman)
p_load(plotly)
p_load(dplyr)
## ##############################
## prepare data set
load(paste('estimate.RData', sep='/')) ## import data 

```
\n\n

***

\n\n### Tumor purity\n

*E*stimation of *ST*romal amd *I*mmune cells in *MA*lignant *T*umor tissues ussing *E*xpression data (*ESTIMATE*) is a computational tool to predict presence of infiltrating stromal/immune cells in tumor tissues using expression data. The tool is based on single sample Gene Set Enrichment Analysis (ssGSEA) and scores gene set signatures of stroma, immune cells and tumor purity. More information about the tool can be found [here](https://bioinformatics.mdanderson.org/main/ESTIMATE:Overview).

The box-whisker plots depict distributions of *ESTIMATE* scores in the different *omic* expression data.


"))
  
    ## boxplots for each signature
    for(i in 1:length(signatures)){
    ## ###########################################
    ## ESTIMATE
    rmd <- paste( rmd, "

```{r echo=F, warning=F}
## ##############################
## box plots
signature.idx <- grep(signatures[",i,"], colnames(est))
p.est <- plot_ly(est, y=est[, signature.idx[1]],  type='box', name=colnames(est)[signature.idx[1]])\n
for(i in 2:length(signature.idx))
   p.est <- p.est %>% add_trace(y=est[, signature.idx[i]], name=colnames(est)[signature.idx[i]])
p.est <- p.est %>% layout(title=signatures[",i,"])
p.est
```\n
\n\n
***
\n\n
", sep='')

    
    }
    
    ## #############################################################
    ##  correlations
    rmd <- paste( rmd, "\n\n
### Sample correlations between *omes*

The heatmap below depicts correlations (*Pearson\'s r*) calculated for each pair of *ome* data available for each sample. The main purpose for this analysis is to identify potential sample mix-ups during sample preparation.

\n
### RNA-PROT correlation
```{r echo=F, fig.height=10, fig.width=10}
    morpheus(cm.pome.rna, symm = TRUE, colorScheme = list(scalingMode='fixed', stepped=FALSE, values=c(-1,0,1), colors=c('blue', 'white', 'red')), columnSize = 'fit', rowSize = 'fit')
```\n
<br>
<br>
<br>
\n
***
\n

### CNA-RNA correlation
```{r echo=F, fig.height=10, fig.width=10}
    morpheus(cm.rna.cna, symm = TRUE, colorScheme = list(scalingMode='fixed', stepped=FALSE, values=c(-1,0,1), colors=c('blue', 'white', 'red')), columnSize = 'fit', rowSize = 'fit')
```\n
<br>
<br>
<br>
\n
***
\n
### CNA-PROT correlation
```{r echo=F, fig.height=10, fig.width=10}
    morpheus(cm.pome.cna, symm = TRUE, colorScheme = list(scalingMode='fixed', stepped=FALSE, values=c(-1,0,1), colors=c('blue', 'white', 'red')), columnSize = 'fit', rowSize = 'fit')
```\n
<br>
<br>
<br>
\n
***
\n
\n
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
    writeLines(rmd, con=paste(tmp.dir, paste('sample-qc_', label,'.rmd', sep=''), sep='/'))
    rmarkdown::render(paste(tmp.dir, paste('sample-qc_', label,'.rmd', sep=''), sep='/'))
    file.copy(paste(tmp.dir, paste('sample-qc_', label,'.html', sep=''), sep='/'), paste(wd, paste('sample-qc_', label,'.html', sep=''), sep='/'))
    
    cat('## Output written to:', paste(wd, paste('sample-qc_', label,'.html', sep=''), sep='/'), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )


}


## ################################################
## run
rmd_sample_qc(tar.file=tar.file, label=label, tmp.dir=tmp.dir, type=type)

