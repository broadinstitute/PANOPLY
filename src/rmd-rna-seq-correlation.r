#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
tar.file <- args[1]
label <- args[2]
type <- args[3]
tmp.dir <- args[4]

## local testing
tar.file <- '/media/sf_Dropbox/Devel/PGDAC/test/input/RNAcorr-output.tar'
tmp.dir <- 'tmp2'
type='proteome'
label='pipeline-test'
##source('/media/sf_Karsten/R-code/gct-io.r')

##source('/prot/proteomics/Projects/R-utilities/gct-io.r')

## ###################################################################
##      create a Rmarkdown report for RNA-proteome correlation
##
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output
## type       - pr
rmd_rna_seq_correlation <- function(tar.file, label='pipeline-test', type='proteome', tmp.dir){    

    wd <- getwd()
    
    ## prepare log file
    logfile=paste(label, '_rmd-rna-seq-correlation.log', sep='')
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_rna_seq_correlation\'--\n\n', file=logfile ) 
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
    ## table with rna-protein correlations
    rna.cor.str <- paste(tmp.dir, '/', label, '/rna/', type, '-mrna-cor.tsv', sep='')
    cat('## Importing file ', rna.cor.str, '\n', file=logfile, append=T)
    rna.cor <- data.frame( read.delim(rna.cor.str, stringsAsFactors = F), stringsAsFactors = F)
    save(rna.cor, file=paste(tmp.dir, 'rna_cor.RData', sep='/'))                                  
    
    
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- c(paste("\n# RNA-`r type` correlation\n
tets
```{r echo=F, warning=F, message=F}\n
require(plotly)
## ##############################
## prepare data set
load(paste('rna_cor.RData', sep='/')) ## import data 

## ##############################
## histogram
pcorr <- plot_ly(rna.cor, x=~correlation,  type='histogram', name='Correlation')
pcorr
##for(i in 2:ncol(du))
##   pu <- pu %>% add_trace(y=du[, i], name=colnames(du)[i])
##pu
```\n
\n\n

***
", sep=''))

    cat('## Rendering Rmarkdown file\n', file=logfile, append=T)
    writeLines(rmd, con=paste(tmp.dir, 'rna-corr.rmd', sep='/'))
    rmarkdown::render(paste(tmp.dir, 'rna-corr.rmd', sep='/'))
    file.copy(paste(tmp.dir, 'rna-corr.html', sep='/'), paste(wd, 'rna-corr.html', sep='/'))
    
    cat('## Output written to:', paste(wd, 'rna-corr.html', sep='/'), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )


}


## ################################################
## run
rmd_normalize(tar.file=tar.file, label=label, type=type, tmp.dir=tmp.dir)

