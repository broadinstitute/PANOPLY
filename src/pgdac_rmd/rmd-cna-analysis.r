#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
tar.file <- args[1]
label <- args[2]
type <- args[3]
fdr.sig <- as.numeric(args[4])
tmp.dir <- args[5]

require(pacman)
p_load(rmarkdown)
p_load(cmapR)
p_load(morpheus)

## local testing
#tar.file <- '~/pgdac_test_data/acetylproteome-data-v1.03011-final.tar'
##tmp.dir <- 'tmp2'
#type <- 'acetylproteome'
##type='proteome'
#label='acetylproteome-data-v1.03011'


## ###################################################################
##      create a Rmarkdown report for CNA analysis
##
## tar.file   - url of tar file
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output

rmd_cna_analysis <- function(tar.file, label='pipeline-test', type, fdr.sig=0.05, tmp.dir, cna.dir='cna'){    

    wd <- getwd()
    
    ## prepare log file
    logfile=paste(label, '_rmd-cna-analysis.log', sep='')
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_cna_analysis\'--\n\n', file=logfile ) 
    cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
    
 
    ## #################################
    ## extract tar ball
    if(!dir.exists(tmp.dir))
        dir.create(tmp.dir)
    cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
    untar(tar.file, exdir=tmp.dir)
  
    
    ## ##################################
    ## copy CNA image
    png.str <- paste(tmp.dir, label, cna.dir, 'all-cna-plot.png', sep='/')
    file.copy( png.str, tmp.dir)
    
    ## ##################################
    ## gather data
    rna.cna <- read.csv(paste(tmp.dir, label, cna.dir, 'all-mrna-vs-cna-sigevents.csv', sep='/'))
    pome.cna <- read.csv(paste(tmp.dir, label, cna.dir, 'all-pome-vs-cna-sigevents.csv', sep='/'))
    save(rna.cna, pome.cna, file=paste(tmp.dir, 'cna.RData', sep='/')) 
    
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    #rmd <- '\n---\ntitle: "CNA correlation"\noutput:\n\thtml_document:\n\t\tdf_print: paged\n---\n'
    rmd <- ''
    rmd <- paste(rmd, '\n
                 
\n<script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>

\n<style>
\n.zoomDiv {
\n opacity: 0;
\n position:absolute;
\n top: 50%;
\n left: 50%;
\n z-index: 50;
\n transform: translate(-50%, -50%);
\n box-shadow: 0px 0px 50px #888888;
\n max-height:100%; 
\n overflow: scroll;
\n}

\n.zoomImg {
\n width: 100%;
\n}
\n</style>
                 
                 
\n<script type="text/javascript">
\n$(document).ready(function() {
\n $(\'slides\').prepend("<div class=\"zoomDiv\"><img src=\"\" class=\"zoomImg\"></div>");
\n // onClick function for all plots (img\'s)
\n $(\'img:not(.zoomImg)\').click(function() {
\n $(\'.zoomImg\').attr(\'src\', $(this).attr(\'src\'));
\n $(\'.zoomDiv\').css({opacity: \'1\', width: \'60%\'});
\n});
\n // onClick function for zoomImg
\n $(\'img.zoomImg\').click(function() {
\n $(\'.zoomDiv\').css({opacity: \'0\', width: \'0%\'});
\n });
\n});
\n</script>
\n')
    
    rmd <- c(paste(rmd, "\n# CNA analysis - ", label, "\n
This document describes the results of the CNA analysis.
```{r echo=F, warning=F, message=F}\n
library(pacman)
p_load(plotly)
p_load(dplyr)
## ##############################
## prepare data set
load(paste('cna.RData', sep='/')) ## import data 
```
\n\n

***

\n\n#### CNA *cis* and *trans* correlations\n

![Figure 1. Correlations of CNA (*x* axes) to ```RNA``` and ```r type``` expression levels (*y* axes). Significant (FDR < ```r fdr.sig```) positive (red) and negative (green) are indicated. CNA *cis*-effects appear as a red diagonal line, CNA *trans*-effects as vertical stripes. The histograms show the number of significant CNA *trans*-effects for each CNA gene.
](all-cna-plot.png)
\n

***


\n\n#### Summary table\n


Type  | No. of genes | No. significant genes  | No. significant events | No. significant *cis* effects
------------ | --------------- | --------------- | --------------------- | ------------------------
RNA | `r length(unique(rna.cna$HGNCsymbol))` | `r sum(rna.cna$SignificantEvents > 0, na.rm=T)`     | `r sum(rna.cna$SignificantEvents, na.rm=T)`     | `r sum(rna.cna$SignificantCisEffect > 0, na.rm=T)`
`r type` | `r length(unique(rna.cna$HGNCsymbol))`  | `r sum(pome.cna$SignificantEvents > 0, na.rm=T)`     | `r sum(pome.cna$SignificantEvents, na.rm=T)`     | `r sum(pome.cna$SignificantCisEffect > 0, na.rm=T)`

\n\n
                   

"))
    
    
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
    writeLines(rmd, con=paste(tmp.dir, paste('cna-analysis_', label,'.rmd', sep=''), sep='/'))
    rmarkdown::render(paste(tmp.dir, paste('cna-analysis_', label,'.rmd', sep=''), sep='/'))
    file.copy(paste(tmp.dir, paste('cna-analysis_', label,'.html', sep=''), sep='/'), paste(wd, paste('cna-analysis_', label,'.html', sep=''), sep='/'))
    
    cat('## Output written to:', paste(wd, paste('cna-analysis_', label,'.html', sep=''), sep='/'), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )


}


## ################################################
## run
rmd_cna_analysis(tar.file=tar.file, label=label, tmp.dir=tmp.dir, type=type)

