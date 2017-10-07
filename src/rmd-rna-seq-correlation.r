#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
tar.file <- args[1]
label <- args[2]
type <- args[3]
fdr.sig <- as.numeric(args[4]) 
tmp.dir <- args[5]

## local testing
#tar.file <- '/media/sf_Dropbox/Devel/PGDAC/test/input/RNAcorr-output.tar'
#tmp.dir <- 'tmp2'
#type='proteome'
#label='proteome-medullo'
##source('/media/sf_Karsten/R-code/gct-io.r')

## ###################################################################
##      create a Rmarkdown report for RNA-proteome correlation
##
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output
## type       - protoeme/phospho
## fdr.sig    - significance level for correlations
rmd_rna_seq_correlation <- function(tar.file, label='pipeline-test', type='proteome', fdr.sig=0.05, tmp.dir){    

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
    save(rna.cor, fdr.sig, file=paste(tmp.dir, 'rna_cor.RData', sep='/'))                                  
    
    
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- c(paste("\n# RNA-`r type` correlation\n
This report summarizes the correlation analysis between ```RNA``` and ``` `r type` ``` expression. Correlations are calculated per gene (across samples) and per sample (across genes) using Pearson's correlation coefficient.

```{r echo=F, warning=F, message=F}\n
require(plotly)
## ##############################
## prepare data set
load(paste('rna_cor.RData', sep='/')) ## import data 

## calculate sample correlations
samp <- unique(sub('\\\\.pome|\\\\.mrna','',grep('pome|mrna',colnames(rna.cor), value=T)))
samp.cor <- sapply(samp, function(x){
  p=rna.cor[, grep(paste(x, 'pome', sep='.'), colnames(rna.cor))]
  r=rna.cor[, grep(paste(x, 'mrna', sep='.'), colnames(rna.cor))]
  cor(p, r, use='pairwise.complete')
} )
```\n

***

### Gene-wise correlations

The histogram shows all gene-wise correlations (blue) across `r length(samp)` samples. Significant correlations (FDR adjusted p-values < `r fdr.sig`) are depicted in red. You can enable/disable layers of the histogram by clicking on the respective legend. The dashed line depicts the median correlation based on all genes.


```{r echo=F, warning=F, message=F}\n
## ##############################
## histogram
rna.cor.plot <- rna.cor[c('geneSymbol', 'correlation', 'p.value.fdr')]
rna.cor.plot <- data.frame(rna.cor.plot, Pearson=rep('all correlations', nrow(rna.cor.plot)))
rna.cor.plot.sig <- rna.cor.plot[which(rna.cor.plot$p.value.fdr < fdr.sig), ]
rna.cor.plot.sig$Pearson <- rep(paste('FDR ',round(fdr.sig*100),'%', sep=''), nrow(rna.cor.plot.sig))
rna.cor.plot <- rbind(rna.cor.plot, rna.cor.plot.sig)

#pcorr <- plot_ly(rna.cor, x=~correlation,  type='histogram', name='Pearson\\'s r', alpha=0.7)
#pcorr <- pcorr %>% add_trace(x=rna.cor$correlation[rna.cor$p.value.fdr < 0.05], name='FDR < 5%') %>% layout(barmode='overlay')
#pcorr <- pcorr %>% add_trace(x=median(rna.cor$correlation, na.rm=T), y=1e4, mode='lines')
#pcorr

Median <- median(rna.cor.plot$correlation, na.rm=T)
pcorr <- ggplot(rna.cor.plot, aes(x=correlation, fill=Pearson, alpha=0.7)) + 
geom_histogram(binwidth=.05, position='identity') + scale_fill_manual(values=c('blue', 'red') ) + 
geom_vline(aes(xintercept=Median), color='black', linetype='dashed')
ggplotly(pcorr)
rm(Median)

```\n
\n\n

***

#### Correlation-rank-plot

Genes are ranked according to their ```RNA``` / ``` `r type` ``` correlation (x-axis) and plotted against the actual correlation value (y-axis). Gene names can be inferred by hovering over the data points.

```{r echo=F, warning=F, message=F}\n

## ##############################
## rank plot
rna.cor <- data.frame(rank=rank(-rna.cor$correlation), rna.cor, stringsAsFactors=F)
crank <- plot_ly(rna.cor, x=~rank, y=~correlation,  type='scatter', mode='markers', text=~geneSymbol)
crank

```\n
\n\n

***
#### Summary table
The table summarizes the data shown above.

No. of correlations | all | negative | positive
---- | --- | -------- | ---------
all | `r sum(!is.na(rna.cor.plot$correlation))` | `r sum(rna.cor.plot$correlation < 0)` | `r sum(rna.cor.plot$correlation >= 0)`
FDR `r round(fdr.sig*100)`% |`r sum(!is.na(rna.cor.plot.sig$correlation))` | `r sum(rna.cor.plot.sig$correlation < 0)` | `r sum(rna.cor.plot.sig$correlation >= 0)`

***

### Sample-wise correlations
For each of the `r length(samp)` samples Pearson correlations between ```RNA``` / ``` `r type` ``` across all genes were calculated and are depicted as a bubble chart. The size of each bubble scales with correlation value. 
```{r echo=F, warning=F, message=F}\n
samp.cor=data.frame(sample=samp, correlation=samp.cor, stringsAsFactors=F)
Median <- median(samp.cor$correlation, na.rm=T)

plot_ly(samp.cor, x=~sample, y=~correlation,  type='scatter', mode='markers', text=~sample, size=~correlation+1.5) %>% 
layout( shapes=list(type='line', x0=0, x1=length(samp)-1, y0=Median, y1=Median, line=list(dash='dot', width=1)), title='Sample correlations')
#%>% add_lines(y=median(~correlation, na.rm=T), colors='black')
rm(Median)
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
rmd_rna_seq_correlation(tar.file=tar.file, label=label, type=type, fdr.sig=fdr.sig, tmp.dir=tmp.dir)

