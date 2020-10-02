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
fdr.sig <- as.numeric(args[4])
tmp.dir <- args[5]

## ###################################################################
##      create a Rmarkdown report for RNA-proteome correlation
##
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
## tmp.dir    - folder used to untar and write output
## type       - protoeme/phospho
## fdr.sig    - significance level for correlations
rmd_rna_seq_correlation <- function(tar.file, label='pipeline-test', type='proteome', fdr.sig=0.05, tmp.dir){
    
    topN <- 10
    gene.col <- 'geneSymbol'
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
    save(gene.col, topN,label, rna.cor, fdr.sig, file=paste(tmp.dir, 'rna_cor.RData', sep='/'))



    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)

    #rmd <- '\n---\ntitle: "RNA correlation"\noutput:\n\thtml_document:\n\t\tdf_print: paged\n---\n'
    rmd <- ''
    rmd <- c(paste(rmd, "\n# RNA-`r type` correlation - ", label,"\n
This report summarizes the correlation analysis between ```RNA``` and ``` `r type` ``` expression. Correlations are calculated per gene (across samples) and per sample (across genes) using Pearson's correlation coefficient.

```{r echo=F, warning=F, message=F}\n
require(plotly)
require(knitr)
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

#### Summary table

No. of correlations | all | negative | positive | median cor.
---- | --- | -------- | --------- | -------------------
all | `r sum(!is.na(rna.cor$correlation))` | `r sum(rna.cor$correlation < 0, na.rm=T)` | `r sum(rna.cor$correlation >= 0, na.rm=T)` | `r median(rna.cor$correlation, na.rm=T)`
FDR `r round(fdr.sig*100)`% |`r sum(!is.na(rna.cor.plot.sig$correlation), na.rm=T)` | `r sum(rna.cor.plot.sig$correlation < 0, na.rm=T)` | `r sum(rna.cor.plot.sig$correlation >= 0, na.rm=T)` | `r median(rna.cor.plot.sig$correlation, na.rm=T)`

***
                   
#### Correlation-rank-plot

Genes are ranked according to their ```RNA``` / ``` `r type` ``` correlation (x-axis) and plotted against the actual correlation value (y-axis). Gene names can be inferred by hovering over the data points.

```{r echo=F, warning=F, message=F}\n

## ##############################
## rank plot
rna.cor <- data.frame(rank=rank(-rna.cor$correlation), rna.cor, stringsAsFactors=F)
crank <- plot_ly(rna.cor, x=~rank, y=~correlation,  type='scatter', mode='markers', text=~geneSymbol)
crank
#save(rna.cor, file='test.RData')
top.pos <- data.frame( rna.cor[order(rna.cor$correlation, decreasing=T)[1:topN], c('", gene.col, "', 'correlation' ,'p.value', 'p.value.fdr')], rep('    ', topN) )
colnames(top.pos)[ncol(top.pos)] <- '      '
top.neg <- rna.cor[order(rna.cor$correlation, decreasing=F)[1:topN], c('", gene.col, "', 'correlation' ,'p.value', 'p.value.fdr')]
```\n
\n\n

#### Top ```r topN``` correlated genes
```{r echo=F, warning=F, message=F, results='as.is'}\n
kable(list(top.pos, top.neg), row.names=F, caption='Left: positive correlation; right: negative correlation ')
```

***

### Sample-wise correlations
For each of the `r length(samp)` samples Pearson correlations between ```RNA``` / ``` `r type` ``` across all genes were calculated and are depicted as bubble chart. The size of each bubble scales with correlation value.
```{r echo=F, warning=F, message=F}\n
samp.cor=data.frame(sample=samp, correlation=samp.cor, stringsAsFactors=F)
Median <- median(samp.cor$correlation, na.rm=T)

plot_ly(samp.cor, x=~sample, y=~correlation,  type='scatter', mode='markers', text=~sample, size=~correlation+1.5) %>%
layout( shapes=list(type='line', x0=0, x1=length(samp)-1, y0=Median, y1=Median, line=list(dash='dot', width=1)), title=paste('Sample correlations ( median: ', round(Median,3),')', sep='' ))
#%>% add_lines(y=median(~correlation, na.rm=T), colors='black')
rm(Median)
```\n
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
    writeLines(rmd, con=paste(tmp.dir, paste('rna-corr_', label,'.rmd', sep=''), sep='/'))
    rmarkdown::render(paste(tmp.dir, paste('rna-corr_', label,'.rmd', sep=''), sep='/'))


    file.copy(paste(tmp.dir, paste('rna-corr_', label,'.html', sep=''), sep='/'), paste(wd, paste('rna-corr_', label,'.html', sep=''), sep='/'))
    cat('## Output written to:', paste(wd, paste('rna-corr_', label,'.html', sep=''), sep='/'), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )


}

## ################################################
## run
rmd_rna_seq_correlation(tar.file=tar.file, label=label, type=type, fdr.sig=fdr.sig, tmp.dir=tmp.dir)

