#!/usr/bin/env Rscript
options( warn = -1 )
args <- commandArgs(trailingOnly=T)

## get arguments
tar.file <- args[1]
label <- args[2]

require(pacman)
p_load(rmarkdown)
p_load(cmapR)
p_load(morpheus)
p_load(dplyr)
p_load(yaml)
p_load(knitr)
p_load(readxl)
p_load(DT)
p_load(tibble)

## ###################################################################
##      create a Rmarkdown report for the mo_nmf module
## tar.file   - url of tar file created by task 'parse_sm_table'
## label      - character, name of folder in the tarball
rmd_mo_nmf <- function(tar.file, label='pipeline-test', fn.ws='workspace_after_NMF.RData', tmp.dir=tempdir()){    

    wd <- getwd()
    #tmp.dir <- tempdir()
      
    ## prepare log file
    logfile=paste0('rmd-mo-nmf-', label, '.log')
    start.time <- Sys.time()
    cat(paste(rep('#', 40), collapse=''),'\n##', paste0(start.time), '--\'rmd_', label,'\'--\n\n', file=logfile ) 
    cat('## parameters\ntar file:', tar.file, '\ntmp dir:', tmp.dir, '\nlabel:', label, '\nlog file:', logfile, '\n', file=logfile, append=T)
    
 
    ## #################################
    ## extract tar ball
    if(!dir.exists(tmp.dir))
        dir.create(tmp.dir)
    cat('\n## Extracting tar file to', tmp.dir, '\n', file=logfile, append=T)
    untar(tar.file, exdir=tmp.dir)
  
    ####################################
    ## get directory name
    dir_idx <- sapply(dir(tmp.dir, full.names = T), dir.exists)
    dir_name <- names(dir_idx)[which(dir_idx)]
    
    ## If - for whatever reason - there are more than one directories after
    ## extracting the tar ball, pick the one that has the file 'workspace_after_NMF.RData'.
    if(length(dir_name) > 1){
      idx <- sapply(dir_name, function(x) sum(grepl(fn.ws,  dir(x))))
      
      if(sum(idx) == 0){
        cat('## ERROR - creating dummy html\n', file=logfile, append=T)
        rmd <- '\n# Error\n'
        writeLines(rmd, con=paste(tmp.dir, paste0( label, '.rmd'), sep='/'))
        rmarkdown::render(paste(tmp.dir, paste(label, '.rmd', sep=''), sep='/'))
        file.copy(paste(tmp.dir, paste(label, '.html', sep=''), sep='/'), paste(wd, paste(label, '.html', sep=''), sep='/'))
        return(1)
      }
      dir_name <- dir_name[which(idx == 1)]
      #dir_name <- file.path(tmp.dir, dir_name)
    }
    
    ## get subdirectory, e.g. 'K_4/'  
    dir_k <- dir(dir_name, pattern = '^K_[0-9]{1,2}')
    
    ## ##################################################################
    ## gather files and data for markdown file
    ## 
    
    ## path to workspace file
    ws_str <- file.path(dir_name, fn.ws) %>% gsub('\\\\', '/', .)
    
    ## path to PNG figures
    fig.str <- c(
      ## heatmap of sample weights
      coef_map=dir(file.path(dir_name, dir_k), pattern=paste0('2.2_coefmap_pheatmap_sorted.png'), full.names = T),
      ## consensus heatmap
      cons_map=dir(file.path(dir_name, dir_k), pattern=paste0('3.1_consensusmap_nrun_[0-9]{1,5}_pheatmap.png'), full.names = T),
      ## feature barplot
      feat_barplot=dir(file.path(dir_name, dir_k), pattern=paste0('4.2_barplot_features_per_.*\\.png'), full.names = T),
      ## feature heatmap, concatenated
      feat_hm_concat=dir(file.path(dir_name, dir_k), pattern=paste0('6.0_ComplexHeatmap_ALL_features-concat.png'), full.names = T),
      ## feature heatmap
      feat_hm=dir(file.path(dir_name, dir_k), pattern=paste0('6.1_ComplexHeatmap_ALL_features.png'), full.names = T),
      ## pca, all features
      pca_all=dir(file.path(dir_name, dir_k), pattern=paste0('9.0_PCA.png'), full.names = T),
      ## pca, nmf features
      pca_nmf=dir(file.path(dir_name, dir_k), pattern=paste0('9.1_PCA_NMF_features.png'), full.names = T),
      ## silhouette plot
      sil_plot=dir(file.path(dir_name, dir_k), pattern=paste0('1.0_silhouette_K_.*.png'), full.names = T)
      )
    ## path to tables
    tab.str <- c(
      ## overrepresentation results
      enrich_all=dir(file.path(dir_name, dir_k), pattern=paste0('cluster-enrichment.txt'), full.names = T),
      ## cluster labels and sample annotation
      clin_anno=dir(file.path(dir_name, dir_k), pattern=paste0('clin_anno_nmf.txt'), full.names = T),
      ## Excel sheet with all features
      nmf_feat_xlsx= dir(file.path(dir_name, dir_k), pattern=paste0('NMF_features_N.*'), full.names = T)
    )
    
    
    ######################################
    save(label, fig.str, tab.str,
           file=file.path( tmp.dir, 'data.RData'))                                  
    load(ws_str)
    
    ## ####################################################
    ## Rmarkdown
    cat('## Generating Rmarkdown file\n', file=logfile, append=T)
    
    rmd <- paste0("---
title: NMF clustering results - ", label,"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
---

# Overview

This document describes the results of the multi-omics non-negative factorization (NMF)-based clustering module. For more information about the module please visit the [PANOPLY Wiki page](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_mo_nmf).
```{r echo=F, warning=F, message=F}\n
options(knitr.table.format = 'pipe', stringsAsFactors=F)
library(pacman)
p_load(plotly)
p_load(dplyr)
p_load(NMF)
p_load(kableExtra)
p_load(readr)

## ##############################
## prepare data set
load('data.RData') ## import data 
load('", ws_str, "')

# counter variable for figures and tables
fig_count <- 1
tab_count <- 1
```
\n\n
***
\n")
    
    ########################################
    ## input data
    rmd <- paste0(rmd, "\n
\n\n# Input data matrix\n

The data matrix subjected to NMF analysis contained ```r nrow(expr)``` features measured across ```r ncol(expr)``` samples. **Table `r tab_count`** summarizes the number of features used in the clustering and their dataype(s).

\n```{r tab_feat, echo=F, warning=F, message=F}
tab_feat <- table( sub('^(.*?)-.*', '\\\\1' ,rownames(expr)))
tab_feat <- data.frame(Type=names(tab_feat), Number=as.vector(tab_feat))
colnames(tab_feat) <- c('Type', 'Number of features')
tab_feat %>%
  kbl(caption=paste0('**Table ', tab_count, '**: Number of features used for clustering.')) %>%
  kable_paper('hover', full_width = F) %>%
  column_spec(2, width = '15em')
  ## increment
  tab_count <- tab_count + 1
\n```

***
\n", sep="")
    
    
    ########################################
    ## cluster metrics
    if(length(res.rank) > 1){
        rmd <- paste0(rmd, "\n
\n\n# Determining the number of clusters\n

To determine an optimal value **k** for the number of clusters, a range of **k** between ```r opt$kmin``` and ```r opt$kmax``` was evaluated using several metrics:

* <u>Cophenetic correlation coefficient</u> (**coph**) measuring how well the intrinsic structure of the data is recapitulated after clustering.
* <u>Dispersion coeffiecient</u> (**disp**) of the consensus matrix as defined in [Kim and Park, 2007](https://pubmed.ncbi.nlm.nih.gov/17483501/) measuring the reproducibility of the clustering across ```r opt$nrun``` random iterations.
* <u>Silhouette score</u> (**sil**) measuring how similar a sample is to its own cluster (cohesion) compared to other clusters (separation) and thus is defined for each sample. The average silhoutte score across all samples serves is calcualted for each cluster number **k**.

The metrics are summarized in **Figure `r fig_count`**. The optimal number of clusters is defined as the maximum of the product of **coph** and  **disp** between **k=`r ifelse(opt$exclude_2, max(c(opt$kmin, 3)), opt$kmin)`** and **k=`r opt$kmax`**.

\n```{r cluster_metrics, echo=F, fig.cap=paste0('**Figure ', fig_count,'**: Cluster metrics as a function of cluster numbers.')}

## cophenetic correlation
rank.coph <- sapply(res.rank, cophcor)
## dispersion of consensus matrix
rank.disp <- sapply(res.rank, dispersion)
## combine
rank.coph.disp <- rank.coph * rank.disp
## silhouette
rank.sil <- lapply(res.rank, silhouette)
rank.sil.avg <- lapply(rank.sil, function(x) tapply( x[,3], x[, 1], mean))
rank.sil.avg <- sapply(rank.sil.avg, mean)      

## plot
dat <- data.frame(coph=rank.coph, disp=rank.disp, coph.disp=rank.coph.disp, sil.avg=rank.sil.avg,
         k=opt$kmin:opt$kmax)
plot_ly(x=dat$k, y=dat$coph, type='scatter', mode='markers+lines', name='coph') %>% 
    add_trace(x=dat$k, y=dat$disp, name='disp' ) %>%
    add_trace(x=dat$k, y=dat$coph.disp, name='coph * disp' ) %>%
    add_trace(x=dat$k, y=dat$sil.avg, name='sil' ) %>%
    add_segments(x=as.numeric(rank.top) , xend =as.numeric(rank.top) , y = 0, yend = 1, name='k_opt') %>%
    layout(xaxis=list(title='Number of clusters / factorization rank'), yaxis=list(title='Score'))
\n```
\n\n
***
\n", sep="")
      rmd <- paste0(rmd, "\n        
```{r inc_fig_1, echo=F}
## increment
fig_count <- fig_count + 1
```")
    }  

    ########################################
    ## clustering results
    rmd <- paste0(rmd, "\n
\n\n# Clustering results

\n The ```r ncol(expr)``` samples were separated into ```r rank.top``` clusters. **Table `r tab_count`** summarizes the number of samples in each cluster.

\n```{r tab_clust, echo=F, warning=F, message=F}
clin_anno <- read_tsv(tab.str[['clin_anno']]) %>% select(one_of(c('Sample.ID', 'NMF.consensus', 'NMF.consensus.core')))

clust_tab <- table(clin_anno$NMF.consensus) 
clust_tab <- data.frame(clust=paste0('C', names(clust_tab)), as.vector(clust_tab))

clust_core_tab <- table(clin_anno$NMF.consensus.core)
clust_core_tab <- data.frame(clust=paste0('C', names(clust_core_tab)), as.vector(clust_core_tab))

clust_tab <- full_join(clust_tab, clust_core_tab)
colnames(clust_tab) <- c('Cluster', '# samples', '# core samples')

clust_tab %>%
  kbl(caption=paste0('**Table ', tab_count,'**: Cluster composition. The cluster core is defined as samples with a cluster membership score > ', opt$core_membership)) %>%
  kable_paper('hover', full_width = F) %>%
   column_spec(2:ncol(clust_tab), width = '10em')
   
## increment
tab_count <- tab_count + 1
\n```

\n\n
***
\n")
    
    ########################################
    ## sample coefficient heatmap
    rmd <- paste0(rmd, "\n
\n## Sample coefficient matrix\n

The heatmap shown in **Figure `r fig_count`** is a visualization of the meta-feature matrix derived from decomposing the input matrix, normalized per column by the maximum entry. The matrix presents one of the main results of NMF as it provides the basis of assigning samples to clusters.  

```{r, include=TRUE, fig.align='left', fig.cap=paste0('**Figure ', fig_count,'**: Heatmap depicting the relative contributions of each sample (x-axis) to each cluster (y-axis). Samples are ordered by cluster and cluster membership score in decreasing order.'), echo=FALSE}
knitr::include_graphics(fig.str[['coef_map']])
```
\n\n
***
\n")
    
    rmd <- paste0(rmd, "\n        
```{r inc_fig_2, echo=F}
## increment
fig_count <- fig_count + 1
```\n")

    ########################################
    ## enrichment results
    rmd <- paste0(rmd, "\n
\n## Overrepresentation analysis\n

**Table `r tab_count`** summarizes the results of an overpresentation analysis of sample metadata terms (e.g. clinial annotation, inferred phenotypes, etc.) in each cluster. Shown are nominal p-values derived from a Fisher's exact test (<span style=\"background-color:#90ee90\">p<0.01</span>, <span style=\"background-color:#ffff00\">0.01<p<0.02</span>, <span style=\"background-color:#ffa500\">0.02<p<0.05</span>). Only samples with cluster memebrship score > ```r opt$core_membership``` were used to characterize the clusters.

```{r, include=TRUE, echo=FALSE, warning=T, message=F}
tab_enrich <- read_tsv(tab.str[['enrich_all']])
colnames(tab_enrich)[1] <- ''
colnames(tab_enrich)[2:ncol(tab_enrich)] <- paste0('C',colnames(tab_enrich)[2:ncol(tab_enrich)])
## remove terms with p>0.05
keep <- which( apply(tab_enrich[, 2:ncol(tab_enrich)], 1, function(x) sum(x < 0.05)/length(x)) > 0)
if(length(keep) > 0){
    if(length(keep) == 1){
      tab_enrich_2 <- matrix(tab_enrich[keep, ], nrow=1)
      colnames(tab_enrich_2) <- colnames(tab_enrich)
      tab_enrich <- tab_enrich_2
    } else {
      tab_enrich <- tab_enrich[keep, ]
    }
    ## colors for table cells
    ## p<0.01: green, 0.01<p<0.02: yellow, 0.02<p<0.05: orange
    tab_enrich_cols <- apply(tab_enrich[, 2:ncol(tab_enrich)], 2, function(x){
      res=rep('white', length(x))
      res[x < 0.01]='lightgreen'
      res[x >= 0.01 & x < 0.02]='yellow'
      res[x >= 0.02 & x < 0.05]='orange'
      res
    })
    ## insert table
    p <- tab_enrich %>%
      kbl(caption=paste0('**Table ', tab_count,'**: Overrepresentation analysis of sample metadata terms in each cluster.')) %>%
      kable_paper('hover', full_width = F) 
    for(i in 2:ncol(tab_enrich))  
      p <- p %>% column_spec(i, background = tab_enrich_cols[, (i-1)])
} else {
p <- 'Nothing to show.'
}
## increment
tab_count <- tab_count + 1
p
```
\n\n
***
\n", sep="")
    
    
    
    ########################################
    ## cluster-specific features
    rmd <- paste0(rmd, "\n# Cluster-specifc features
    
Matrix W containing the weights of each feature in a certain cluster was used to derive a list of ```r ``` representative features separating the clusters using the method proposed in ([Kim and Park, 2007](https://pubmed.ncbi.nlm.nih.gov/17483501/)). In order to derive a p-value for each cluster-specific feature, a 2-sample moderated t-test ([Ritchie et al., 2015](https://pubmed.ncbi.nlm.nih.gov/25605792/)) was used to compare the abundance of the features between the respective cluster and all other clusters. Derived p-values were adjusted for multiple hypothesis testing using the methods proposed in ([Benjamini and Hochberg, 1995](https://www.jstor.org/stable/2346101?seq=1)). Features with FDR <`r opt$feat_fdr`are used in subsequent analyses.   
")
    
    ## all features, concatenated and ordered by ccluster
    if('feat_hm_concat' %in% names(fig.str)){
      rmd <- paste0(rmd, "\n
```{r, include=TRUE, fig.align='left', fig.cap=paste0('**Figure ', fig_count,'**: Heatmap depicting abundances of cluster specific features defined as descibed above. Samples are ordered by cluster and cluster membership score in decreasing order.'), echo=FALSE}
knitr::include_graphics(fig.str[['feat_hm_concat']])
```
```{r inc_fig_3, echo=F}
## increment
fig_count <- fig_count + 1
```
\n\n
***
\n")} else if('feat_hm' %in% names(fig.str)){
  ## all features
  rmd <- paste0(rmd, "\n
```{r, include=TRUE, fig.align='left', fig.cap=paste0('**Figure ', fig_count,'**: Heatmap depicting abundances of cluster specific features defined as descibed above. Samples are ordered by cluster and cluster membership score in decreasing order.'), echo=FALSE}
knitr::include_graphics(fig.str[['feat_hm']])
```
```{r inc_fig_3, echo=F}
## increment
fig_count <- fig_count + 1
```
\n\n
***
\n")
} else {
  rmd <- paste0(rmd, "\n

**No cluster-specific fetaureas found.**

")
  }
    
  ########################################
  ## import feature xlsx
  if('nmf_feat_xlsx' %in% names(tab.str)){   
      
     rmd <- paste0(rmd, "
```{r import_feat_xlsx, echo=F, message=F}
## get all sheets
nmf_xlsx_sheets <- excel_sheets(tab.str[['nmf_feat_xlsx']])
## import all sheets
nmf_xlsx <- lapply(nmf_xlsx_sheets, function(x) read_excel(tab.str[['nmf_feat_xlsx']], sheet=x) )
names(nmf_xlsx) <- nmf_xlsx_sheets

##
cols <- c('Accession', 'Type', 'Direction','SYMBOL', 'ENZYME', 'CYTOBAND', 'NMF.Score')
nmf_xlsx_slim <- lapply(nmf_xlsx, function(x) dplyr::select(x, one_of(cols)))
nmf_xlsx_slim <- lapply(names(nmf_xlsx), function(x) add_column(nmf_xlsx_slim[[x]], Cluster=paste0('C', x), .before=1))
nmf_xlsx_comb <- Reduce('rbind', nmf_xlsx_slim)
write.table(nmf_xlsx_comb, sep='\t', file='debug.txt')
```
\n\n
***
\n")
     
  }
    
  ## feature barchart
  if('feat_barplot' %in% names(fig.str)){ 
        
    #######################################
    ## barchart of NMF features
    rmd <- paste0(rmd, "

In total ```r nrow(nmf_xlsx_comb)``` features separating the clusters have been detected using the method descibed above. The distribution of features across the different clusters are shown in **Figure `r fig_count`**. 
    
```{r feat_barplot, include=TRUE, fig.align='left', fig.cap=paste0('**Figure ', fig_count,'**: Barpchart depicting the number of cluster specific features'), echo=FALSE, out.width='50%'}
knitr::include_graphics(fig.str[['feat_barplot']])
```

```{r inc_fig_4, echo=F}
## increment
fig_count <- fig_count + 1
```\n

\n\n
***
\n")
    } ## end if feature barplot is present
  
  
  if('nmf_feat_xlsx' %in% names(tab.str)){
    ######################################
    ## display xlsx    
    rmd <- paste0(rmd, "
The data table below depicts all cluster specific features. The table is interactive and can be sorted and filtered. Please note that the table represents a condensed verison of the entire table which can be found the Excel sheet ```r sub('.*/', '',tab.str[['nmf_feat_xlsx']])```

```{r display_feat_xlsx, echo=F}
DT::datatable(nmf_xlsx_comb, width='1000', escape=F, filter='top', rownames=FALSE,
                  options = list( pageLength = 10, scrollX = T, selection='none', 
                                  autoWidth = F, paging=T, searchHighlight = TRUE,
                                  initComplete = JS(
                                      \"function(settings, json) {\",
                                      \"$(this.api().table().header()).css({'font-size': '80%'});\",
                                  \"}\"))
) %>% DT::formatStyle(colnames(nmf_xlsx_comb), fontSize='80%')
```
\n\n
***
\n")
    
  } ## end if features were found   
    
    
    
    ########################################
    ## consensus matrix
    rmd <- paste0(rmd, "\n# Cluster stability
\n\n## Consensus matrix\n

 The entries in the sample-by-samle matrix shown in **Figure `r fig_count`** depict the relative frequences with which two samples were assigned to the same cluster across ```r opt$nrun``` iterations.


```{r, include=TRUE, fig.align='center', fig.cap=paste0('**Figure ', fig_count, '**: Consensus matrix derived from ', opt$nrun,' randomly initialized iterations.'), echo=FALSE}
knitr::include_graphics(fig.str[['cons_map']])
```


```{r inc_fig_5, echo=F}
## increment
fig_count <- fig_count + 1
```\n


\n\n
***
\n")
    
    ###########################
    ## silhouette plot
    if('sil_plot' %in% names(fig.str)){
      rmd <- paste0(rmd, "\n## Silhouette plot
      
      The silhoutte plot shown in **Figure `fig_count`** depicts the consistency of the derived clusters. Samples with negative silhouette score indicate outliers in the respectvie cluster. 
      
      ```{r feat_barplot, include=TRUE, fig.align='left', fig.cap=paste0('**Figure ', fig_count,'**: Silouette plot.'), echo=FALSE, out.width='50%'}
      knitr::include_graphics(fig.str[['sil_plot']])
      ```
    
      ```{r inc_fig_4, echo=F}
      ## increment
      fig_count <- fig_count + 1
      ```\n
      \n\n
      ***
        \n")
 
      
    }
    
    ##########################################
    ## parameters
    rmd <- paste0(rmd, "\n# Parameters
    
Details about the parameters listed in **Table `r tab_count`** can be found in the [PANOPLY WIKI](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_mo_nmf).    
    
```{r params, include=TRUE, echo=FALSE, warning=T, message=F}
tab_param <- data.frame(param=names(opt), value=unlist(opt))
## insert table
tab_param %>%
      kbl(row.names =FALSE, caption=paste0('**Table ', tab_count,'**: List of paramaeters.')) %>%
      kable_paper('hover', full_width = F) 
````

```{r inc_tab_5, echo=F}
## increment
tab_count <- tab_count + 1
```\n

\n\n
***
\n")
    
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
    
    fn_out <- paste0( 'mo-nmf-report-', label, '.rmd')
    
    cat('## Rendering Rmarkdown file\n', file=logfile, append=T)
    writeLines(rmd, con=file.path(tmp.dir, fn_out))
    rmarkdown::render(file.path(tmp.dir, fn_out), output_format = "html_document")
    file.copy(file.path(tmp.dir, sub('\\.rmd', '.html', fn_out)), wd, overwrite = T)
    
    cat('## Output written to:', file.path(wd, sub('\\.rmd', '.html', fn_out)), file=logfile, append=T)
    cat('\n\n## all done. ', format(Sys.time() - start.time), file=logfile, append=T )
}


## ################################################
## run
rmd_mo_nmf(tar.file=tar.file, label=label)
#rmd_mo_nmf(tar.file=tar.file, label=label, tmp.dir = '.')
