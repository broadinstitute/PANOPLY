#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for consensus clustering module ###

# tar_file  - URL of tar file created by task panoply_cons_clust
# yaml_file - URL of master parameters yaml file including output from startup notebook
# label     - character, name of folder in tarball
# type      - character, data type

args = commandArgs(TRUE)

tar_file = args[1]
yaml_file = args[2]
label = args[3]
type = args[4]

library(rmarkdown)
library(yaml)
library(dplyr)
library(stringr)
library(glue)
library(ggplot)
library(cmapR)
library(factoextra)

rmd_cons_clust <- function(tar_file, label='pipeline-test', label.rmd='cons-clust', type, tmp.dir, res.dir='clustering'){

  # extract files from tarball
  untar(tar_file)
  
  clust_dir = "clustering"
  
  # source config.R
  #source(file.path(label,clust_dir, "config.r"))
  clustering.sd.threshold <- 2
  
  rmd = paste0('---
title: "Consensus clustering results for ', type, '"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
---\n
## Overview
This report summarizes the results of the consensus clustering module, which utilizes resampling-based consensus clustering ([Monti et al, 2003](https://doi.org/10.1023/A:1023949509487)) to derive robust proteome clusters. 
1000 bootstrap sample data sets are clustered into *K* clusters using k-means, and a consensus matrix is constructed whose entries (i, j) record the number of times items i and j are assigned to the same cluster divided by the total number of times both items are selected.
A range of possible cluster numbers *K* between 2 and 10 are evaluated and the best *K* is determined by comparing the *empirical cumulative distribution (CDF)* of the resulting consensus matrices. 
To compare the clusterings, the increase of CDF area *K<sub>delta</sub>* is evaluated and the *K* with the largest *K<sub>delta</sub>* is defined as best *K*.
               ')
  
  cm.all = read.csv(file.path(label, clust_dir, paste0(type, "_cluster_metrics.csv"))) %>%
    rename(clust = X)
  cm.best = read.csv(file.path(label, clust_dir, paste0(type, "_best_K.csv")), row.names = 1)
  
  k.max = rownames(cm.all) %>% as.numeric() %>% max()
  k.min = rownames(cm.all) %>% as.numeric() %>% min()
  n.clust = length(rownames(cm.all))
  k.all = rownames(cm.all)
  
  ## number of rows in plot
  #nr <- ceiling((nrow(cm.best)-1)/2)
  
  # recreate metrics plots as png
  png(paste0(type, "_cluster_metrics.png"))
  #par(mfrow=c(nr ,2))
  plots = list()
  loop = colnames(cm.all)[-1]
  count = 1
  for(i in loop){
    print(i)
    #best.idx <- cm.best[i, 'best.k.idx' ]
    best.k.plot = cm.best[i, 'best.k' ]
    metric <- cm.all[, i]
    
    if(i == 'delta.auc.diff'){
      metric <- cm.all[, 'delta.auc']
    }
    print(metric)
    cm.all2 = cm.all %>%
      mutate(highlight = ifelse(cm.all$clust == best.k.plot, "best","not"))
    
    if( !(i %in% c('delta.auc', 'cophenetic.correlation.diff') ) ) {
      p = ggplot(data=cm.all2, aes(x=clust, y = metric, text = metric)) +
        geom_line() +
        geom_point(size = 3, aes(colour = highlight)) +
        scale_color_manual(values=c("red", "black"), guide = FALSE) +
        scale_x_continuous(breaks = cm.all2$clust) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(title = i, x = "Number of clusters", y = "Score") +
        geom_text(aes(label = glue('k={cm.best[i, "best.k"]}'), x = Inf, y = Inf, hjust = 1, vjust = 1))
      print(p)
      plots[count] = p
      count = count + 1
    }
  }


  
  save(cm.best, cm.all, file = "test.Rdata")
  
  rmd = '
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("test.RData")
for (i in 1:length(plots)){
  ggplotly(plots[i], tooltip = "text")
}
```
'

rmd = '
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("test.RData")
for(i in colnames(cm.all)[-1]){
  best.k.plot = cm.best[i, "best.k" ]
  metric <- cm.all[, i]

  if(i == "delta.auc.diff"){
    metric <- cm.all[, "delta.auc"]
  }
  cm.all2 = cm.all %>%
    mutate(highlight = ifelse(cm.all$clust == best.k.plot, "best","not"))
  
  if( !(i %in% c("delta.auc", "cophenetic.correlation.diff") ) ) {
    p = ggplot(data=cm.all2, aes(x=clust, y = metric, text = metric)) +
      geom_line() +
      geom_point(size = 3, aes(colour = highlight)) +
      scale_color_manual(values=c("red", "black"), guide = FALSE) +
      scale_x_continuous(breaks = cm.all2$clust) +
      theme_classic() +
      theme(legend.position = "none") +
      labs(title = i, x = "Number of clusters", y = "Score") +
      geom_text(aes(label = paste0("k= ", cm.best[i, "best.k"]), x = Inf, y = Inf, hjust = 1, vjust = 1))
    print(p)
  }
}
```
'
  
  dev.off()
  
  best_pdf = grep(paste0(type, "_consensus_matrix.*\\.pdf"), list.files(file.path(label, clust_dir)), value = TRUE)
  best_png = gsub("pdf", "png", best_pdf) %>%
    tolower()
  best_k = str_extract(best_png,"[:digit:]") %>% as.numeric()
  
  best_k_data = read.csv(file.path(label, clust_dir, paste0(type, "_best_K.csv"))) %>%
    select(-best.k.idx) %>%
    na.omit() %>%
    rename(algorithm = X)
  
  save(best_k_data, file = "best_k_data.Rdata")
  
  rmd = paste0(rmd, '\n## Consensus matrix for best *K*: ', best_k, '
![**Figure**: Consensus matrix for *K* = 3, determined by several algorithms to be the best *K*.](', file.path(label, clust_dir, best_png), ') \n

**Table**: Best *K* and score for different clustering algorithms.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("best_k_data.RData")
library(DT)
datatable(best_k_data, rownames = FALSE, width = "500px")
```
\n
               ')
  
  ###############################################
  # PCA plot for best K
  data = parse.gctx(file.path(label, clust_dir, paste0(type, "-Cluster.gct")))
  mat = data.frame(data@mat)
  
  membership_full = read.csv(file.path(label, clust_dir, assoc.subgroups), row.names = 1) 
  membership = pull(membership_full, Cluster)
  
  # eliminate features with not enough variation
  min.feat.sd = 500
  feature.sd <- apply (mat, 1, sd, na.rm=TRUE)
  keep <- which( feature.sd > clustering.sd.threshold )
  
  if(length(keep) < min.feat.sd){ # if SD filtering returned too few features for clustering, switch to upper decile mode keeping the top 10 percent
    upper.decile <- quantile(feature.sd, c(.9))
    keep <- which( feature.sd >  upper.decile)
  } 
  mat2 <- mat[keep, ]
  
  p <- fviz_cluster(list(cluster=membership, data=t(na.omit(mat2))), palette='Dark2', text = colnames(mat), geom = "point", labelsize = 8)
  ggplotly(p, tooltip = "text")
  
  p2 <- fviz_cluster(list(cluster=membership, data=t(na.omit(mat2))), palette='Dark2', text = colnames(mat), labelsize = 8)
  ggplotly(p2, tooltip = "text")
  
  #########################################################
  # heatmaps for marker selection in each cluster of best K
  contrast_files = grep(paste0(type, "-Cluster-class\\..*-analysis-markers-fdr", assoc.fdr, "\\.csv"), list.files(file.path(label, clust_dir), full.names = TRUE), value = TRUE)
  yaml_params = read_yaml(yaml_file)
  
  for (contrast in contrast_files){
    ids = read.csv(contrast, stringsAsFactors = FALSE) %>%
      pull(Gene.ID)
    mat_filt = mat[ids,]
    cluster_num = str_extract(contrast,"[:digit:]") %>% as.numeric()
    
    ## heatmap annotations
    if (!is.null (cluster.enrichment.subgroups)) {
      subgroup.table <- read.csv (file.path(label, clust_dir, basename(cluster.enrichment.subgroups)))
      rownames (subgroup.table) <- subgroup.table[,'Sample.ID']
      
      samples <- as.character (subgroup.table[,'Sample.ID'])
      samples <- samples [ samples %in% colnames (mat_filt) ]
      annot <- subgroup.table [samples,]
    } else {
      samples <- colnames (mat_filt)
      annot <- data.frame (Sample.ID=samples)
    }
    rownames (annot) <- samples
    
    annot <- annot %>%
      mutate_at (names (.[, sapply (., function (x) length(unique(x)) <= 5)]), funs (factor(.))) %>%
      mutate (Cluster =  as.character(membership_full[samples, "Cluster"])) %>%
      column_to_rownames ("Sample.ID")
    
    color = lapply(yaml_params$groups.colors, unlist)
    
    annotation <- HeatmapAnnotation (df=annot, annotation_height = 0.25, annotation_width = 0.5,
                                     show_annotation_name=TRUE, col = color)
    
    heatmap <- Heatmap (mat_filt, top_annotation=annotation, 
                        row_names_gp = gpar (fontsize=8), column_names_gp = gpar(fontsize=8), 
                        clustering_distance_columns = 'pearson',
                        clustering_distance_rows='pearson',
                        show_row_names = F,
                        column_title = paste0("Marker selection results for Cluster ", cluster_num),
                        row_title = paste0( dim(mat_filt)[1], " features"),
                        width = 15, height = 10)
    
    png (paste0("marker_hm_clust", cluster_num, ".png"), width = 1500, height = 1000)
    draw (heatmap)
    dev.off()
    
  }
  
  rmd_name = paste(label, type, "cons_clust_rmd.rmd", sep = "_")
  
  # write .rmd file
  writeLines(rmd, con = rmd_name)
  
  # render .rmd file
  rmarkdown::render(rmd_name)

}

rmd_cons_clust(tar.file=tar.file, label=label, tmp.dir=tmp.dir, type=type)

