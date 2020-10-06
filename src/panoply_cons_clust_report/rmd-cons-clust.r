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
library(tidyr)
library(tibble)
library(stringr)
library(glue)
library(ggplot2)
library(plotly)
library(cmapR)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)

rmd_cons_clust <- function(tar_file, yaml_file, label, type){

  # extract files from tarball
  untar(tar_file)
  
  clust_dir = "clustering"
  
  # source config.R
  #source(file.path(label,clust_dir, "config.r"))
  clustering.sd.threshold <- 2
  cluster.enrichment.subgroups <- "/cromwell_root/fc-eaf370c1-9c61-4cd8-b54f-cccd2e316032/sample_sets/all/groups-aggregate.csv"
  assoc.subgroups <- 'proteome-bestclus.csv'
  assoc.fdr <- 0.01
  
  #############################
  # Overview section
  
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
  
  #############################
  # Metrics plots
  
  cm.all = read.csv(file.path(label, clust_dir, paste0(type, "_cluster_metrics.csv"))) %>%
    rename(clust = X) %>%
    select(-c(delta.auc.diff, cophenetic.correlation.diff)) %>%
    rename(delta.auc.diff = delta.auc)
  cm.best = read.csv(file.path(label, clust_dir, paste0(type, "_best_K.csv")), row.names = 1) %>%
    na.omit()
  
  k.min = cm.all$clust %>% min()
  k.max = cm.all$clust %>% max() 
  
  # recreate metrics plots 
  cm.all2 = cm.all %>% 
    gather(key = "algorithm", value = "value", -clust) %>%
    mutate(highlight = "none") %>%
    mutate(value = signif(value, 3))
  
  for (i in rownames(cm.best)){
    best.k.plot = cm.best[i, 'best.k' ]
    index = which(cm.all2$clust == best.k.plot & cm.all2$algorithm == i)
    cm.all2[index, "highlight"] = "best"
  }
  
  metrics_plot = ggplot(data=cm.all2, aes(x=clust, y = value, text = value)) +
    geom_line() +
    geom_point(size = 3, aes(colour = highlight)) +
    facet_wrap(~algorithm, ncol = 2) +
    scale_color_manual(values=c("red", "black"), guide = FALSE) +
    scale_x_continuous(breaks = cm.all2$clust) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = "Clustering metrics", x = "Number of clusters", y = "Score") +
    geom_text(aes(label = glue('k={cm.best[i, "best.k"]}'), x = Inf, y = Inf, hjust = 1, vjust = 1))
  
  best_k_data = cm.best %>%
    rownames_to_column("algorithm") %>%
    select(-best.k.idx) %>%
    mutate(best.k.score = signif(best.k.score, 3))
  
  save(metrics_plot, best_k_data, file = "metrics_plot.Rdata")
  
  rmd = paste0(rmd, '\n## Metrics
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("metrics_plot.RData")
ggplotly(metrics_plot, tooltip = "text")
```
**Figure**: Interactive plots of clustering metrics for *K* = ', k.min, ' through *K* = ', k.max, '. The best *K* is highlighted in red. Hover over each point to see the score value it corresponds to.


**Table**: Best *K* and score for different clustering algorithms.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("metrics_plot.RData")
library(DT)
datatable(best_k_data, rownames = FALSE, width = "500px")
```
\n
              ')
  
  #######################
  # Consensus matrix for best K
  
  best_pdf = grep(paste0(type, "_consensus_matrix.*\\.pdf"), list.files(file.path(label, clust_dir)), value = TRUE)
  best_png = gsub("pdf", "png", best_pdf) %>%
    tolower()
  best_k = str_extract(best_png,"[:digit:]") %>% as.numeric()
  
  rmd = paste0(rmd, '\n## Best *K* = ', best_k, '\n
### Consensus matrix for *K* = ', best_k, '\n
![**Figure**: Consensus matrix for *K* = 3, determined by several algorithms to be the best *K*.](', file.path(label, clust_dir, best_png), ') \n
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
  
  pca <- fviz_cluster(list(cluster=membership, data=t(na.omit(mat2))), palette='Dark2', text = colnames(mat), geom = "point", labelsize = 8)
  
  pca2 <- fviz_cluster(list(cluster=membership, data=t(na.omit(mat2))), palette='Dark2', text = colnames(mat), labelsize = 8)

  save(pca, pca2, file = "pca_data.Rdata")
  
  rmd = paste0(rmd, '\n### PCA plot for best *K* = ', best_k, '\n
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("pca_data.RData")
ggplotly(pca, tooltip = "text")
ggplotly(pca2, tooltip = "text")
```
**Figure**: PCA plot for *K* = 3. Clusters are separated by colors. Hover over each point to see the sample name it corresponds to.
               ')
  
  #########################################################
  # heatmaps for marker selection in each cluster of best K
  # contrast_files = grep(paste0(type, "-Cluster-class\\..*-analysis-markers-fdr", assoc.fdr, "\\.csv"), list.files(file.path(label, clust_dir), full.names = TRUE), value = TRUE)
  yaml_params = read_yaml(yaml_file)
  full_MS_file = file.path(label, clust_dir, paste0(type, "-Cluster-analysis-markers-fdr", assoc.fdr, ".csv"))
  ids = read.csv(full_MS_file, stringsAsFactors = FALSE) %>%
        pull(Gene.ID)
  mat_filt = mat[ids,]
  
  # for (contrast in contrast_files){
  #   ids = read.csv(contrast, stringsAsFactors = FALSE) %>%
  #     pull(Gene.ID)
  #   mat_filt = mat[ids,]
  #   cluster_num = str_extract(contrast,"[:digit:]") %>% as.numeric()
    
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
    arrange(Cluster) %>%
    column_to_rownames ("Sample.ID")
  
  column_order = rownames(annot)
  
  color = lapply(yaml_params$groups.colors, unlist)
  clust.col = brewer.pal(n = best_k, "Dark2")
  names(clust.col) = 1:best_k
  color["Cluster"] = list(clust.col)
  
  mat_filt2 = mat_filt[,column_order]
  
  annotation <- HeatmapAnnotation (df=annot, annotation_height = 0.25, annotation_width = 0.5,
                                   show_annotation_name=TRUE, col = color)
  
  heatmap <- Heatmap (as.matrix(mat_filt2), top_annotation=annotation, 
                      row_names_gp = gpar (fontsize=8), column_names_gp = gpar(fontsize=8), 
                      cluster_columns = FALSE,
                      clustering_distance_rows='pearson',
                      show_row_names = F,
                      column_title = paste0("Marker selection results, ", dim(mat_filt)[1], " features"),
                      width = 15, height = 10)
    
  png ("marker_hm_clust.png", width = 1000, height = 1500)
  draw (heatmap)
  dev.off()
  
  rmd = paste0(rmd, '\n### Marker selection results for *K* = ', best_k, '\n
![**Figure**: Heatmap showing marker selection results for *K* = 3 clusters. Columns (samples) are sorted by cluster number and rows are clustered by hierarchical clustering using the Pearson correlation method. Sample annotations, including cluster number, are labeled at the top of the heatmap.](marker_hm_clust.png)
               ')
  
  rmd_name = paste(label, type, "cons_clust_rmd.rmd", sep = "_")
  
  # write .rmd file
  writeLines(rmd, con = rmd_name)
  
  # render .rmd file
  rmarkdown::render(rmd_name)

}

rmd_cons_clust(tar_file = tar_file, yaml_file = yaml_file, label=label, type=type)

