#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for consensus clustering module ###

# tar_file  - URL of tar file created by task panoply_cons_clust
# yaml_file - URL of master parameters yaml file including output from startup notebook
# label     - character, name of folder in tarball
# data_type - character, data type

args = commandArgs(TRUE)

tar_file = args[1]
yaml_file = args[2]
label = args[3]
data_type = args[4]

library(rmarkdown)
library(yaml)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(plotly)
library(cmapR)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)

hallmark_process_category <- c(
  HALLMARK_TNFA_SIGNALING_VIA_NFKB='signaling',
  HALLMARK_HYPOXIA='pathway',
  HALLMARK_CHOLESTEROL_HOMEOSTASIS='metabolic',
  HALLMARK_MITOTIC_SPINDLE='proliferation',
  HALLMARK_WNT_BETA_CATENIN_SIGNALING='signaling',
  HALLMARK_TGF_BETA_SIGNALING='signaling',
  HALLMARK_IL6_JAK_STAT3_SIGNALING='immune',
  HALLMARK_DNA_REPAIR='DNA damage',
  HALLMARK_G2M_CHECKPOINT='proliferation',
  HALLMARK_APOPTOSIS='pathway',
  HALLMARK_NOTCH_SIGNALING='signaling',
  HALLMARK_ADIPOGENESIS='development',
  HALLMARK_ESTROGEN_RESPONSE_EARLY='signaling',
  HALLMARK_ESTROGEN_RESPONSE_LATE='signaling',
  HALLMARK_ANDROGEN_RESPONSE='signaling',
  HALLMARK_MYOGENESIS='development',
  HALLMARK_PROTEIN_SECRETION='pathway',
  HALLMARK_INTERFERON_ALPHA_RESPONSE='immune',
  HALLMARK_INTERFERON_GAMMA_RESPONSE='immune',
  HALLMARK_APICAL_JUNCTION='cellular component',
  HALLMARK_APICAL_SURFACE='cellular component',
  HALLMARK_HEDGEHOG_SIGNALING='signaling',
  HALLMARK_COMPLEMENT='immune',
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE='pathway',
  HALLMARK_PI3K_AKT_MTOR_SIGNALING='signaling',
  HALLMARK_MTORC1_SIGNALING='signaling',
  HALLMARK_E2F_TARGETS='proliferation',
  HALLMARK_MYC_TARGETS_V1='proliferation',
  HALLMARK_MYC_TARGETS_V2='proliferation',
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION='development',
  HALLMARK_INFLAMMATORY_RESPONSE='immune',
  HALLMARK_XENOBIOTIC_METABOLISM='metabolic',
  HALLMARK_FATTY_ACID_METABOLISM='metabolic',
  HALLMARK_OXIDATIVE_PHOSPHORYLATION='metabolic',
  HALLMARK_GLYCOLYSIS='metabolic',
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY='pathway',
  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY='pathway',
  HALLMARK_P53_PATHWAY='proliferation',
  HALLMARK_UV_RESPONSE_UP='DNA damage',
  HALLMARK_UV_RESPONSE_DN='DNA damage',
  HALLMARK_ANGIOGENESIS='development',
  HALLMARK_HEME_METABOLISM='metabolic',
  HALLMARK_COAGULATION='immune',
  HALLMARK_IL2_STAT5_SIGNALING='signaling',
  HALLMARK_BILE_ACID_METABOLISM='metabolic',
  HALLMARK_PEROXISOME='cellular component',
  HALLMARK_ALLOGRAFT_REJECTION='immune',
  HALLMARK_SPERMATOGENESIS='development',
  HALLMARK_KRAS_SIGNALING_UP='signaling',
  HALLMARK_KRAS_SIGNALING_DN='signaling',
  HALLMARK_PANCREAS_BETA_CELLS='development'
)
hallmark_category = as.data.frame(hallmark_process_category)
names(hallmark_category) = "category"
hallmark_category = hallmark_category %>%
  rownames_to_column("pathway")

rmd_cons_clust <- function(tar_file, yaml_file, label, data_type){

  # extract files from tarball
  untar(tar_file)
  
  clust_dir = "clustering"
  
  # source config.R
  source(file.path(label,clust_dir, "config.r"))
  
  #############################
  # Overview section
  
  rmd = paste0('---
title: "Consensus clustering results for ', data_type, '"
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
Detailed documentation for the consensus clustering module can be found [here](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_cons_clust).
This report shows metrics comparing the clusterings used to determine the best *K*, the consensus matrix for the best *K*, principal component analysis for the best *K* clusters, and marker selection & GSEA results for each cluster.
               ')
  
  #############################
  # Metrics plots
  
  cm.all = read.csv(file.path(label, clust_dir, paste0(data_type, "_cluster_metrics.csv"))) %>%
    rename(clust = X) %>%
    select(-c(delta.auc.diff, cophenetic.correlation.diff)) %>%
    rename(delta.auc.diff = delta.auc)
  cm.best = read.csv(file.path(label, clust_dir, paste0(data_type, "_best_K.csv")), row.names = 1) %>%
    na.omit()
  
  k.min = cm.all$clust %>% min()
  k.max = cm.all$clust %>% max() 
  
  # metrics plots 
  cm.all2 = cm.all %>% 
    gather(key = "metric", value = "value", -clust) %>%
    mutate(highlight = "none") %>%
    mutate(value = signif(value, 3))

  for (i in rownames(cm.best)){
    best.k.plot = cm.best[i, 'best.k' ]
    index = which(cm.all2$clust == best.k.plot & cm.all2$metric == i)
    cm.all2[index, "highlight"] = "best"
    cm.all2$metric[cm.all2$metric == i] = paste0(i, ", K = ", best.k.plot)
  }

  metrics_plot = ggplot(data=cm.all2, aes(x=clust, y = value, text = value)) +
    geom_line() +
    geom_point(size = 3, aes(colour = highlight)) +
    facet_wrap(~metric) +
    scale_color_manual(values=c("red", "black"), guide = FALSE) +
    scale_x_continuous(breaks = cm.all2$clust) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = "Plots of clustering metrics", x = "\r\nNumber of clusters", y = "Score\n\r") 
  
  metrics_plotly = ggplotly(metrics_plot, tooltip = "text", width = 800, height = 400)
  metrics_plotly$x$layout$margin$b = metrics_plotly$x$layout$margin$b + 5

  best_k_data = cm.best %>%
    rownames_to_column("metric") %>%
    select(-best.k.idx) %>%
    mutate(best.k.score = signif(best.k.score, 3))
  
  save(metrics_plotly, best_k_data, file = "metrics_plot.RData")
  
  rmd = paste0(rmd, '\n## Clustering metrics
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("metrics_plot.RData")
metrics_plotly
```
**Figure**: Interactive plots of clustering metrics for *K* = ', k.min, ' through *K* = ', k.max, '. The best *K* is highlighted in red for each metric. Hover over each point to see the score value it corresponds to.


**Table**: Best *K* and score for different clustering metrics.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("metrics_plot.RData")
library(DT)
datatable(best_k_data, rownames = FALSE, width = "500px")
```
\n
              ')
  
  #######################
  # Consensus matrix for best K
  
  best_pdf = grep(paste0(data_type, "_consensus_matrix.*\\.pdf"), list.files(file.path(label, clust_dir)), value = TRUE)
  best_png = gsub("pdf", "png", best_pdf) %>%
    tolower()
  best_k = str_extract(best_png,"[:digit:]") %>% as.numeric()
  
  rmd = paste0(rmd, '\n## Results for best *K* = ', best_k, '\n
### Consensus matrix for *K* = ', best_k, '\n
![**Figure**: Consensus matrix for best *K* = 3.](', file.path(label, clust_dir, best_png), ') \n
               ')
  
  ###############################################
  # PCA plot for best K
  data = parse.gctx(file.path(label, clust_dir, paste0(data_type, "-Cluster.gct")))
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
  
  pca <- fviz_cluster(list(cluster=membership, data=t(na.omit(mat2))), palette='Dark2', repel=TRUE, labelsize = 8, main = paste("PCA, K =", best_k))
  png ("pca_plot.png", width = 800, height = 800)
  print(pca)
  dev.off()
  
  rmd = paste0(rmd, '\n### PCA plot for best *K* = ', best_k, '\n
![**Figure**: Plot of principal component analysis for *K* = 3. Clusters are separated by colors and shapes. Sample names are labeled for each point.](pca_plot.png)
               ')
  
  #########################################################
  # heatmaps for marker selection in each cluster of best K
  yaml_params = read_yaml(yaml_file)
  full_MS_file = file.path(label, clust_dir, paste0(data_type, "-Cluster-analysis-markers-fdr", assoc.fdr, ".csv"))
  ids = read.csv(full_MS_file, stringsAsFactors = FALSE) %>%
        pull(Gene.ID)
  mat_filt = mat[ids,]
    
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
  clust.col = brewer.pal(n = best_k, "Set3")
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
![**Figure**: Heatmap showing marker selection results for *K* = 3 clusters. Columns (samples) are sorted by cluster number and rows (features) are clustered by hierarchical clustering using the Pearson correlation method. Sample annotations, including cluster number, are labeled at the top of the heatmap.](marker_hm_clust.png)
               ')
  
  #############################
  # Hallmark categories for each cluster
  
  gsea_dirs = grep(paste0(data_type, "-Cluster-class\\..*-analysis-gsea-analysis"), list.dirs(file.path(label, clust_dir), full.names = FALSE), value = TRUE)
  gsea_pos = "gsea.SUMMARY.RESULTS.REPORT.0.txt"
  gsea_neg = "gsea.SUMMARY.RESULTS.REPORT.1.txt"
  
  nperm = 1000
  
  for (dir_name in gsea_dirs){
    cluster_num = str_extract(dir_name,"[:digit:]") %>% as.numeric()
    hallmark_neg = read.delim(file.path(label, clust_dir, dir_name, gsea_neg))
    hallmark = read.delim(file.path(label, clust_dir, dir_name, gsea_pos)) %>%
      rbind(hallmark_neg) %>%
      select(GS, NES, FDR.q.val) %>%
      rename(pathway = GS, fdr = FDR.q.val) %>%
      mutate(pathway = as.character(pathway)) %>%
      left_join(hallmark_category) 
    
    hallmark_volc = hallmark %>%
      mutate(group = ifelse(fdr<assoc.fdr, "Significant", "Not significant")) %>%
      mutate(pathway = gsub("^HALLMARK_", "", pathway)) %>%
      mutate(fdr = ifelse(fdr == 0, 1/nperm, fdr))
    
    hallmark_table = hallmark %>%
      mutate(fdr = signif(fdr, 3),
             NES = signif(NES, 3)) %>%
      filter(fdr < assoc.fdr) %>%
      arrange(fdr)
    
    save(hallmark_volc, hallmark_table, file = paste0("hallmark", cluster_num, ".RData"))
    
    rmd = paste0(rmd, '\n#### **GSEA results for Cluster ', cluster_num, '**
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("hallmark', cluster_num, '.RData")
p = ggplot(hallmark_volc, aes(x = NES, y = -log10(fdr), colour = group, text = pathway, group = category)) +
  geom_point() +
  geom_hline(yintercept = -log10(', assoc.fdr, '), linetype = "dashed") +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "NES.Rest.over.Cluster', cluster_num, '")
ggplotly(p, tooltip = c("text", "category"))
```
**Figure**: Interactive volcano plot summarizing GSEA pathway results for cluster ', cluster_num, '. X axis represents the Normalized Enrichment Score (NES); negative NES values indicate enrichment in cluster ', cluster_num, ', and positive NES values indicate enrichment in the rest of the clusters. Y axis represents the -log10 of the FDR value. Results above the dashed line are significant at FDR cutoff = ', assoc.fdr , '. Hover over each point to see which pathway and category it corresponds to.
                 ')
    
    if (dim(hallmark_table)[1]>= 1){
      rmd = paste0(rmd, '\n
**Table**: ', dim(hallmark_table)[1], ' significantly enriched pathways in cluster ', cluster_num, ' with FDR < ', assoc.fdr, '.
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("hallmark', cluster_num, '.RData")
library(DT)
datatable(hallmark_table, rownames = FALSE, width = "500px")
```
\n
                 ')
    } else {
      rmd = paste0(rmd, '\n
No significantly enriched pathways in cluster ', cluster_num, ' with FDR < ', assoc.fdr, '.  
\n
                 ')
    }
      
  }
  
  rmd_name = paste(label, data_type, "cons_clust_rmd.rmd", sep = "_")
  
  # write .rmd file
  writeLines(rmd, con = rmd_name)
  
  # render .rmd file
  rmarkdown::render(rmd_name)

}

rmd_cons_clust(tar_file = tar_file, yaml_file = yaml_file, label = label, data_type = data_type)

