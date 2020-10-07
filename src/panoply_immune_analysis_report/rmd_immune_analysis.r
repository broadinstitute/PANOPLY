#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
### Create Rmarkdown report for immune analysis module ###

# tar_file  - URL of tar file created by task panoply_immune_analysis
# yaml_file - URL of master parameters yaml file including output from startup notebook
# label     - character, name of folder in tarball

args = commandArgs(TRUE)

tar_file = args[1]
yaml_file = args[2]
label = args[3]

library(rmarkdown)
library(yaml)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(tidyr)
library(ggplot2)
library(plotly)
library(RColorBrewer)

rmd_immune = function(tar_file, yaml_file, label){
  
  # extract files from tarball
  untar(tar_file)
  
  # read yaml file
  yaml_params = read_yaml(yaml_file)
  
  # source config.R
  source(file.path(label, "config.r"))
  
  immune_dir = "immune-analysis"
  
  # read in all immune results
  rna.ES = read.csv(file.path(label, immune_dir, "estimate-scores.csv"), row.names = 1)
  rna.XC = read.csv(file.path(label, immune_dir, "xcell-scores.csv"), row.names = 1)
  rna.subtype = read.csv(file.path(label, immune_dir, "immune-subtype.csv"), row.names = 1)
  
  # make xCell heatmap (including annotations from ESTIMATE & immune subtyping)
  results <- t (rna.XC [, -(1:3 + ncol (rna.XC))])

  # assign annotations depending on groups file  
  if (!is.null (immune.enrichment.subgroups)) {
    subgroup.table <- read.csv (file.path(label, immune_dir, basename(immune.enrichment.subgroups)))
    rownames (subgroup.table) <- subgroup.table[,'Sample.ID']
    
    samples <- as.character (subgroup.table[,'Sample.ID'])
    samples <- samples [ samples %in% colnames (results) ]
    annot <- subgroup.table [samples,]
  } else {
    samples <- colnames (results)
    annot <- data.frame (Sample.ID=samples)
  }
  rownames (annot) <- samples
  
  annot <- annot %>%
    mutate_at (names (.[, sapply (., function (x) length(unique(x)) <= 5)]), funs (factor(.))) %>%
    mutate (xCell.ImmuneScore = rna.XC[samples, 'ImmuneScore']) %>%
    mutate (xCell.StromaScore = rna.XC[samples, 'StromaScore']) %>% 
    mutate (xCell.MicroenvironmentScore = rna.XC[samples, 'MicroenvironmentScore']) %>%
    mutate (ESTIMATE.TumorPurity = rna.ES[samples, 'TumorPurity']) %>%
    mutate (Immune.Subtype = factor (rna.subtype[samples, 'Immune.Subtype'])) %>%
    select (-Sample.ID)
  
  color = lapply(yaml_params$groups.colors, unlist)
  subtype_col = brewer.pal(n=6, "Set3")
  names(subtype_col) = 1:6
  color["Immune.Subtype"] = list(subtype_col)
  
  annotation <- HeatmapAnnotation (df=annot, annotation_height = 0.5, annotation_width = 0.5,
                                   show_annotation_name=TRUE, col = color)
  
  heatmap <- Heatmap (results, top_annotation=annotation, 
                      row_names_gp = gpar (fontsize=8), column_names_gp = gpar(fontsize=8),
                      clustering_distance_columns='pearson', clustering_distance_rows='pearson',
                      width=immune.heatmap.width, height=immune.heatmap.height)
  
  png ('xcell-scores-heatmap.png', 
       height=1800, width=1200)
  draw (heatmap)
  dev.off()
  
  # xCell vs ESTIMATE plots
  plot.data <- data.frame (id=samples, xCell.ImmuneScore=rna.XC[samples,'ImmuneScore'], 
                           xCell.StromalScore=rna.XC[samples,'StromaScore'],
                           ESTIMATE.ImmuneScore=rna.ES[samples,'ImmuneScore'], 
                           ESTIMATE.StromalScore=rna.ES[samples,'StromalScore'])
  
  scatter.data <- gather (plot.data, key, value, -id) %>%
    separate (key, into=c('type', 'score'), sep='\\.') %>%
    spread (type, value)
  save(scatter.data, file = "XC_ES_scatterdata.Rdata")
  
  # read immune subtype enrichment results filtered for pval cutoff
  subtype_data = read.csv(file.path(label, immune_dir, paste0("immune-subtype-enrichment-pval", immune.enrichment.fdr, ".csv"))) %>%
    mutate(fisher.test.pvalue = signif(fisher.test.pvalue, 3), 
           adj.pvalue = signif(adj.pvalue, 3))
  
  subtype = data.frame(Immune.Subtype = rep(1:6), 
                             Immune.Subtype.Description = c("Wound healing",
                                                            "IFN-gamma dominant",
                                                            "Inflammatory",
                                                            "Lymphocyte depleted",
                                                            "Immunologically quiet",
                                                            "TGF-beta dominant")) %>%
    right_join(subtype_data)

  
  # write rmd
  rmd = paste0('---
title: "Immune analysis results"
output: 
  html_document:
   toc: true
   toc_depth: 3
   toc_float:
     collapsed: true
---\n
## Overview
This report summarizes the results of the immune analysis module, which runs several algorithms assigning immune scores for understanding the tumor microenvironment. **E**stimation of **ST**romal and **I**mmune cells in **MA**lignant **T**umor tissues using **E**xpression data (**ESTIMATE**, [Yoshihara et al., 2013](https://doi.org/10.1038/ncomms3612)) uses gene expression signatures to calculate the fraction of stromal and immune cells and infer tumor purity. **xCell** ([Aran et al., 2017](https://doi.org/10.1186/s13059-017-1349-1)) uses a gene signature-based method to infer immune and stromal cell types. **ImmuneSubtypeClassifier** ([Thorrson et al., 2018](https://doi.org/10.1016/j.immuni.2018.03.023)) uses immune gene expression signatures to classify tumor samples into one of 6 immune subtypes. Enrichment analysis (Fisher\'s exact test) is performed on immune subtypes and significant results at p-value < ', immune.enrichment.fdr, ' reported.

## xCell scores

![**Figure**: Heatmap showing xCell scores for each sample. Rows (cell types) and columns (samples) are clustered by hierarchical clustering using the Pearson correlation method. Sample annotations, including key immune scores from all three algorithms, are labeled at the top of the heatmap.](xcell-scores-heatmap.png)

## Comparison of xCell and ESTIMATE scores
```{r echo=FALSE, warning=FALSE, message=FALSE}
load("XC_ES_scatterdata.Rdata")
p = ggplot (aes (x=ESTIMATE, y=xCell, group=score, color=score, text = id), data=scatter.data) + geom_point() +
  facet_wrap(~score, scales="free")
ggplotly(p, tooltip = "text", height = 500, width = 900)
```
**Figure**: Interactive scatter plots comparing Immune and Stromal Scores calculated by xCell and ESTIMATE. Hovering over each point reveals which sample it corresponds to.


## Immune subtype enrichment analysis

Number of significant enrichments between immune subtypes and annotation groups: ', dim(subtype)[1]
               )
  
  if (dim(subtype)[1] >= 1){
    save(subtype, file = "subtype_enrichment.Rdata")
    rmd = paste0(rmd, '\n
**Table**: significant (p-value < ', immune.enrichment.fdr, ') enrichment analyses between immune subtypes and annotation groups.\n
* Subtype 1 = Wound healing
* Subtype 2 = IFN-gamma dominant
* Subtype 3 = Inflammatory
* Subtype 4 = Lymphocyte depleted
* Subtype 5 = Immunologically quiet
* Subtype 6 = TGF-beta dominant \n

```{r echo=FALSE, warning=FALSE, message=FALSE}
load("subtype_enrichment.Rdata")
library(DT)
datatable(subtype, rownames = FALSE, width = "500px")
```               
           ')
  }
  
  rmd_name = paste(label, "immune_rmd.rmd", sep = "_")
  
  # write .rmd file
  writeLines(rmd, con = rmd_name)
  
  # render .rmd file
  rmarkdown::render(rmd_name)
  
}

# run rmd_association function to make rmd report
rmd_immune(tar_file = tar_file, yaml_file = yaml_file, label = label)
