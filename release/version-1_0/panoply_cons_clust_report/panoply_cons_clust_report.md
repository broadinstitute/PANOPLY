# ```panoply_cons_clust_report```

## Description

This module creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the [panoply_cons_clust](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_cons_clust) module. Clustering metrics used to define the best *K* are plotted, and the consensus matrix and principal component analysis for the best *K* are displayed. Marker selection results for the best *K* are shown in a heatmap, and interactive volcano plots show the GSEA results for each cluster with significant pathways, if any, additionally described in an interactive table.

## Input

### Required inputs:

* ```tar_file```: (`.tar` file) output tar file containing results from ```panoply_cons_clust```
* ```yaml_file```: (`.yaml` file) master parameters file
* ```label```: (String) directory name within the tar file
* ```type```: (String) data (-ome) type

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.