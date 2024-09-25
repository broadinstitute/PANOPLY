# ```panoply_nmf_report```

## Description

This module creates an [R markdown](https://rmarkdown.rstudio.com/) report for the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module. The report provides:

* Table summarizing the number of features per data type used in the clustering.

* Interactive line graph illustrating different cluster metrics as a function of the number of clusters tested. The optimal number of clusters is indicated as vertical line.

* Table summarizing the number of samples in each cluster.

* Heatmap depicting the relative contributions of each sample to each cluster.

* Table summarizing the results of an overrepresentation analysis of sample metadata terms (e.g. clinical annotation, inferred phenotypes, etc.) in each cluster. 

* Heatmap depicting the abundances of cluster specific features across clusters.

* Barpchart depicting the number of cluster specific features.

* Browsable table of all cluster-specific features.

* Heatmap depicting the consensus matrix across multiple random restarts of NMF.

* Silhouette plot depicting cluster stability.

* Table listing all parameters used in [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) and [panoply_nmf_postprocess](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess).

## Input

### Required inputs:
* ```nmf_results```: (`.tar` file) tarball from ```panoply_nmf```
* ```postprocess_tarball```: (`.tar` file) tarball from ```panoply_nmf_postprocess```
* ```nclust```: (Int) Factorization rank used in ```panoply_nmf_postprocess```
* ```label```: (String) label used in the report title


## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
