# ```panoply_so_nmf_sankey_report```

## Description

This module creates an [R markdown](https://rmarkdown.rstudio.com/) report for the [panoply_so_nmf_sankey](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_so_nmf_sankey) module. The report provides tabset of sankey-diagrams for each -ome. Each -ome's tabset compares the sample-composition of clusters, produced with those single-omic features, to the sample-composition of clusters produced by each of the other single-omes.

If multi-omics data was provided to ```panoply_so_nmf_sankey```, an additional tabset of sankey-diagrams will be created, comparing the sample composition of clusters produced using multi-omic features, to those produced by each single-ome.

## Input

### Required inputs:
* ```tarball```: (`.tar` file) tarball from ```panoply_so_nmf_sankey```
* ```label```: (String) label used in the report title


## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
