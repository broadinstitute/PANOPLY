Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_sampleqc_report/panoply_sampleqc_report.md

# ```panoply_sampleqc_report```

## Description

This module creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report for the [panoply_sampleqc](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_sampleqc) module. 

The report provides:

* Interactive box-whisker plots visualizing the distributions of stromal, immune and tumor purity scores inferred by  [ESTIMATE](https://bioinformatics.mdanderson.org/public-software/estimate/) across data types.

* Interactive heatmaps visualizing pairwise sample correlations between RNA, CNA and proteome data to expore potential sample swaps. Heatmaps are created with the [MORPHEUS](https://software.broadinstitute.org/morpheus/) [R-package](https://github.com/cmap/morpheus.R). 


## Input

### Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_sampleqc```
* ```label```: (String) label
* ```type```: (String) data (-ome) type
* ```tmpDir```: (String, default="tmp") name of directory used to extract the tarball

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.