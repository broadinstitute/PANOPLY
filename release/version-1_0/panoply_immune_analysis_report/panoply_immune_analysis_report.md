# ```panoply_immune_analysis_report```

## Description

This module creates an interactive  [R Markdown](https://rmarkdown.rstudio.com/) report of the [panoply_immune_analysis](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_immune_analysis) module. xCell scores are shown in an annotated heatmap and compared to ESTIMATE scores in interactive scatter plots. Significant (defined by FDR value given in ```panoply_immune_analysis```) results, if any, of enrichment analysis using Fisher\'s exact test between immune subtypes and annotation groups are shown in an interactive table.

## Input

### Required inputs:

* ```tar_file```: (`.tar` file) output tar file containing results from ```panoply_immune_analysis```
* ```yaml_file```: (`.yaml` file) master parameters file
* ```label```: (String) directory name within the tar file

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.