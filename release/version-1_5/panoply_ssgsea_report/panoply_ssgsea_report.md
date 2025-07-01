# ```panoply_ssgsea_report```

## Description

This module creates an [R Markdown](https://rmarkdown.rstudio.com/) report to provide a high-level summary of the [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) module results.

The report provides:

* A heatmap depicting the normalized enrichment scores (NES) of gene sets significant in at least one data column. Significance is defined by parameter `fdr` in section `panoply_ssgsea_report` of the `cfg_yaml` file.
.
* List of parameters used in the [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) module.


## Input

### Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_ssgsea```
* ```cfg_yaml```: (`.yaml` file) master-parameters.yaml
* ```label```: (String) label used in report title

## Output

* ```report```: (`.html` file) [R Markown](https://rmarkdown.rstudio.com/) report.
