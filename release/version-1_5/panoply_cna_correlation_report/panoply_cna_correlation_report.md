# ```panoply_cna_correlation_report```

## Description

This module creates an [R Markdown](https://rmarkdown.rstudio.com/) report for the [panoply_cna_correlation](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_cna_correlation) module.

The report provides:

* A genome-wide heatmap depicting significant (FDR < `fdr`) _cis_ and _trans_ correlation between CNA-RNA and CNA-`type`.

* A table providing a high level summary about the total number of genes that were used the analysis, number of significant genes (FDR < `fdr`), number of significant events (_cis_ + _trans_) and the number of significant _cis_ effects. 

## Input

### Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_cna_correlation```
* ```config_yaml```: (`.yaml` file) master-parameters.yaml
* ```label```: (String) label used in the report title
* ```type```: (String) data (-ome) type
* ```tmpDir```: (String, default="tmp") name of directory used to extract the tarball


### yaml-exclusive parameters:

`panoply_cna_correlation_report`

* ```fdr```: (Float, default = 0.05) FDR value used to define signifcance.

## Output

* ```report```: (`.html` file) [R Markown](https://rmarkdown.rstudio.com/) report.
