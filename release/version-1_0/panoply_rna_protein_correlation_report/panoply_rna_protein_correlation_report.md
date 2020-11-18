Documentation at https://github.com/broadinstitute/PANOPLY/blob/version-1_0/release/version-1_0/panoply_rna_protein_correlation_report/panoply_rna_protein_correlation_report.md

# ```panoply_rna_protein_correlation_report```

## Description

This module creates an interactive [R markdown](https://rmarkdown.rstudio.com/) report for the [panoply_rna_protein_correlation](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_rna_protein_correlation). The report provides:

* Interactive histogram visualizing the distribution of gene-level RNA-protein correlation.

* High-level summary table about the number total number and number of significant correlations (FDR < `$fdr`). 

* Interactive correlation rank plot depicting gene-level correlations ordered by rank.

* Table summarizing the top 10 most positively and negatively correlated protein/RNA pairs.

* Bubble chart showing the sample-level RNA-protein correlations.


## Input

### Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_rna_protein_correlation```
* ```config_yaml```: (`.yaml` file) master-parameters.yaml
* ```label```: (String) label used in the report title
* ```type```: (String) data (-ome) type
* ```tmpDir```: (String, default="tmp") name of directory used to extract the tarball

### yaml-exclusive parameters:

`panoply_rna_protein_correlation_report`

* ```fdr```: (Float, default = 0.05) FDR value used to define signifcance.

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.