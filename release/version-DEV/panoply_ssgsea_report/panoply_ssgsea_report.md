# ```panoply_ssgsea_report```

## Description

This module creates an [R Markdown](https://rmarkdown.rstudio.com/) report to provide a high-level summary of the [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) module results.

The report provides:

* A heatmap depicting the normalized enrichment scores (NES) of gene sets significant in at least one data column. Significance is defined by parameter `fdr` in section `panoply_ssgsea_report` of the `cfg_yaml` file.

* List of parameters used in the [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea) module.


## Input

### Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_ssgsea```
* ```cfg_yaml```: (`.yaml` file) master-parameters.yaml
* ```label```: (String) label used in report title

### Optional inputs:

* Heatmap Parameters:
	* ```fdr```: (Float, default=0.05) FDR threshold for heatmap significance stars.
	* ```top_n```: (Int, default=20) max number of significant hits to be plotted on heatmap
	* ```cluster_rows```: (Boolean, default=TRUE) If TRUE, rows will be clustered using distance matrix, with `ser_meth` as seriation method
	* ```ser_meth```: (String, default='ARSA') Seriation method used for clustering heatmap rows.
	* ```split_by_prefix```: (Boolean) If TRUE, separate heatmaps will be created for pathways with unique prefixes '<prefix>-<pathway_name>'. Automatically turned on by default for PTM-SEA pathways, but can be manually overridden with this parameter.
	* ```geneset_groups_file```: (`.csv` file) Optional file used to group related pathways. File must have two columns, with genesets in the first column and geneset-groupings in the second. Geneset names must match the `.gmt` (`gene_set_database`) used in [panoply_ssgsea](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea).

## Output

* ```report```: (`.html` file) [R Markown](https://rmarkdown.rstudio.com/) report.
