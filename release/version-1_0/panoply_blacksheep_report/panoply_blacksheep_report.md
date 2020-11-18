Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_blacksheep_report/panoply_blacksheep_report.md

# ```panoply_blacksheep_report```

## Description

This module creates an interactive [R Markown](https://rmarkdown.rstudio.com/) report of the significant results of the [panoply_blacksheep](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_blacksheep) module, if a groups file was provided for outlier analysis. Signficance is determined by the ```fdr_value``` input to [panoply_blacksheep](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_blacksheep). For each comparison group, the number of genes significantly enriched for outliers is listed. If there are one or more significant genes, an interactive, searchable table showing the gene names and their corresponding FDR values is displayed, as well as the heatmap showing the fraction of outliers in these genes for visualization of the enrichment of outliers in the group of interest.

## Input

### Required inputs:

* ```input_tar```: (`.tar` file) output tar file from ```panoply_blacksheep```
* ```output_prefix```: (String) prefix used to name the output tar file in ```panoply_blacksheep```
* ```type```: (String) data (-ome) type

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.