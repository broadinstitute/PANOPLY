# ```panoply_association_report```

## Description

This module creates an interactive [R Markown](https://rmarkdown.rstudio.com/) report of the GSEA results of the [panoply_association](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_association) module. For each contrast, an interactive volcano plot shows the Normalized Enrichment Score (NES) vs -log10 of the FDR value for each database pathway analyzed. Hovering over a point reveals which pathway it corresponds to. Any significantly enriched pathways (determined by the ```fdr_value``` input) are highlighted in red and described in an interactive, searchable table. For multi-contrast categories, an overview heatmap shows NES values for all contrasts, with asterisks denoting signficant enrichment.

## Input

## Required inputs:

* ```input_tar```: (`.tar` file) output tar file containing results from ```panoply_association```
* ```master_yaml```: (`.yaml` file) master parameters file
* ```label```: (String) directory name within the tar file
* ```type```: (String) data (-ome) type

### Optional inputs:

* ```fdr_value```: (Float, default = 0.01) FDR value cutoff to be considered significant.

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
