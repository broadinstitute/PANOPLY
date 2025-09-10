# ```panoply_metaboanalyst_report```

## Description

This module creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the metaboanalyst results of the [panoply_metaboanalyst](./Data-Analysis-Modules%3A-panoply_metaboanalyst) module. It contains the summary heatmaps produced in `panoply_metaboanalyst`, as well as interactive network viewers for all pathways, colorcoded according to significant features.

This report is designed for the joint-pathway analysis version of `panoply_metaboanalyst`, and assumes that Overrepresentation Analysis (`anal_type="ORA"`) has been run. If quantitative enrichment analysis or metabolite-only analysis has been run instead, figures may not behave as expected.

## Input

### Required inputs:

* ```metaboanalyst_results```: (`.tar` file) output tar file(s) containing results from ```panoply_metaboanalyst```
* ```label```: (String) prefix used to name the `.html` report file


## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
