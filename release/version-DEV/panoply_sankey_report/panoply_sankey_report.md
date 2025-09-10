# ```panoply_sankey_report```

## Description

This module creates an [R markdown](https://rmarkdown.rstudio.com/) report for the [panoply_sankey](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_sankey) module. The report provides tabset of sankey-diagrams for each sample-annotation file. Each annotation-file's tabset compares the sample-composition of the chosen annotation, to the sample-composition of each other sample-annotation file.

If a primary annotation-file data was provided to ```panoply_sankey```, this annotation file's tabset will be placed at the top of the HTML report.

## Input

### Required inputs:
* ```annot_of_comparison```: (String) annotation value compared in ```panoply_sankey```
* ```sankey_tar```: (`.tar` file) tarball from ```panoply_sankey```
* ```label```: (String) label used in the report title


## Output

* ```report_out```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
