# ```panoply_clumps_ptm_report```

## Description

This module creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the clumps_ptm results of the [panoply_clumps_ptm](./Data-Analysis-Modules%3A-panoply_clumps_ptm) module. The report contains summary figures listing the top significant proteins for each test, as well as PyMol figures showing the locations of significant PTM-sites on protein structures.

This module can take multiple [panoply_clumps_ptm_postprocess](./Data-Analysis-Modules%3A-panoply_clumps_ptm_postprocess) results files, displaying results from multiple annotations in a single interactive report.

## Input

### Required inputs:

* ```postprocess_results```: (`.tar` file(s)) output tar file(s) containing results from ```panoply_clumps_ptm_postprocess```
* ```label```: (String) prefix used to name the output files

### Optional inputs:

* ```mapping_params```: (`.yaml` file) optional file from [panoply_clumps_ptm_mapping](./Data-Preparation-Modules%3A-panoply_clumps_ptm_mapping), containing parameters used for PTM-mapping; allows mapping parameters to be listed in the report. Automatically provided if `panoply_clumps_ptm_mapping` was run in the workflow; must be manually provided if mapping was skipped/overridden.


## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
