Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_normalize_ms_data_report/panoply_normalize_ms_data_report.md

# ```panoply_normalize_ms_data_report```

## Description

This module creates an [R markdown](https://rmarkdown.rstudio.com/) report for the [panoply_normalize_ms_data](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_normalize_ms_data) module.

The report provides:

* Table summarizing the number of features, number of samples and normalization method applied.

* Interactive box-whisker plots visualizing the data distributions **before** normalization.

* Interactive box-whisker plots visualizing the data distributions **after** normalization.

* Interactive line plot depicting normalization coefficients (center, scale) across samples.

## Input

### Required inputs:

* ```tarball```: (`.tar` file) tarball from ```panoply_normalize_ms_data```
* ```label```: (String) label used in the report title
* ```type```: (String) data type
* ```tmpDir```: (String, default="tmp") name of directory used to extract the tarball
* ```yaml```: (File) updated yaml parameter file output from ```panoply_normalize_ms_data``` which indicates if normalization was performed

## Output

* ```report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report. If no normalization has been applied (e.g. the data was already normalized), a placeholder file with a respective message will be returned.
