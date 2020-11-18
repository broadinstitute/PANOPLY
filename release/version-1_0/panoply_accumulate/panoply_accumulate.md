Documentation at https://github.com/broadinstitute/PANOPLY/blob/version-1_0/release/version-1_0/panoply_accumulate/panoply_accumulate.md

# ```panoply_accumulate```

## Description
This module assembles results from association analysis (primarily `panoply_association` and clustering) to run `panoply_ssgsea` on  the fold change(s) for each class vector, to provide additional biological insight into the differential markers for the class.

## Input

Required inputs:
* ```input_tar```: (`.tar` file) tarball from `panoply_association`
* ```module```: (String) panoply analysis module name (usually `panoply_association`)

Optional inputs:
* ```analysisDir```: (String, default = 'input_tarball') name of analysis directory
* ```output_tar```: (String, default = "panoply_contrasts.tar") output `.tar` file name

## Output
Tarball containing results of `module` assembled into appropriate `gct` files for executing `panoply_ssgsea`.

