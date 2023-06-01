
# ```panoply_so_nmf_pre```

## Description

This module assembles the input for the [panoply_mo_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf) module, with the expectation that only a single-ome is being provided.

## Input

### Required inputs:

* ```label```: (String) name for output tar file

* ```ome```: (File) single-ome data matrix `.gct` file
* ```ome_type```: (String) Type of -omics data provided (**options:** *'RNA', 'CNA', 'prot', 'pSTY', 'acK', 'ubK', 'nglyco'*)


## Output

* ```tar```: (File) Tarball including:

  - `.gct` files defined in input parameter `omes`
  - `nmf.conf` - mapping `.gct` files to omic data type