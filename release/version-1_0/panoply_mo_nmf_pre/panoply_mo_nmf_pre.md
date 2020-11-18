Documentation at https://github.com/broadinstitute/PANOPLY/blob/release-1_0/release/version-1_0/panoply_mo_nmf_pre/panoply_mo_nmf_pre.md


# ```panoply_mo_nmf_pre```

## Description

This module assembles the input for the [panoply_mo_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Analysis-Modules%3A-panoply_mo_nmf) module.

## Input

### Required inputs:

* ```label```: (String) name for output tar file
* ```omes```: (Array[File]) array of normalized ome (proteome, phosphoproteome, ubiqitylome, and/or acetylome) data matrix `.gct` files
* ```rna_ome```: (File) rna data matrix `.gct` file
* ```cna_ome```: (File) cna data matrix `.gct` file

## Output

* ```tar```: (File) Tarball including:

  - `.gct` files defined in input parameter `omes`
  - `nmf.conf` - mapping `.gct` files to omic data type