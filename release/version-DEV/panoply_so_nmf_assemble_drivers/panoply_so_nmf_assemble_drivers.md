
# ```panoply_so_nmf_assemble_drivers```

## Description

This module concatenates the driver features of multiple NMF clusters into a single GCT file. This file can be used as an alternative method of balancing the number of features in each single-ome input to multi-omic clustering. The module takes the tar output of [panoply_mo_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf) as input, with the expectation that only a single-ome was clustered.

## Input

### Required inputs:

* ```nmf_tar```: (`.tar` file) Tar file output of [NMF-based clustering](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf) applied to a single-ome
* ```ome_type```: (String) Type of -omics data (**options:** *'RNA', 'CNA', 'prot', 'pSTY', 'acK', 'ubK', 'nglyco'*)


## Output

* ```driver_pair```: Pair\[File, string\] Pair object including:

  - `ome_type` Type of -omics data
  - `gct` - `.gct` file containing concatenated driver features
  