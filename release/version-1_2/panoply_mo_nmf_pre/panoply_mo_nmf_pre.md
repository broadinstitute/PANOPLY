
# ```panoply_mo_nmf_pre```

## Description

This module assembles the input for the [panoply_mo_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf) module.

## Input

### Required inputs:

* ```label```: (String) name for output tar file
* ```omes```: (Array[File]) array of normalized data matrices (proteome, phosphoproteome, ubiqitylome, and/or acetylome) in `.gct` format.
* ```rna_ome```: (File) rna data matrix `.gct` file
* ```cna_ome```: (File) cna data matrix `.gct` file

* ```balance_omes```:  (Boolean, default = TRUE) If TRUE, the contributions of the different data types will be balanced.  
* ```tol```: (Float, default = 0.01) Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.
* ```var```: (Float, default = 0.9) Explained variance by PCA (between 0-1). Used to extract the number of PCs explaining the specified fraction of variance in the multi-omics data matrix.
* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns).


## Output

* ```tar```: (File) Tarball including:

  - `.gct` files defined in input parameter `omes`
  - `nmf.conf` - mapping `.gct` files to omic data type

* ```pdf```: (File) Visualization of the filtering approach to balance the contribution of the data types.   
  