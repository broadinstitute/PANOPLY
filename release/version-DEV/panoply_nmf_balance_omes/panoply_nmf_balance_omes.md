
# ```panoply_nmf_balance_omes```

## Description

This module can optionally be run before the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module, to balance the number of features from input-omes.

## Input

### Required inputs:

* ```label```: (String) name for output tar file
* ```ome_gcts```: (Array[File]+) array of normalized data matrices (e.g. proteome, phosphoproteome, RNA, CNA, etc.) in `.gct` format.
* ```ome_labels```: (Array[String]+) array of labels associated with each gct file (e.g. "prot", "pSTY", "rna', "cna", etc.). Must match the length and order of `ome_gct` exactly.

* ```tol```: (Float, default = 0.01) Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.
* ```var```: (Float, default = 0.9) Explained variance by PCA (between 0-1). Used to extract the number of PCs explaining the specified fraction of variance in the multi-omics data matrix.
* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns). Note that z-scoring can also be performed directly in the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module.


## Output

* ```ome_gcts_balanced```: Array[File]+ array of balance data-matrices in `.gct` format for input into the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module. Mmaintains order of `ome_gcts`.
* ```pdf```: (File) Visualization of the filtering approach to balance the contribution of the data types.   
  