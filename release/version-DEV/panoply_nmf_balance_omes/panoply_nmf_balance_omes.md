# ```panoply_nmf_balance_omes```

## Description

To mitigate the impact of a potential bias towards a particular data type in the multi-omics clustering (i.e. vastly different number of genomic and proteomic features), the following filtering approach is applied:

* Concatenate data matrices and remove all rows containing missing values
* Standardize the resulting matrix by z-scoring the rows followed by z-scoring of columns 
* Apply principal component analysis (PCA) to the resulting standardized multi-omic data matrix.
  + Based on the _factors matrix_, determine the number of principle components (PCs) explaining 90% of total variance in the data matrix (PCs<sub>90</sub>)
  + Based on the _loadings-matrix_, calculate the relative contribution of each feature to each PCs<sub>90</sub> (equivalent to _squared cosine_ described in (Abdi and Williams, 2010)
  + For each feature calculate relative, cumulative contributions across all PCs<sub>90</sub>
* The resulting vector of relative contributions of each feature (i.e. vector sums up to 1) is then used to balance the contribution of the different data types using the following procedure:
  1. For each data type sum up the contributions of all features; this determines the overall contribution of each data type, which ideally should be equal across the data types within a given tolerance (parameter <code>$tol</code>), i.e.:  
            <code>sum<sub>ome</sub>&asymp;1/(No. data types)</code> 
  2. Remove the feature with the **lowest contribution** that belongs to the data type with the **largest overall contribution**
  3. Recalculate the overall contributions of each data type and repeat steps 1-2 until the deviation is within the specified tolerance (default: <code>tol=0.01</code>). 

The results of this balancing approach are visualized in the file <code>balance_omes_pdf</code> returned by module <code>panoply_nmf_balance_omes</code>.

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
  