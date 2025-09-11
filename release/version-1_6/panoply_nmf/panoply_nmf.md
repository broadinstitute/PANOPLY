# ```panoply_nmf```

## Description

This module performs unsupervised non-negative matrix factorization (NMF)-based clustering, on single or multi-omic data. Before runing analysis, data is preprocessed: multi-omic datasets are combined, data are filtered and normalized, and the dataset is transformed into a non-negative matrix. After running, an appropriate factorization rank is chosen based on the cophenetic correlation and dispersion coefficients.

#### Pre-processing and non-negative transformation of the input data matrix

In order to merge the input-array of data matrices into a single non-negative fully-quantified input-matrix, the following transformations are applied (in order):

* Combine array data-matrices into a single matrix, retaining ome-type as an `rdesc` annotation
* Remove all features which are not fully quantified
* (optional) Apply standard-deviation filter, to remove features with low variance. For multiomics data, the following filtering methods are available:
  - ```global```: apply filter globally to the multi-omics data matrix
  - ```separate```: apply filter to each data type separately
  - ```equal```: filter the multi-omics data matrix such that each data type will be represented by the same number of 
* (optional) Apply z-scoring by  `row` (z-score rows), `col` (z-score columns), or `rowcol` (z-score rows and then columns).
* Transform signed matrix into a non-negative matrix
  - Create a data matrix with all negative numbers zeroed. 
  - Create another data matrix with all positive numbers zeroed and the signs of all negative numbers removed.
  - Concatenate both matrices resulting in a data matrix twice as large as the original, but with positive values only and zeros and hence appropriate for NMF.

### Non-negative matrix factorization (NMF)

Given a factorization rank _k_ (where _k_ is the number of clusters), NMF decomposes a non-negative <code>p x n</code> data matrix _V_ into two matrices _W_ and _H_ such that multiplication of _W_ and _H_ approximates _V_. Matrix _H_ is a <code>k x n</code> matrix whose entries represent weights for each sample (_1_ to _N_) to contribute to each cluster (_1_ to _k_). Matrix _W_ is a <code>p x k</code> matrix representing weights for each feature (_1_ to _p_) to contribute to each cluster (_1_ to _k_). Matrix _H_ is used to assign samples to clusters by choosing the _k_ with maximum score in each column of _H_. Matrix _W_ containing the weights of each feature in a certain cluster is used to derive a list of representative features separating the clusters using the method proposed in (Kim and Park, 2007). Cluster-specific features are further subjected to a 2-sample moderated t-test (Ritchie et al., 2015) comparing the feature abundance between the respective cluster and all other clusters. Derived p-values are adjusted for multiple hypothesis testing using the methods proposed in (Benjamini and Hochberg, 1995).


#### Determination of the factorization rank

To determine the optimal factorization rank _k_ (number of clusters) for the input data matrix, a range of clusters between _k_=`$kmin` and `$kmax` is tested. For each _k_ the matrix _V_ gets factorized  using `$nrun` iterations with random initialization of _W_ and _H_. To determine the optimal factorization rank the pipeline calculates two metrics for each value of  k: 1) cophenetic correlation coefficient measuring how well the intrinsic structure of the data is recapitulated after clustering (`coph`) and 2) the dispersion coefficient of the consensus matrix as defined in (Kim and Park, 2007) measuring the reproducibility of the clustering across `$nrun` iterations (`disp`). The optimal k is defined as the maximum of `disp^(1-coph)` for cluster numbers between k=`$kmin` and `$kmax`.

## Input

### Required inputs:

* ```ome_gcts```: (Array[File]+) Array of GCT files, upon which NMF analysis is run. GCT files should have share common sample IDs; samples which do not appear in all GCTs will be excluded from analysis.
* ```ome_labels```: (Array[String]+) Array of data-labels corresponding to each GCT file (e.g. `ome_gcts: ["proteome-subset.gct", "phosphoproteome-subsest.gct"]`, `ome_labels: ["proteome", "pSTY"]`). Must match the length and order of `ome_gct` exactly.

* ```kmin```: (Int, default = 2) Minimal factorization rank.
* ```kmax```: (Int, default = 8) Maximal factorization rank.
* ```exclude_2```: (Boolean, default = TRUE) If TRUE, 'k=2' will be excluded from calculation of the optimal rank.
* ```nmf_method```: (String, default = "lee"). NMF method supported by the [NMF R-package](https://cran.r-project.org/web/packages/NMF/index.html). Controlled by `method` in `.yaml` file
* ```nrun```: (Int, default = 50) Number of NMF runs with different starting seeds.
* ```seed```: (Int, default = "random") Seed for NMF factorization. To set the seed explicitly, provide a **numeric** value. Providing the string "random" will result in a random seed.

* ```sd_filt_min```: (Float, default = 0.05) Lowest percentile of standard deviation (SD) across row to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest SD will be removed. Will be applied before z-scoring. Controlled by `sd_filt` in `.yaml` file
* ```sd_filt_mode```: (String, default = "global") Determines how the SD filter will be applied to the multi-omics data matrix. Controlled by `filt_mode` in `.yaml` file.
  - `global`: apply filter globally to the multi-omics data matrix
  - `separate`: apply filter to each data type separately
  - `equal`: filter the multi-omics data matrix such that each data type will be represented by the same number of features. The number of features used to represent each data type _N<sub>feat</sub>_ is determined by the data type with smallest number of features. Other data types will be filtered to retain the _N<sub>feat</sub>_ most variable features.
* ```z_score```: (Boolean, default = TRUE) If TRUE, the data matrix will be z-scored according to yaml-exclusive paramter `z_szore_mode`.
* ```z_score_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns).

* ```gene_column```: (String, default = "geneSymbol") (optional) Column name in rdesc in the GCT file that contains gene symbols, used for adding additional feature-annotations. Controlled by global-parameter `gene_id_col` in `.yaml` file.
* ```organism_id```: (String, default = "human") (optional) Organism type, used for gene-mapping if gene_col is provided. Support for \'human\' (Hs), \'mouse\' (Mm), or \'rat\' (Rn). Controlled by global-parameter `organism` in `.yaml` file.
* ```output_prefix```: (String, default = "results_nmf") name of the output ```.tar``` file.
* ```yaml_file```: (`.yaml` file) master-parameters.yaml

 
## Output

* ```results```: (`${label}_NMF_results.tar.gz`) Tar file containing combined pre-processed expression GCTs (signed and non-negative), `.Rdata` file with res.rank output of `nmf()`, and `.Rdata` file with `opt` object containing parameters.
* ```nclust```: Best factorization rank.
* ```preprocess_figs```: Tar file containing preprocessing figures with the results of sd-filtering and z-scoring (if applicable).



## References

1. Kim, H. & Park, H. **Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.** Bioinformatics 23, 1495-1502 (2007).

1. Ritchie, M. E. et al. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** Nucleic Acids Res. 43, e47-e47 (2015).

1. Benjamini, Y. & Hochberg, Y. **Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society.** Series B (Methodological) 57, 289-300 (1995).
