# ```panoply_mo_nmf_gct```

## Description

PANOPLY workflow to perform unsupervised non-negative matrix factorization (NMF)-based multi-omics clustering. The workflow expects data matrices (in [GCT format](https://clue.io/connectopedia/gct_format) derived from the same samples (columns) with consistent sample identifers and imported into a Terra workspace by the [PANOPLY startup notebook](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial#3-run-panoply-startup-notebook)).


The workflow consists of four modules: 

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_mo_nmf_pre</code>](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_mo_nmf_pre)    |   filtering of GCT files and creation of a single <code>.tar</code>-file |
| [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf)         |   multi-omics NMF and functional chracterization of derived clusters |
| [<code>panoply_mo_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_mo_nmf_report)   |   RMarkown report for `panoply_mo_nmf` |
| [<code>panoply_ssgsea</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea)          |   single sample Gene Set Enrichment Analysis applied to derived clusters  |
| [<code>panoply_ssgsea_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_ssgsea_report)   |   RMarkown report for `panoply_ssgsea` |


### Preprocessing

#### Balancing the contribution of data types

Module(s): [<code>panoply_mo_nmf_pre</code>](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_mo_nmf_pre)

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

The results of this balancing approach are visualized in the file <code>balance_omes_pdf</code> returned by module <code>panoply_mo_nmf_pre</code>.


#### Non-negative transformation of the input data matrix

Module(s): [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf)  

In order to convert a input matrix of log ratios or z-scores to a non-negative input matrix, the following transformation is applied:

* Create one data matrix with all negative numbers zeroed. 
* Create another data matrix with all positive numbers zeroed and the signs of all negative numbers removed.
* Concatenate both matrices resulting in a data matrix twice as large as the original, but with positive values only and zeros and hence appropriate for NMF.

After NMF the matrix _W_ of feature weights contains two separate weights for positive and negative value (e.g. z-scores) of each feature, respectively. In order to reverse the non-negative transformation and to derive a single signed weight for each feature, each row in matrix _W_ is normalized by dividing by the sum of feature weights in each row. Weights per feature and cluster were then aggregated by keeping the maximal normalized weight and multiplying with the sign of the z-score from the initial data matrix. Thus, the resulting transformed version of matrix _W_<sub>signed</sub> contains signed cluster weights for each feature present in the input matrix.


### Non-negative matrix factorization (NMF)

Module(s): [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf), [<code>panoply_mo_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_mo_nmf_report) 

Given a factorization rank _k_ (where _k_ is the number of clusters), NMF decomposes a non-negative <code>p x n</code> data matrix _V_ into two matrices _W_ and _H_ such that multiplication of _W_ and _H_ approximates _V_. Matrix _H_ is a <code>k x n</code> matrix whose entries represent weights for each sample (_1_ to _N_) to contribute to each cluster (_1_ to _k_). Matrix _W_ is a <code>p x k</code> matrix representing weights for each feature (_1_ to _p_) to contribute to each cluster (_1_ to _k_). Matrix _H_ is used to assign samples to clusters by choosing the _k_ with maximum score in each column of _H_. Matrix _W_ containing the weights of each feature in a certain cluster is used to derive a list of representative features separating the clusters using the method proposed in (Kim and Park, 2007). Cluster-specific features are further subjected to a 2-sample moderated t-test (Ritchie et al., 2015) comparing the feature abundance between the respective cluster and all other clusters. Derived p-values are adjusted for multiple hypothesis testing using the methods proposed in (Benjamini and Hochberg, 1995).


#### Determination of the factorization rank

Module(s): [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf), [<code>panoply_mo_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_mo_nmf_report) 

To determine the optimal factorization rank _k_ (number of clusters) for the multi-omic data matrix, a range of clusters between _k_=`$kmin` and `$kmax` is tested. For each _k_ the matrix _V_ gets factorized  using `$nrun` iterations with random initialization of _W_ and _H_. To determine the optimal factorization rank the pipeline calculates two metrics for each value of  k: 1) cophenetic correlation coefficient measuring how well the intrinsic structure of the data is recapitulated after clustering (`coph`) and 2) the dispersion coefficient of the consensus matrix as defined in (Kim and Park, 2007) measuring the reproducibility of the clustering across `$nrun` iterations (`disp`). The optimal k is defined as the maximum of `disp^(1-coph)` for cluster numbers between k=`$kmin` and `$kmax`.

#### Cluster membership scores

Module(s): [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf), [<code>panoply_mo_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_mo_nmf_report) 

For each sample, a cluster membership score is calculated indicating how representative a sample is to each cluster. Based on this score a set of "core samples" is defined that is most representative for a given cluster. 

Two approaches to define the _"cluster core"_ are currently implemented: 

1) `legacy`-mode: membership score > `$core_membership` (default: 0.5)
2) `mindiff`-mode: minimal membership score difference between all cluster pairs > 1/K, where _K_ is the total number of clusters. 


#### Functional charcterization of NMF clusters

Module(s): [<code>panoply_ssgsea</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea), [<code>panoply_ssgsea_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_ssgsea_report), [<code>panoply_mo_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_mo_nmf), [<code>panoply_mo_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_mo_nmf_report)

Functional characterization of resulting NMF clusters is performed by projecting the matrix of signed multi-omic feature weights (W<sub>signed</sub>) onto gene sets in `$gene_set_database` via the ```panoply_ssgsea``` module. To derive a single weight for each gene measured across multiple omics data types (e.g. protein, RNA, phosphorylation site, acetylation site) the weight with maximal absolute amplitude is retained.

To test for overrepresentation of categorical variables defined under `group.cols` in `$yaml_file` in the resulting clusters, a Fisher's exact test (R function <code>fisher.test</code>) is used in the set of samples defining the _"cluster core"_ as described above. For continuous variables defined under ```groups.col.continuous``` in ```$yaml_file``` a Wilcoxon rank-sum test (`ggpubr` R-package) used to assess whether the continuous values are differentially distributed between any pair of clusters.


## Input

### Required inputs:

* ```cna_ome```: (`cna-aggregate.gct`)  CNA data natrix created by the [PANOPLY startup notebook](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial#3-run-panoply-startup-notebook)
* ```rna_ome```: (`rna-aggregate.gct`) RNA data matrix created by the [PANOPLY startup notebook](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial#3-run-panoply-startup-notebook)
* ```omes```: (Array[File]) One or more data matrices other that CNA and RNA created by the [PANOPLY startup notebook](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial#3-run-panoply-startup-notebook). Currently supported are
     + `proteome-aggregate.gct`
     + `phosphoproteome-aggregate.gct`
     + `acetylproteome-aggregate.gct`
     + `ubiquitylproteome-aggregate.gct`
     + `glycoproteome-aggregate.gct` 
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```gene_set_database```: (`.gmt` file) gene set database
* ```label```: (String) label

### Optional inputs:

#### panoply_mo_nmf_pre
* ```balance_omes```:  (Boolean, default = TRUE) If TRUE, the contributions of the different data types will be balanced.  
* ```tol```: (Float, default = 0.01) Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.
* ```var```: (Float, default = 0.9) Explained variance by PCA (between 0-1). Used to extract the number of PCs explaining the specified fraction of variance in the multi-omics data matrix.
* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns).

#### panoply_mo_nmf
* ```kmin```: (Int, default = 2) Minimal factorization rank.
* ```kmax```: (Int, default = 8) Maximal factorization rank.
* ```nrun```: (Int, default = 50) Number of NMF runs with different starting seeds.
* ```seed```: (String, default = "random") Seed method for NMF factorization.
* ```gene_column```: (String, default = "geneSymbol") Column name in rdesc in the GCT file that contains gene symbols.
* ```z_score```: (Boolean, default = TRUE) If TRUE, the data matrix will be z-scored according to yaml-exclusive paramter `z_szore_mode`.
* ```impute```: (Boolean, default = FALSE) If TRUE, the data matrix will be first filtered for features present in >70% of samples. The remaining missing values will be imputed via KNN (R package `impute`).'
* ```exclude_2```: (Boolean, default = TRUE) If TRUE, 'k=2' will be excluded from calculation of the optimal rank.
* ```sd_min```: (Float, default = 0.05) Lowest percentile of standard deviation (SD) across row to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest SD will be removed. Will be applied before z-scoring (optional).
* ```filt_mode```: (String, default = "global") Determines how the SD filter will be applied to the multi-omics data matrix:
  - ```global```: apply filter globally to the multi-omics data matrix
  - ```separate```: apply filter to each data type separately
  - ```equal```: filter the multi-omics data matrix such that each data type will be represented by the same number of features. The number of features used to represent each data type _N<sub>feat</sub>_ is determined by the data type with smallest number of features. Other data types will be filtered to retain the _N<sub>feat</sub>_ most variable features.
* ```no_plot```: (Boolean, default=FALSE) If TRUE no downstream analyses of the clustering results will be performed and the R-object after NMF clustering will be returned.
* ```output_prefix```: (String, default = "results_nmf") name of the output ```.tar``` file.


### yaml-exclusive parameters:

* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns).
* ```core_membership_mode```: (String, default="mindiff") Method to define core samples.
    + `legacy`  - membership score > `core_membership`
    + `mindiff` - minimal difference between fractional scores of cluster pairs `>1/_K_` 
* ```core_membership```: (Float, default = 0.5) Minimal cluster membership score to define the cluster core. Only used in 'core_membership_mode = "legacy"`.
* ```sample_order```: (String, default = "clust_then_core_then_score") Determines how samples in the heatmaps are ordered.
    + "clust_then_score" - Order by cluster and membership score (use if `$core_membership_mode = "legacy"`).
    + "clust_then_core_then_score" - Order by cluster, cluster core and membership score (use if `$core_membership_mode = "mindiff"`).
* ```method```: (String, default = "lee"). NMF method supported by the [NMF R-package](https://cran.r-project.org/web/packages/NMF/index.html)
* ```impute_k```: (Int, default = 5)
* ```feature_fdr```: (Float, default = 0.01) Maximal FDR for feature selection (2-sample t-test).
* ```ora_pval```: (Float, default = 0.01) Maximal p-value (Bonferroni-corrected) for overrepresentation analysis (Fisher's exact test).
* ```ora_max_categories```: (Int, default = 10) Maximal number of levels in an annotation category. Categories with more levels will be excluded from the overrepresentation analysis.
 
 
 
## Output

* ```results_nmf.tar```: Results of NMF clustering pipeline.
* ```NMF clustering results - ${label}.html```: Summary report of NMF clustering pipeline.
* ```results_ssgsea.tar```: Results of ```panoply_ssgsea``` applied to the NMF results.
* ```ssGSEA report - ${label}.html```: Summary report of ssGSEA applied to the NMF results.


## References

1. Kim, H. & Park, H. **Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.** Bioinformatics 23, 1495-1502 (2007).

1. Ritchie, M. E. et al. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** Nucleic Acids Res. 43, e47-e47 (2015).

1. Benjamini, Y. & Hochberg, Y. **Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society.** Series B (Methodological) 57, 289-300 (1995).

1. Abdi, H., &#38; Williams, L. J. (2010). **Principal component analysis.** Wiley Interdisciplinary Reviews: Computational Statistics, 2(4), 433-459. https://doi.org/10.1002/wics.101
