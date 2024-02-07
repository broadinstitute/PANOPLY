# ```panoply_nmf_internal_workflow```

## Description

PANOPLY workflow to perform unsupervised non-negative matrix factorization (NMF)-based clustering, on single or multi-omic data. The workflow expects data matrices (in [GCT format](https://clue.io/connectopedia/gct_format) derived from the same samples (columns) with consistent sample identifers and imported into a Terra workspace by the [PANOPLY startup notebook](https://github.com/broadinstitute/PANOPLY/wiki/PANOPLY-Tutorial#3-run-panoply-startup-notebook)). 


The workflow consists of six modules: 

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_nmf_balance_omes</code>](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_nmf_balance_omes)    |   optional filtering of GCT files, to balance feature-input across multiple omes |
| [<code>panoply_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf)         |   data pre-processing and NMF analysis |
| [<code>panoply_nmf_postprocess</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess)         |   functional chracterization of derived NMF clusters |
| [<code>panoply_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_nmf_report)   |   RMarkown report for `panoply_nmf` and `panoply_nmf_postprocess` |
| [<code>panoply_ssgsea</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea)          |   single sample Gene Set Enrichment Analysis applied to derived clusters  |
| [<code>panoply_ssgsea_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_ssgsea_report)   |   RMarkown report for `panoply_ssgsea` |


### Preprocessing

#### Balancing the contribution of data types

Module(s): [<code>panoply_nmf_balance_omes</code>](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_nmf_balance_omes)

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


#### Pre-processing and non-negative transformation of the input data matrix

Module(s): [<code>panoply_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf)  

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

After NMF the matrix _W_ of feature weights contains two separate weights for positive and negative value (e.g. z-scores) of each feature, respectively. In order to reverse the non-negative transformation and to derive a single signed weight for each feature, each row in matrix _W_ is normalized by dividing by the sum of feature weights in each row. Weights per feature and cluster were then aggregated by keeping the maximal normalized weight and multiplying with the sign of the z-score from the initial data matrix. Thus, the resulting transformed version of matrix _W_<sub>signed</sub> contains signed cluster weights for each feature present in the input matrix.


### Non-negative matrix factorization (NMF)

Module(s): [<code>panoply_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf), [<code>panoply_nmf_postprocess</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess), [<code>panoply_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_nmf_report) 

Given a factorization rank _k_ (where _k_ is the number of clusters), NMF decomposes a non-negative <code>p x n</code> data matrix _V_ into two matrices _W_ and _H_ such that multiplication of _W_ and _H_ approximates _V_. Matrix _H_ is a <code>k x n</code> matrix whose entries represent weights for each sample (_1_ to _N_) to contribute to each cluster (_1_ to _k_). Matrix _W_ is a <code>p x k</code> matrix representing weights for each feature (_1_ to _p_) to contribute to each cluster (_1_ to _k_). Matrix _H_ is used to assign samples to clusters by choosing the _k_ with maximum score in each column of _H_. Matrix _W_ containing the weights of each feature in a certain cluster is used to derive a list of representative features separating the clusters using the method proposed in (Kim and Park, 2007). Cluster-specific features are further subjected to a 2-sample moderated t-test (Ritchie et al., 2015) comparing the feature abundance between the respective cluster and all other clusters. Derived p-values are adjusted for multiple hypothesis testing using the methods proposed in (Benjamini and Hochberg, 1995).


#### Determination of the factorization rank

Module(s): [<code>panoply_nmf</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf), [<code>panoply_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_nmf_report) 

To determine the optimal factorization rank _k_ (number of clusters) for the input data matrix, a range of clusters between _k_=`$kmin` and `$kmax` is tested. For each _k_ the matrix _V_ gets factorized  using `$nrun` iterations with random initialization of _W_ and _H_. To determine the optimal factorization rank the pipeline calculates two metrics for each value of  k: 1) cophenetic correlation coefficient measuring how well the intrinsic structure of the data is recapitulated after clustering (`coph`) and 2) the dispersion coefficient of the consensus matrix as defined in (Kim and Park, 2007) measuring the reproducibility of the clustering across `$nrun` iterations (`disp`). The optimal k is defined as the maximum of `disp^(1-coph)` for cluster numbers between k=`$kmin` and `$kmax`.

#### Cluster membership scores

Module(s): [<code>panoply_nmf_postprocess</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess), [<code>panoply_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_nmf_report) 

For each sample, a cluster membership score is calculated indicating how representative a sample is to each cluster. This score is used to define a set of "core samples" that is most representative for a given cluster, as follows:

For each sample, the difference between its highest cluster membership score and all other cluster membership scores is calculated. If the minimum of these differences exceeds 1/K, where _K_ is the total number of clusters, a sample is considered a core-member.


#### Functional charcterization of NMF clusters

Module(s): [<code>panoply_ssgsea</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_ssgsea), [<code>panoply_ssgsea_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_ssgsea_report), [<code>panoply_nmf_postprocess</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf_postprocess), [<code>panoply_nmf_report</code>](https://github.com/broadinstitute/PANOPLY/wiki/Report-Modules%3A-panoply_nmf_report)

Functional characterization of resulting NMF clusters is performed by projecting the matrix of signed multi-omic feature weights (W<sub>signed</sub>) onto gene sets in `$gene_set_database` via the ```panoply_ssgsea``` module. To derive a single weight for each gene measured across multiple omics data types (e.g. protein, RNA, phosphorylation site, acetylation site) the weight with maximal absolute amplitude is retained.

To test for overrepresentation of categorical variables defined under `group.cols` in `$yaml_file` in the resulting clusters, a Fisher's exact test (R function <code>fisher.test</code>) is used in the set of samples defining the _"cluster core"_ as described above. For continuous variables defined under ```groups.col.continuous``` in ```$yaml_file``` a Wilcoxon rank-sum test (`ggpubr` R-package) used to assess whether the continuous values are differentially distributed between any pair of clusters.


## Input

### Required inputs:

* ```ome_gcts```: (Array[File]+) array of normalized data matrices (e.g. proteome, phosphoproteome, RNA, CNA, etc.) in `.gct` format.
* ```ome_labels```: (Array[String]+) array of labels associated with each gct file (e.g. "prot", "pSTY", "rna', "cna", etc.). Must match the length and order of `ome_gct` exactly.
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```gene_set_database```: (`.gmt` file) gene set database
* ```label```: (String) label

### Optional inputs:

#### panoply_nmf_balance_omes
* ```balance_omes```:  (Boolean, default = TRUE) If TRUE, the contributions of the different data types will be balanced.  
* ```tol```: (Float, default = 0.01) Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.
* ```var```: (Float, default = 0.9) Explained variance by PCA (between 0-1). Used to extract the number of PCs explaining the specified fraction of variance in the multi-omics data matrix.
* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns). Note that z-scoring can also be performed directly in the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module.

#### panoply_nmf
* ```kmin```: (Int, default = 2) Minimal factorization rank.
* ```kmax```: (Int, default = 10) Maximal factorization rank.
* ```exclude_2```: (Boolean, default = TRUE) If TRUE, 'k=2' will be excluded from calculation of the optimal rank.
* ```nmf_method```: (String, default = "lee"). NMF method supported by the [NMF R-package](https://cran.r-project.org/web/packages/NMF/index.html). Controlled by `method` in `.yaml` file
* ```nrun```: (Int, default = 50) Number of NMF runs with different starting seeds.
* ```seed```: (Int, default = "random") Seed for NMF factorization. To set the seed explicitly, provide a **numeric** value. Providing the string "random" will result in a random seed.

* ```sd_filt_min```: (Float, default = 0.05) Lowest percentile of standard deviation (SD) across row to remove from the data. 0 means all data will be used, 0.1 means 10 percent of the data with lowest SD will be removed. Will be applied before z-scoring. Controlled by `sd_filt` in `.yaml` file
* ```sd_filt_mode```: (String, default = "global") Determines how the SD filter will be applied to the multi-omics data matrix. Controlled by `filt_mode` in `.yaml` file.
  - ```global```: apply filter globally to the multi-omics data matrix
  - ```separate```: apply filter to each data type separately
  - ```equal```: filter the multi-omics data matrix such that each data type will be represented by the same number of features. The number of features used to represent each data type _N<sub>feat</sub>_ is determined by the data type with smallest number of features. Other data types will be filtered to retain the _N<sub>feat</sub>_ most variable features.
* ```z_score```: (Boolean, default = TRUE) If TRUE, the data matrix will be z-scored according to yaml-exclusive paramter `z_szore_mode`.
* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns).

* ```gene_column```: (String, default = "geneSymbol") (optional) Column name in rdesc in the GCT file that contains gene symbols, used for adding additional feature-annotations. Controlled by global-parameter `gene_id_col` in `.yaml` file.
* ```organism_id```: (String, default = "human") (optional) Organism type, used for gene-mapping if gene_col is provided. Support for \'human\' (Hs), \'mouse\' (Mm), or \'rat\' (Rn). Controlled by global-parameter `organism` in `.yaml` file.
* ```output_prefix```: (String, default = "results_nmf") name of the output ```.tar``` file.



#### panoply_nmf_postprocess
* ```groups_file```: (`.csv` file, default = NULL) subset of sample annotations, to be used for calculating enrichment of clusters and in the generation of figures. If no groups file is provided, all annotations from the original `cdesc` will be used. Please note that discrepencies between the `cdesc` of different input data matrices may result in unexpected behavior.
* ```feature_fdr```: (Float, default = 0.01) Maximal FDR for feature selection (2-sample t-test).
* ```pval_signif```: (Float, default = 0.01) Maximal p-value for overrepresentation analysis (Fisher's exact test). Controlled by `ora_pval` in `.yaml` file.
* ```max_annot_levels```: (Int, default = 10) Maximal number of levels in an annotation category. Categories with more levels will be excluded from figures and overrepresentation analysis. Controlled by `ora_max_categories` in `.yaml` file.
* ```top_n_features```: (Int, default = 25) Maximal number of driver features, per cluster, to create boxplots / heatmaps for visualizing expression.
 
 
 
## Output

* ```${label}_NMF_results.tar.gz```: Tar file containing combined pre-processed expression GCTs (signed and non-negative), `.Rdata` file with res.rank output of `nmf()`, and `.Rdata` file with `opt` object containing parameters.
* ```nclust```: Best factorization rank.
* ```${label}_NMF_postprocess.tar.gz```: Tar file containing figures and analyses from post-processing.
* ```"${label}_K${nclust}_clusterMembership.tsv"```: TSV file with sample membership scores, consensus mapping, and core-membership.
* ```${label}_nmf_report.html```: Summary report of NMF clustering pipeline.
* ```results_ssgsea.tar```: Results of ```panoply_ssgsea``` applied to the NMF results.
* ```ssGSEA report - ${label}.html```: Summary report of ssGSEA applied to the NMF results.


## References

1. Kim, H. & Park, H. **Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.** Bioinformatics 23, 1495-1502 (2007).

1. Ritchie, M. E. et al. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** Nucleic Acids Res. 43, e47-e47 (2015).

1. Benjamini, Y. & Hochberg, Y. **Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society.** Series B (Methodological) 57, 289-300 (1995).

1. Abdi, H., &#38; Williams, L. J. (2010). **Principal component analysis.** Wiley Interdisciplinary Reviews: Computational Statistics, 2(4), 433-459. https://doi.org/10.1002/wics.101
