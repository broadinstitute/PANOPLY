# ```panoply_nmf_postprocess```

## Description

This module visualizes and characterizes the clustering results from the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module for the selected factorization-rank `nclust`.


### Non-negative matrix factorization (NMF)

Given a factorization rank _k_ (where _k_ is the number of clusters), NMF decomposes a non-negative <code>p x n</code> data matrix _V_ into two matrices _W_ and _H_ such that multiplication of _W_ and _H_ approximates _V_. Matrix _H_ is a <code>k x n</code> matrix whose entries represent weights for each sample (_1_ to _N_) to contribute to each cluster (_1_ to _k_). Matrix _W_ is a <code>p x k</code> matrix representing weights for each feature (_1_ to _p_) to contribute to each cluster (_1_ to _k_). Matrix _H_ is used to assign samples to clusters by choosing the _k_ with maximum score in each column of _H_. Matrix _W_ containing the weights of each feature in a certain cluster is used to derive a list of representative features separating the clusters using the method proposed in (Kim and Park, 2007). Cluster-specific features are further subjected to a 2-sample moderated t-test (Ritchie et al., 2015) comparing the feature abundance between the respective cluster and all other clusters. Derived p-values are adjusted for multiple hypothesis testing using the methods proposed in (Benjamini and Hochberg, 1995).

### Signed W-Matrix

After NMF the matrix _W_ of feature weights contains two separate weights for positive and negative value (e.g. z-scores) of each feature, respectively. In order to reverse the non-negative transformation and to derive a single signed weight for each feature, each row in matrix _W_ is normalized by dividing by the sum of feature weights in each row. Weights per feature and cluster were then aggregated by keeping the maximal normalized weight and multiplying with the sign of the z-score from the initial data matrix. Thus, the resulting transformed version of matrix _W_<sub>signed</sub> contains signed cluster weights for each feature present in the input matrix.

### Cluster membership scores

For each sample, a cluster membership score is calculated indicating how representative a sample is to each cluster. This score is used to define a set of "core samples" that is most representative for a given cluster, as follows:

For each sample, the difference between its highest cluster membership score and all other cluster membership scores is calculated. If the minimum of these differences exceeds 1/K, where _K_ is the total number of clusters, a sample is considered a core-member.


### Functional charcterization of NMF clusters

To test for overrepresentation of categorical variables defined under `group.cols` in `$yaml_file` in the resulting clusters, a Fisher's exact test (R function <code>fisher.test</code>) is used in the set of samples defining the _"cluster core"_ as described above. For continuous variables defined under ```groups.col.continuous``` in ```$yaml_file``` a Wilcoxon rank-sum test (`ggpubr` R-package) used to assess whether the continuous values are differentially distributed between any pair of clusters.


## Input

### Required inputs:

* ```nmf_results```: (`.tar` File) Results from [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module
* ```nclust```: (Int) factorization rank / number of clusters to generate analyses for. Defaults to 'optimal' nclust if run as part of [panoply_nmf_workflow](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf), but can be manually overriden to visualize any `nclust` between `kmin` and `kmax`.
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```output_prefix```: (String) label

### Optional inputs:

* ```groups_file```: (`.csv` file, default = NULL) subset of sample annotations, to be used for calculating enrichment of clusters and in the generation of figures. If no groups file is provided, all annotations from the original `cdesc` will be used. Please note that discrepencies between the `cdesc` of different input data matrices may result in unexpected behavior.
* ```feature_fdr```: (Float, default = 0.01) Maximal FDR for feature selection (2-sample t-test).
* ```pval_signif```: (Float, default = 0.01) Maximal p-value for overrepresentation analysis (Fisher's exact test). Controlled by `ora_pval` in `.yaml` file.
* ```max_annot_levels```: (Int, default = 10) Maximal number of levels in an annotation category. Categories with more levels will be excluded from figures and overrepresentation analysis. Controlled by `ora_max_categories` in `.yaml` file.
* ```top_n_features```: (Int, default = 25) Maximal number of driver features, per cluster, to create boxplots / heatmaps for visualizing expression.
* ```gene_column```: (String, default = "geneSymbol") (optional) Column name in rdesc in the GCT file that contains gene symbols, used for adding additional feature-annotations. Controlled by global-parameter `gene_id_col` in `.yaml` file.
 
 
## Output

* ```${output_prefix}_NMF_postprocess.tar.gz```: Tar file containing figures and analyses from post-processing.
* ```"${output_prefix}_K${nclust}_clusterMembership.tsv"```: TSV file with sample membership scores, consensus mapping, and core-membership.
* ```${output_prefix}_K${nclust}_W_rowNorm_combined_signed_n*.gct```: GCT file containing signed W-Matrix for GSEA analysis


## References

1. Kim, H. & Park, H. **Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.** Bioinformatics 23, 1495-1502 (2007).

1. Ritchie, M. E. et al. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** Nucleic Acids Res. 43, e47-e47 (2015).

1. Benjamini, Y. & Hochberg, Y. **Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society.** Series B (Methodological) 57, 289-300 (1995).
