# ```panoply_mo_nmf```

## Description

Workflow to perform unsupervised non-negative matrix factorization (NMF)-based multi-omics clustering. 


### Non-negative matrix factorization (NMF)

Given a factorization rank k (where k is the number of clusters), NMF decomposes a non-negative p x n data matrix V into two matrices W and H such that multiplication of W and H approximates V. Matrix H is a k x n matrix whose entries represent weights for each sample (1 to N) to contribute to each cluster (1 to k), whereas matrix W is a p x k matrix representing weights for each feature (1 to p) to contribute to each cluster (1 to k). Matrix H is used to assign samples to clusters by choosing the k with maximum score in each column of H. For each sample, a cluster membership score is calculated as the maximal fractional score of the corresponding column in matrix H. A "cluster core" is defined as the set of samples with cluster membership score > 0.5. Matrix W containing the weights of each feature in a certain cluster was used to derive a list of representative features separating the clusters using the method proposed in (Kim and Park, 2007). Cluster-specific features are further subjected to a 2-sample moderated t-test (Ritchie et al., 2015) comparing the feature abundance between the respective cluster and all other clusters. Derived p-values are adjusted for multiple hypothesis testing using the methods proposed in (Benjamini and Hochberg, 1995).

### Non-negative transformation of the input data matrix

In order to convert a input matrix of log ratios or z-scores to a non-negative input matrix, we applied the followings transformation:

* Create one data matrix with all negative numbers zeroed. 
* Create another data matrix with all positive numbers zeroed and the signs of all negative numbers removed.
* Concatenate both matrices resulting in a data matrix twice as large as the original, but with positive values only and zeros and hence appropriate for NMF.

After NMF the matrix W of feature weights contains two separate weights for positive and negative value (e.g. z-scores) of each feature, respectively. In order to reverse the non-negative transformation and to derive a single signed weight for each feature, each row in matrix W is first normalized by dividing by the sum of feature weights in each row. Weights per feature and cluster were then aggregated by keeping the maximal normalized weight and multiplying with the sign of the z-score from the initial data matrix. Thus, the resulting transformed version of matrix W<sub>signed</sub> contains signed cluster weights for each feature present in the input matrix.


### Determination of the factorization rank

To determine the optimal factorization rank k (number of clusters) for the multi-omic data matrix, a range of clusters between k=`$kmin` and `$kmax` is tested. For each k the matrix V gets factorized  using `$nrun` iterations with random initializations of W and H. To determine the optimal factorization rank the pipeline calculates two metrics for each value of  k: 1) cophenetic correlation coefficient measuring how well the intrinsic structure of the data is recapitulated after clustering and 2) the dispersion coefficient of the consensus matrix as defined in (Kim and Park, 2007) measuring the reproducibility of the clustering across `$nrun` iterations. The optimal k is defined as the maximum of the product of both metrics for cluster numbers between k=`$kmin` and `$kmax`.

### Functional charcterization of NMF clusters

Functional characterization of resulting NMF clusters is performed by projecting the matrix of signed multi-omic feature weights (W<sub>signed</sub>) onto gene sets in `$gene_set_database` via the ```panoply_ssgsea``` module. To derive a single weight for each gene measured across multiple omics data types (e.g. protein, RNA, phosphorylation site, acetylation site) the weight with maximal absolute amplitude is retained.

To test for overrepresentation of categorial variables defined under `group.cols` in `$yaml_file` in the resulting clusters, a Fisher's exact test (R function fisher.test) was used in the set of samples defining the cluster core as described above. For continous variables defined under ```groups.col.continuous``` in ```$yaml_file``` a Wilcoxon rank-sum test (ggpubr R-package) used to assess whether the continuous values are differentially distributed between any pair of clusters.


## Input

### Required inputs:

* ```tar_file```: (`.tar` file) Contains all data tables in GCT v1.3 format as well as a ```.tsv``` file ```nmf.conf``` with no header line and two columns: first column contains the data type (e.g. prot, pSTY, CNV) whereas the seconf column contaons the filenae of the respective GCT file.
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```gene_set_database```: (`.gmt` file) gene set database
* ```label```: (String) label

### Optional inputs:

#### panoply_mo_nmf
* ```kmin```: (Int, default = 2) Minimal factorization rank.
* ```kmax```: (Int, default = 8) Maximal factorization rank.
* ```nrun```: (Int, default = 50) Number of NMF runs with different starting seeds.
* ```seed```: (String, default = 'random') Seed method for NMF factorization.
* ```gene_column```: (String, default = 'geneSymbol') Column name in rdesc in the GCT file that contains gene symbols.
* ```z_score```: (Boolean, default = TRUE) If TRUE, rows in the data matrix will be z-scored.
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
* ```core_membership```: (Float, default=0.5) Minimal cluster membership score to define the cluster core. 
* ```method```: (String, default=brunet). NMF method
* ```impute_k```: (Int, default=5)

## Output

* ```results_nmf.tar```: Results of NMF clustering pipeline.
* ```NMF clustering results - ${label}.html```: Summary report of NMF clustering pipeline.
* ```results_ssgsea.tar```: Results of ```panoply_ssgsea``` applied to the NMF results.
* ```ssGSEA report - ${label}.html```: Summary report of ssGSEA applied to the NMF results.


## References

1. Kim, H. & Park, H. **Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.** Bioinformatics 23, 1495-1502 (2007).

1. Ritchie, M. E. et al. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** Nucleic Acids Res. 43, e47-e47 (2015).

1. Benjamini, Y. & Hochberg, Y. **Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society.** Series B (Methodological) 57, 289-300 (1995).
