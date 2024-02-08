# ```panoply_nmf_workflow```

## Description
Performs unsupervised non-negative matrix factorization (NMF)-based clustering on multi-omic data. Can be used to run the [panoply_nmf_internal_workflow](./Data-Analysis-Modules%3A-panoply_nmf_internal_workflow), to perform multi-omic NMF on all input data matrices combined and/or single-omic NMF on each data matrix independently. Multi-omic and Single-omic NMF results are compared 

This pipeline executes the following workflows:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_nmf_internal_workflow</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf)         |   performs NMF analysis, post-processing, and characterization |
| [<code>panoply_sankey_workflow</code>](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_sankey_workflow)         |   generates sankey-diagram comparisons of SO-NMF results, requires `run_so_nmf=TRUE` |
| [<code>panoply_nmf_assemble_results</code>](https://github.com/broadinstitute/PANOPLY/wiki/Support-Modules%3A-panoply_nmf_assemble_results)   |   compiles all results and reports into single tar |


## Input

### Required inputs:

* ```ome_pairs```: (Array[Pair[String?,File?]]+) Input data matrices for NMF clustering, formatted as an array of pairs. Each pair must contain an identifying label (left) and a data matrix file in `.gct` format (right), e.g. `[ { "left":"prot", "right":this.proteome_ss }, ... ]`.
* ```yaml_file```: (`.yaml` file) master-parameters.yaml
* ```gene_set_database```: (`.gmt` file) gene set database
* ```label```: (String) label
* ```run_mo_nmf```:  (Boolean) Toggle for running multi-omic NMF on all data matrices combined
* ```run_so_nmf```:  (Boolean) Toggle for running single-omic NMF on each data matrix independently
* ```run_sankey```:  (Boolean) Toggle for comparing single-omic NMF results via sankey-diagrams. If multi-omic NMF has been run, its results will also be included in comparisons.

### Optional inputs:

#### panoply_nmf_balance_omes
* ```balance_omes```:  (Boolean, default = TRUE) If TRUE, the contributions of the different data types will be balanced.  
* ```tol```: (Float, default = 0.01) Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.
* ```var```: (Float, default = 0.9) Explained variance by PCA (between 0-1). Used to extract the number of PCs explaining the specified fraction of variance in the multi-omics data matrix.
* ```zscore_mode```: (String, default = "rowcol") z-score mode: `row` (z-score rows), `col` (z-score columns), `rowcol` (z-score rows and then columns). Note that z-scoring can also be performed directly in the [panoply_nmf](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_nmf) module.

#### panoply_nmf (preprocessing)

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

#### panoply_nmf
* ```kmin```: (Int, default = 2) Minimal factorization rank.
* ```kmax```: (Int, default = 10) Maximal factorization rank.
* ```exclude_2```: (Boolean, default = TRUE) If TRUE, 'k=2' will be excluded from calculation of the optimal rank.
* ```nmf_method```: (String, default = "lee"). NMF method supported by the [NMF R-package](https://cran.r-project.org/web/packages/NMF/index.html). Controlled by `method` in `.yaml` file
* ```nrun```: (Int, default = 50) Number of NMF runs with different starting seeds.
* ```seed```: (Int, default = "random") Seed for NMF factorization. To set the seed explicitly, provide a **numeric** value. Providing the string "random" will result in a random seed.

#### panoply_nmf_postprocess
* ```groups_file```: (`.csv` file, default = NULL) subset of sample annotations, to be used for calculating enrichment of clusters and in the generation of figures. If no groups file is provided, all annotations from the original `cdesc` will be used. Please note that discrepencies between the `cdesc` of different input data matrices may result in unexpected behavior.
* ```feature_fdr```: (Float, default = 0.01) Maximal FDR for feature selection (2-sample t-test).
* ```pval_signif```: (Float, default = 0.01) Maximal p-value for overrepresentation analysis (Fisher's exact test). Controlled by `ora_pval` in `.yaml` file.
* ```max_annot_levels```: (Int, default = 10) Maximal number of levels in an annotation category. Categories with more levels will be excluded from figures and overrepresentation analysis. Controlled by `ora_max_categories` in `.yaml` file.
* ```top_n_features```: (Int, default = 25) Maximal number of driver features, per cluster, to create boxplots / heatmaps for visualizing expression.


## Output

`panoply_nmf_workflow` produces the follow outputs:

* ```nmf_results```: Tar file containing results, post-processing, and GSEA analysis from all NMF runs (single-omic and multi-omic). Also contains sankey diagrams (if run). Compiled by [panoply_nmf_assemble_results](./Support-Modules%3A-panoply_nmf_assemble_results).
* ```nmf_results```: Tar file containing reports from all NMF runs (single-omic and multi-omic). Also contains sankey report (if run). Compiled by [panoply_nmf_assemble_results](./Support-Modules%3A-panoply_nmf_assemble_results).



