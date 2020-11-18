Documentation at https://github.com/broadinstitute/PANOPLY/blob/version-1_0/release/version-1_0/panoply_association/panoply_association.md

# ```panoply_association```

## Description
Performs **association analysis** to identify differential markers for classes of interest using a moderated t-test (for binary classes) or F-test (for categorical multi-level classes). The significant markers are then ranked using a combination of p-values and variable importance in accurate classifiers. 

Classes used for the analysis are derived from the annotations provided in the ```groups``` file, as long as each level in a class has at least 3 samples.

This module performs the following steps:

* Filters out (removes) features with >50% missing values, and imputes missing values using [k-NN](https://www.rdocumentation.org/packages/bnstruct/versions/1.0.2/topics/knn.impute) imputation.
* Runs marker selection using [LIMMA](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).
* Uses statistically significant markers to build classifiers using [PLS](https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf), [PAM](http://statweb.stanford.edu/~tibs/PAM/pam.pdf), [GLMNET](https://www.rdocumentation.org/packages/glmnet/versions/2.0-16/topics/glmnet), and 
[RF](https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest) models.
* Extracts marker importance (a measure of how crucial the marker was for training an accurate classifier) for input markers and combines all the ranks into a global ranking using [rank aggregation](https://cran.r-project.org/web/packages/RobustRankAggreg/RobustRankAggreg.pdf). 
* If provided with test data, runs prediction using all classifiers.
* Plots a heatmap of training data showing significant markers.
* Runs [GSEA](http://software.broadinstitute.org/gsea/index.jsp) for (binary) classes with only two levels.

## Input

Required inputs:

* ```inputData```: (`.tar` file) tarball from ```panoply_parse_sm_table``` or other PANOPLY module;\
(`.gct` file) normalized/filtered input if ```standalone``` is ```TRUE```
* ```type```: (String) (proteome) data type
* ```standalone```: (String) set to ```TRUE``` to run as a self-contained module;\
if ```TRUE``` the ```analysisDir``` and ```groupsFile``` inputs are required
* ```yaml```: (`.yaml` file) parameters in `yaml` format

Optional inputs:

* ```analysisDir```: (String) name of analysis directory
* ```groupsFile```: (`.csv` file) subset of sample annotations, providing classes for association analysis
* ```fdr_assoc```: (Float, default = 0.01) FDR cutoff value for significance
* ```sample_na_max```: (Float, default = 0.8) maximum allowed fraction of NA values per sample/column; error if violated.
* ```nmiss_factor```: (Float, default = 0.5) features (genes, proteins, PTM sites) with more than ```nmiss_factor``` fraction of NA values will be removed from the analysis
* ```duplicate_gene_policy```: (String, default = 'maxvar') method used to combine duplicate genes (when mapping protein accession or PTM site to gene symbols) for running GSEA; possible options are:
  - maxvar: select row with largest variance
  - union: union of binary (0/1) values in all rows (e.g for mutation status)
  - median: median of values in all rows (for each column/sample)
  - mean: mean of values in all rows (for each column/sample)
  - min: minimum of values in all rows (for each column/sample)
* ```gene_id_col```: (String, default = 'geneSymbol') name of sample annotation column containing gene ids.
* ```outFile```: (String, default = "panoply_association-output.tar") output `.tar` file name

## Output

Tarball of files containing the following in the `association` subdirectory, for *each class vector* considered for association analysis: 

* List of significant differential markers derived using LIMMA (`*-markers-all-fdr*.csv`) and p-values for all input features (`*-markers-all*.csv`)
* Marker importance for significant markers, along with final rank (`*-markerimp-fdr*.csv`)
* Heatmap of significant differential markers (`*-markers-heatmap.pdf`)
* Classifier performance contingency tables (`*-analysis-model-results.txt`)
* Table of prediction results for training data (`*-train-results-*.csv`) and testing data (`*-test-results.csv`) using all classifiers
* GSEA outputs, along with `.gct`. and `.cls` input files, for each binary class (in `*-gsea-analysis/` subdirectory).
