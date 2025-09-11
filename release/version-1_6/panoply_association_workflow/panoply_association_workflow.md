# ```panoply_association_workflow```

## Description

This workflow performs [association](./Data-Analysis-Modules%3A-panoply_association) to identify differential markers for classes of interest, and create an [interactive report](./Report-Modules%3A-panoply_association_report.md). The workflow executes the following modules:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_association</code>](./Data-Analysis-Modules%3A-panoply_association)         |   performs association analysis, to identify marker genes associated with provided annotations |
| [<code>panoply_accumulate</code>](./Support-Modules%3A-panoply_accumulate.md)         |   assembles results from association analysis to run [panoply_ssgsea](./Data-Analysis-Modules%3A-panoply_ssgsea) |
| [<code>panoply_ssgsea</code>](./Data-Analysis-Modules%3A-panoply_ssgsea)         |   Performs ssGSEA on association contrast values  |
| [<code>panoply_association_report</code>](./Report-Modules%3A-panoply_association_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the association results |



### Association Analysis

Module(s): [<code>panoply_association</code>](./Data-Analysis-Modules%3A-panoply_association)

This module performs **association analysis** to identify differential markers for classes of interest using a moderated t-test (for binary classes) or F-test (for categorical multi-level classes). The significant markers are then ranked using a combination of p-values and variable importance in accurate classifiers. 

Classes used for the analysis are derived from the annotations provided in the ```groups``` file, as long as each level in a class has at least 3 samples.


### ssGSEA Analysis

Module(s): [<code>panoply_accumulate</code>](./Support-Modules%3A-panoply_accumulate.md), [<code>panoply_ssgsea</code>](./Data-Analysis-Modules%3A-panoply_ssgsea)

ssGSEA is performed on the contrast values (i.e. coefficients of a limma linear model) for all marker features. This analysis is scattered across the results from every class of interest.


### Interactive Report

Module(s): [<code>panoply_association_report</code>](./Report-Modules%3A-panoply_association_report.md)

This module creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the GSEA results of the [panoply_association](https://github.com/broadinstitute/PANOPLY/wiki/Data-Analysis-Modules%3A-panoply_association) module. For each class, for each contrast, an interactive volcano plot shows the Normalized Enrichment Score (NES) vs -log10 of the FDR value for each database pathway analyzed.


## Input

### Required inputs:

* ```inputData```: (`.tar` file) tarball from ```panoply_parse_sm_table``` or other PANOPLY module;\
(`.gct` file) normalized/filtered input if ```standalone``` is ```TRUE```
* ```association_groups``` (`.csv` file) subset of sample annotations, providing classes for association analysis

* ```job_identifier```: (String) label for job
* ```type```: (String) (proteome) data type
* ```standalone```: (String) set to ```TRUE``` to run as a self-contained module;\
**if ```TRUE``` the ```analysisDir``` and ```groupsFile``` inputs are required**
* ```yaml```: (`.yaml` file) parameters in `yaml` format

* ```geneset_db```: (`.gmt` file) gene set database


### Optional inputs:

#### panoply_association

* ```sample_na_max```: (Float, default = 0.8) maximum allowed fraction of NA values per sample/column; error if violated.
* ```nmiss_factor```: (Float, default = 0.5) features (genes, proteins, PTM sites) with more than ```nmiss_factor``` fraction of NA values will be removed from the analysis
* ```duplicate_gene_policy```: (String, default = 'maxvar') method used to combine duplicate genes (when mapping protein accession or PTM site to gene symbols) for running GSEA; possible options are:
  - maxvar: select row with largest variance
  - union: union of binary (0/1) values in all rows (e.g for mutation status)
  - median: median of values in all rows (for each column/sample)
  - mean: mean of values in all rows (for each column/sample)
  - min: minimum of values in all rows (for each column/sample)
* ```gene_id_col```: (String, default = 'geneSymbol') name of sample annotation column containing gene ids.


#### panoply_association_report

* ```fdr_value```: (Float, default = 0.01) FDR value cutoff to be considered significant.


## Output

* ```outputs```: Tarball of files containing the following in the `association` subdirectory, for *each class vector* considered for association analysis: 
  - List of significant differential markers derived using LIMMA (`*-markers-all-fdr*.csv`) and p-values for all input features (`*-markers-all*.csv`)
  - Marker importance for significant markers, along with final rank (`*-markerimp-fdr*.csv`)
  - Heatmap of significant differential markers (`*-markers-heatmap.pdf`)
  - Classifier performance contingency tables (`*-analysis-model-results.txt`)
  - Table of prediction results for training data (`*-train-results-*.csv`) and testing data (`*-test-results.csv`) using all classifiers
  - GSEA outputs, along with `.gct`. and `.cls` input files, for binary classes (in `*-gsea-analysis/` subdirectory).
* ```contrasts```: Tarfile of `.gct` files containing the contrast values, for each class of interest.
* ```ssgsea_assoc_tars```: Array of tarfiles containing the results of ssGSEA analysis on each of the contrast tarfiles.
* ```report```: Report summarizing the ssGSEA analysis on association results.


