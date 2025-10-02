# ```panoply_blacksheep_workflow```

## Description

This workflow performs [blacksheep outlier analysis](./Data-Analysis-Modules%3A-panoply_blacksheep) and generates an [interactive report](./Report-Modules%3A-panoply_blacksheep_report.md). It executes the following modules:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_blacksheep</code>](./Data-Analysis-Modules%3A-panoply_blacksheep)         |  performs blacksheep analysis  |
| [<code>panoply_blacksheep_report</code>](./Report-Modules%3A-panoply_blacksheep_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the blacksheep results |


## Input

### Required inputs:

* ```input_gct```: (`.gct` file) normalized data table in the form of a GCT. Any ome type is accepted, will be aggregated to gene level.
* ```master_yaml```: (`.yaml` file) master parameters file
* ```output_prefix```: (String) prefix for naming the output tar file
* ```type```: (String) data (-ome) type

### Optional inputs:

#### panoply_blacksheep

* ```groups_file```: (`.csv` file **RECOMENDED**, default = NULL) subset of sample annotations to be used for calculating enrichment of outliers in these groups. If no groups file is provided, enrichment analysis will not be performed and only outlier counts will be output. **Groups file is recommended** when running [panoply_blacksheep_report](./Report-Modules%3A-panoply_blacksheep_report.md).

* ```apply_filtering```: (Boolean, default = FALSE) flag for filtering the data table. If FALSE, no filtering will occur and outlier results for all gene symbols will be output. If TRUE, data table will be filtered according to ```identifiers_file```.
* ```identifiers_file```: (`.txt` file, default = NULL) used in combination with ```apply_filtering``` = TRUE, this is an option to supply a list of gene symbols for filtering the data table. If no file is provided, a list of kinases included in the module will be used for filtering. If a text file is supplied, it will be used for filtering. Gene symbols should be separated by new lines in the user-provided text file.
* ```fraction_samples_cutoff```: (Float, default = 0.3) value for minimum fraction of samples containing outliers in a given group.
* ```fdr_value```: (Float, default = 0.05) FDR value cutoff to be considered significant

## Output

`panoply_blacksheep_workflow` produces the follow outputs:

* ```blacksheep_tar```: (`.tar`) An output tar file that contains the outlier analysis results in a directory called "blacksheep." Subdirectories with the results for "positive" and "negative" outlier analysis include:
  - Outlier count table (`*_outliers_table.csv`): result of initial outlier calculation at the input feature level. Binary table indicating whether or not a given sample contains an outlier for a particular feature; 0 = not outlier, 1 = outlier. **Note**: This is the only output if no `groups_file` is provided.
  - Aggregated outlier count table (`aggregated_*_outliers_table.csv`): outlier count table collapsed to the gene level. Indicates number of outlier features that correspond to a particular gene symbol for a given sample.
  - Fraction of outliers table (`fraction_*_outliers_per_feature.csv`): ratio of features that are outliers to total features for a given gene symbol in a given sample; table of the fraction of features corresponding to a given gene symbol that are outliers.
  - Groups-File Specific Results:
	  - Outlier analysis results (`outlieranalysis_for_*.csv`): results from Fisher's Exact Test for every comparison group given by the groups file. Table includes gene symbols, p-values, fdr values, and breakdown of numbers in the contigency table (ie in-group, out-group, outliers, not outliers).
	  - Heatmaps depicting the fraction of outlier features in significantly (fdr < ```fdr_value```) enriched genes, saved as both pdf and png files (`*.pdf` and `*.png`)
	  - Analysis Log (`outlier_analysis_log.csv`): Log of all attempted comparisons; lists the outlier directionality (`pos_neg`), the comparison annotation-of-interest (`annotation`), the ingroup used for a binary comparison (`ingroup`), the filename prefix associated with that comparison (`binary_annotation`), and whether the comparison produced figures (`has_heatmap`). Used in the report module to retrieve results files by name.

* ```blacksheep_report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.


