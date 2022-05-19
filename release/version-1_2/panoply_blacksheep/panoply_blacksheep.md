# ```panoply_blacksheep```

## Description

This module utilizes the BlackSheep algorithm for differential extreme value analysis of -omic data (Blumenberg et al., 2019, [preprint](https://www.biorxiv.org/content/10.1101/825067v2.full.pdf)) to identify both positive and negative outliers in the data and subsequently perform enrichment analysis across all groups of interest, if a ```groups_file``` is provided. Outlier count and analysis can optionally be filtered by either a user-supplied list of gene symbols or a list of kinases included in the module. 

Using the [blacksheepr package](https://www.bioconductor.org/packages/release/bioc/html/blacksheepr.html), outliers are first identified within each feature (row) as those values more than 1.5 * IQR above (for positive outliers) or below (for negative outliers) the median across the feature. The workflow stops here if no ```groups_file``` is provided. If a ```groups_file``` is provided, features are aggregated into their parent gene symbol, and both the number and fraction of features that are outliers in each gene are tabulated for each sample. Enrichment analysis using Fisher's Exact Test is then performed to identify genes enriched for outliers within each annotation group provided in the ```groups_file```. To avoid one or a small number of strong samples driving the enrichment analysis, genes are filtered out if their fraction of samples containing outliers in a given group is below a minimum value (default = 0.3). Finally, heatmaps visualizing the fraction of outliers within each significantly enriched gene are produced. 

Note: This module currently does not work on ```groups_file``` annotations that contain NA or empty data values.

## Input

Required inputs:

* ```input_gct```: (`.gct` file) normalized data table in the form of a GCT. Any ome type is accepted, will be aggregated to gene level.
* ```master_yaml```: (`.yaml` file) master parameters file
* ```output_prefix```: (String) prefix for naming the output tar file

Optional inputs:

* ```apply_filtering```: (Boolean, default = FALSE) flag for filtering the data table. If FALSE, no filtering will occur and outlier results for all gene symbols will be output. If TRUE, data table will be filtered according to ```identifiers_file```.
* ```identifiers_file```: (`.txt` file, default = NULL) used in combination with ```apply_filtering``` = TRUE, this is an option to supply a list of gene symbols for filtering the data table. If no file is provided, a list of kinases included in the module will be used for filtering. If a text file is supplied, it will be used for filtering. Gene symbols should be separated by new lines in the user-provided text file.
* ```groups_file```: (`.csv` file, default = NULL) subset of sample annotations to be used for calculating enrichment of outliers in these groups. If no groups file is provided, enrichment analysis will not be performed and only outlier counts will be output. Please note that this module currently does not work on annotations that contain NA or empty data values.
* ```fraction_samples_cutoff```: (Float, default = 0.3) value for minimum fraction of samples containing outliers in a given group.
* ```fdr_value```: (Float, default = 0.05) FDR value cutoff to be considered significant

## Output

An output tar file that contains the outlier analysis results in a directory called "blacksheep." Subdirectories with the results for "positive" and "negative" outlier analysis include:

* Outlier count table (*_outliers_table.csv): result of initial outlier calculation at the input feature level. Binary table indicating whether or not a given sample contains an outlier for a particular feature; 0 = not outlier, 1 = outlier. **Note**: This is the only output if no `groups_file` is provided.
* Aggregated outlier count table (aggregated_*_outliers_table.csv): outlier count table collapsed to the gene level. Indicates number of outlier features that correspond to a particular gene symbol for a given sample.
* Fraction of outliers table (fraction_*_outliers_per_feature.csv): ratio of features that are outliers to total features for a given gene symbol in a given sample; table of the fraction of features corresponding to a given gene symbol that are outliers.
* Outlier analysis results (outlieranalysis_for_*.csv): results from Fisher's Exact Test for every comparison group given by the groups file. Table includes gene symbols, p-values, fdr values, and breakdown of numbers in the contigency table (ie in-group, out-group, outliers, not outliers).
* Heatmaps depicting the fraction of outlier features in significantly (fdr < ```fdr_value```) enriched genes, saved as both pdf and png files (*.pdf and *.png)
* List of genes significantly (fdr < ```fdr_value```) enriched in each comparison group, corresponding to rownames of heatmaps (significant_genes_*.txt)