# ```panoply_normalize_filter_workflow```

## Description

This workflow performs [normalization](./Data-Preparation-Modules%3A-panoply_normalize_ms_data) and [filtering](./Data-Preparation-Modules%3A-panoply_filter) of an input GCT file. It additionally produces an [interactive report](./Report-Modules%3A-panoply_normalize_ms_data_report.md) the normalization results. It executes the following modules:

| module                    | description
| ----------------------- | ---------------------------------------------------------------- |
| [<code>panoply_normalize_ms_data</code>](./Data-Preparation-Modules%3A-panoply_normalize_ms_data)         |  performs data normalization  |
| [<code>panoply_normalize_ms_data_report</code>](./Report-Modules%3A-panoply_normalize_ms_data_report.md)   |  creates an interactive [R Markdown](https://rmarkdown.rstudio.com/) report of the normalization results |
| [<code>panoply_filter</code>](./Data-Preparation-Modules%3A-panoply_filter.md)   |  applies various QC filters to data  |


## Input

### Required inputs:

* ```inputData```: un-normalized, unfiltered input data in `gct` format
* ```ome_type```: (String) proteomics data type (e.g. proteome, phosphoproteome, etc)
* ```job_identifier```: (String) label to be associated with run
* ```yaml```: (`.yaml` file) parameters in `yaml` format

### Optional inputs:

* ```normalizeProteomics```: (String) when 'true' normalization will be applied, when 'false' normalization is skipped
* ```filterProteomics```: (String) when 'true' filtering will be applied, when 'false' filtering is skipped. Preprocessing is always applied, regardless of toggle value.


#### panoply_normalize_filter

* ```ndigits```: (Int, default = 5) number of decimal digits to use in output tables
* ```normMethod```: (String, default = '2comp') normalization method; options are '2comp', 'median', 'mean'
* ```altMethod```: (String, default = 'median') alternate normalization method for comparison with `normMethod`; downstream modules typically do not generate analyses for the data normalized using `altMethod`
* ```outTar```: (String, default = "panoply_normalize_ms_data-output.tar") output `.tar` file name
* ```outTable```: (String, default = "normalized_table-output.gct") output `.gct` normalized file name


#### panoply_normalize_filter_report

* ```separateQCTypes```: (String, default = 'false') toggle for generating additional output files, subset to non-`QC.pass` samples (e.g. `*-QC.fail.gct`). Filtering is not applied to these outputs.
* ```geneIdCol```: (String, default = 'geneSymbol') name of (row) annotation column containing gene IDs.
* ```proteinIdCol```: (String, default = 'id') name of (row) annotation column containing protein IDs.
* ```proteinIdType```:  (String, default chosen in startup notebook) keytype of protein IDs in ```proteinIdCol```
* ```combineReplicates```: (String, default = 'mean') method used to combine replicate samples, as are identified by identical values in the `Participant`, `Type` (optional), and `Timepoint` (optional) columns of the sample annotation table. If `null`, replicates will not be combined.
* ```naMax```: (Float, default = 0.7) maximum allowed NA values per row (protein/PTM site); can be fraction between 0-1 or an integer specifying actual number of samples. If `null`, NA values will not be removed.
* ```noNA```: (String, default = 'false') toggle for generating a GCT in which rows (protein/PTM sites) containing any NA values are excluded
* ```sdFilterThreshold```: (Float, default = 0.5) standard deviation (SD) threshold for SD filtering; rows (proteins/PTM sites) with SD less than `sdFilterThreshold` are excluded from the filtered output table. If `null`, sd filtering will not be applied.
* ```ndigits```: (Int, default = 5) number of decimal digits to use in output tables
* ```outTar```: (String, default = "panoply_filter-output.tar") output `.tar` file name
* ```outTable```: (String, default = "filtered_table-output.gct") output `.gct` filtered file name



## Output

`panoply_normalize_filter_workflow` produces the follow outputs:

* ```filtered_data_table```: (`.gct` file) noramlized, filtered data table 
* ```filtered_tar```: (`.tar`) An output tar file that contains the normalizied, filtered results 
* ```normalize_report```: (`.html` file) Interactive [R Markown](https://rmarkdown.rstudio.com/) report.
* ```output_ome_type```: (String) `ome_type` variable; internal variable used in workflows
