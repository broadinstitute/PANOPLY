# ```panoply_filter```

## Description

This module preprocesses and filters proteomics (protein/PTM site) data. It is typically run after the `panoply_normalize_ms_data` module.


Preprocessing (always applied):

* If ```geneIdCol``` with Hugo Gene Symbols missing from sample annotation table, it will be created from the ```proteinIdCol```
* If the sample annotation table contains a `QC.status` column, samples marked `QC.pass` will be retained in the output files. If not, all samples are assumed to be `QC.pass`, and a `QC.status` column is created accordingly.
* If ```separateQCTypes``` is set to 'true', additional output files (e.g. `*-QC.fail.gct`) will be created with non-`QC.pass` samples.


Filters available are:

* If `sdFilterThreshold` is specified, exclude rows with standard deviation less than `sdFilterThreshold`
* If `combineReplicates` replicates is specified and replicates are present in the data (identified by identical values in the `Participant`, `Type` (optional), and `Timepoint` (optional) columns in the sample annotation table), combine values across replicates for each row, using the method specified by `combineReplicates`.
* If `naMax` is specified, exclude rows with more than `naMax` missing values
* If `noNA` is 'true', create an additional table with no missing values


## Input

Required inputs:

* ```inputData```: (`.tar` file) tarball from ```panoply_normalize_ms_data```, or normalized input data in `gct` format (when `standalone` is `true`)
* ```type```: (String) proteomics data type
* ```standalone```: (String) set to ```true``` to run as a self-contained module; if ```true``` the ```analysisDir``` input is required
* ```yaml```: (`.yaml` file) parameters in `yaml` format
* ```analysisDir```: (String) name of analysis directory

Optional inputs:

* ```filterProteomics```: (String, default chosen in startup notebook) when 'true' filtering will be applied, when 'false' filtering is skipped. Preprocessing is always applied, regardless of toggle value.
* ```separateQCTypes```: (String, default = 'false') toggle for generating additional output files, subset to non-`QC.pass` samples (e.g. `*-QC.fail.gct`). Filtering is not applied to these outputs.
* ```geneIdCol```: (String, default = 'geneSymbol') name of (row) annotation column containing gene IDs.
* ```proteinIdCol```: (String, default = 'id') name of (row) annotation column containing protein IDs.
* ```proteinIdType```:  (String, default chosen in startup notebook) keytype of protein IDs in ```proteinIdCol```
* ```combineReplicates```: (String, default = 'mean') method used to combine replicate samples, as are identified by identical values in the `Participant`, `Type` (optional), and `Timepoint` (optional) columns of the sample annotation table. If `null`, replicates will not be combined.
* ```naMax```: (Float, default = 0.7) maximum allowed NA values per row (protein/PTM site); can be fraction between 0-1 or an integer specifying actual number of samples. If `null`, NA values will not be removed.
* ```noMax```: (String, default = 'false') toggle for generating a GCT in which rows (protein/PTM sites) containing any NA values are excluded
* ```sdFilterThreshold```: (Float, default = 0.5) standard deviation (SD) threshold for SD filtering; rows (proteins/PTM sites) with SD less than `sdFilterThreshold` are excluded from the filtered output table. If `null`, sd filtering will not be applied.
* ```ndigits```: (Int, default = 5) number of decimal digits to use in output tables
* ```outTar```: (String, default = "panoply_filter-output.tar") output `.tar` file name
* ```outTable```: (String, default = "filtered_table-output.gct") output `.gct` filtered file name

## Output

Tarball including the following files in the `filtered-data` subdirectory: 

* Filtered data files:
	* data table containing *only* QC-pass samples (`*-ratio-norm.gct`), with no other filters applied
	* filtered data table (`*-ratio-norm-filt.gct`)

* Optional data files:
	* data table containing non-`QC.pass` samples of some {qc.type} (`*-ratio-norm-{qc.type}.gct`), with no other filters applied
	* filtered data table, with rows (protein/PTM sites) containing any NA values excluded (`*-ratio-norm-filt-noNA.gct`)


	
## References

* Mertins, P., Mani, D., Ruggles, K., Gillette, M., Clauser, K., Wang, P., Wang, X., Qiao, J., Cao, S., Petralia, F., et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature*  534(7605), 55 - 62. https://dx.doi.org/10.1038/nature18003.
* Gillette, M., Satpathy, S., Cao, S., Dhanasekaran, S., Vasaikar, S., Krug, K., Petralia, F., Li, Y., Liang, W., Reva, B., et al. (2020). Proteogenomic Characterization Reveals Therapeutic Vulnerabilities in Lung Adenocarcinoma. *Cell*  182(1), 200 - 225.e35. https://dx.doi.org/10.1016/j.cell.2020.06.013
