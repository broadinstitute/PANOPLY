# ```panoply_normalize_ms_data```

## Description

This module normalizes proteomics (protein/PTM site) data. The normalization methods, especially two-component normalization, assume that the input data are log ratios to a common reference. 

Normalization methods available are:

* Z-scoring (`mean`): mean-centering followed by standard deviation scaling
* Median-MAD normalization (`median`): median-centering followed by median absolute deviation (MAD) scaling
* Two-component mixture model-based normalization (`2comp`): In this method, we assume that for every sample there is a set of unregulated proteins or PTM sites. In the normalized sample, these proteins or PTM sites should have a log ratio centered at zero. In addition, there are proteins or PTM sites that are either up- or downregulated. This normalization scheme attempts to identify the unregulated proteins or PTM sites, and centers the distribution of these log-ratios around zero in order to nullify the effect of differential protein loading and/or systematic MS variation. A 2-component Gaussian mixture model-based is used to achieve this effect. The two Gaussians N (&mu;<sub>1</sub>, &sigma;<sub>1</sub>) and N (&mu;<sub>2</sub>, &sigma;<sub>2</sub>) for a sample *i* are fitted and used in the normalization process as follows: the mode *m*<sub>i</sub> of the log-ratio distribution is determined for each sample using kernel density estimation with a Gaussian kernel and Shafer-Jones bandwidth. A two-component Gaussian mixture model is then fit with the mean of both Gaussians constrained to be *m*<sub>i</sub>, i.e., &mu;<sub>i1</sub> = &mu;<sub>i2</sub> = *m*<sub>i</sub>. The Gaussian with the smaller estimated standard deviation &sigma;<sub>i</sub> = min (&sigma;<sub>i1</sub> , &sigma;<sub>i2</sub>) is assumed to represent the unregulated component of proteins/PTM sites, and is used to normalize the sample by subtracting the mean *m*<sub>i</sub> from each protein/PTM site and dividing by the standard deviation &sigma;<sub>i</sub>. See (Mertins et al., 2016) and (Gillette et al, 2020).

**CAVEAT:** The two-component mixture model-based normalization (`2comp`) method has been tuned for log-transformed ratio (to a common reference) data, and hence cannot be used with intensity-based (e.g., label-free) data.

The module also creates

* Profile plots (density plots) for each sample before and after normalization for comparison.
* Plots of summary of normalizaton statistics.


## Input

Required inputs:

* ```inputData```: (`.tar` file) tarball from ```panoply_parse_sm_table``` or un-normalized input data in `gct` format (when `standalone` is `true`)
* ```type```: (String) proteomics data type
* ```standalone```: (String) set to ```true``` to run as a self-contained module; if ```true``` the ```analysisDir``` input is required
* ```yaml```: (`.yaml` file) parameters in `yaml` format
* ```analysisDir```: (String) name of analysis directory

Optional inputs:

* ```normalizeProteomics```: (String, default chosen in startup notebook) when 'true' normalization will be applied, when 'false' normalization is skipped
* ```normMethod```: (String, default = '2comp') normalization method; options are '2comp', 'median', 'mean'
* ```altMethod```: (String, default = 'median') alternate normalization method for comparison with `normMethod`; downstream modules typically do not generate analyses for the data normalized using `altMethod`
* ```ndigits```: (Int, default = 5) number of decimal digits to use in output tables
* ```outTar```: (String, default = "panoply_normalize_ms_data-output.tar") output `.tar` file name
* ```outTable```: (String, default = "normalized_table-output.gct") output `.gct` normalized file name

## Output

Tarball including the following files in the `normalized-data` subdirectory: 

* Normalized data files:
	* normalized data table (`*-ratio-norm.gct`)
	* normalized data table using alternate normalization method specified in `altMethod` (`*-ratio-*-norm.gct`)

* Plots and normalization statistics
	* profile plot (density of log ratio values for each sample) showing distribution for all samples in input data, before normalization (`*-ratio-profile-plot.pdf`)
	* profile plot showing distribution for samples *after* normalization, using primary normalization method (`*-ratio-norm-profile-plot.pdf`) and alternate normalization method (`*-ratio-median-norm-profile-plot.pdf`)
	* normalization statistics table showing centering and scaling factors for each sample, using primary normalization method (`*-ratio-norm-stats.csv`) and alternate normalization method (`*-ratio-median-norm-stats.csv`)
	* boxplot of normalization statistics for `QC.pass` and `QC.fail` samples (`*-ratio-norm-stats.pdf` and `*-ratio-median-norm-stats.pdf`)

	
## References

* Mertins, P., Mani, D., Ruggles, K., Gillette, M., Clauser, K., Wang, P., Wang, X., Qiao, J., Cao, S., Petralia, F., et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature*  534(7605), 55 - 62. https://dx.doi.org/10.1038/nature18003.
* Gillette, M., Satpathy, S., Cao, S., Dhanasekaran, S., Vasaikar, S., Krug, K., Petralia, F., Li, Y., Liang, W., Reva, B., et al. (2020). Proteogenomic Characterization Reveals Therapeutic Vulnerabilities in Lung Adenocarcinoma. *Cell*  182(1), 200 - 225.e35. https://dx.doi.org/10.1016/j.cell.2020.06.013
