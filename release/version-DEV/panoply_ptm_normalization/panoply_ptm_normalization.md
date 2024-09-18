# ```panoply_ptm_normalization```

## Description
In a mass-spectrometry proteomics experiment, the abundances of proteins and PTM sites are typically measured as a ratio between a sample (e.g. patient or cell line) and reference (e.g. control) and are performed as separate experiments on the mass spectrometer (Hogrebe et al., 2018). Ratios observed at the PTM level cannot readily be interpreted as change of PTM site stoichiometry alone and the expression of the cognate protein on which the PTM site was detected has to be taken into account.

This module normalizes PTM abundance data to global proteome data using global linear regression (Mani et al., 2021). Specifically, it takes all PTM log-ratios in all samples and regresses them against the log-ratios of cognate proteins. Then, the resulting residuals are the normalized PTM values.


## Input

Required inputs:

* ```proteome_gct```: (`.gct` file) Input GCT file for proteome ratios
* ```ptm_gct```: (`.gct` file) Input GCT file for PTM ratios
* ```yaml```: (.yaml file) parameters in yaml format

**NOTE:** for the input GCTs proteome and PTM, the sample column names must match exactly.

Optional inputs:

* ```accession_number_col```: (String, default = "accession_number") Name of column with protein or PTM accession number
* ```accession_numbers_col```: (String, default = "accession_numbers") Name of column with accession numbers for the protein or PTM **group**

**NOTE:** The inputs described above are for the primary workflow. Additional optional inputs for tasks constituting the workflow are already set to appropriate defaults and do not need to be modified.


## Output

* PTM GCT data table (`<NAME_OF_INPUT_GCT>-proteome-relative-norm.gct.gct`) with normalized PTM abundances with respect to the global proteome abudances

## References

1. Hogrebe, A., von Stechow, L., Bekker-Jensen, D. B., Weinert, B. T., Kelstrup, C. D., &amp; Olsen, J. V. (2018). **Benchmarking common quantification strategies for large-scale phosphoproteomics.** Nature News. https://doi.org/10.1038/s41467-018-03309-6

2. Mani, D.R., Maynard, M., Kothadia, R. et al. **PANOPLY: a cloud-based platform for automated and reproducible proteogenomic data analysis.** Nat Methods 18, 580–582 (2021). https://doi.org/10.1038/s41592-021-01176-6
